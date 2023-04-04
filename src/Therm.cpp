#include <iostream>
#include <stdexcept>

#include "Therm.hpp"

#include <omp.h>
#include <cmath>

using namespace std;

//Gets pointer to access elements of PyBind arrays
double * get_arr_pointer(py::array_t<double> Y){
	return static_cast<double *>( Y.request().ptr);
} 

//--------------------------------------------------------
//Prepares model for computation by setting the T and mu arrays and those of dependent models.
//Then, allocate arrays for the output of computation
void Therm::set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) {
	if(divs>2 or divs<0)
		throw std::invalid_argument("Invalid divs in Therm.compute");
	py::buffer_info buf = Temps.request();

	if (buf.size != mu_Bs.request().size)
		throw std::runtime_error("Input shapes must match");

	this->Temps=Temps; this->mu_Bs=mu_Bs;
	clear_comp_flag();
	//allocate arrays for output
	Y = py::array_t<double>(buf.size); Y.resize(buf.shape);
	if(divs>0){
		Y_T = py::array_t<double>(buf.size); Y_T.resize(buf.shape);
		Y_m = py::array_t<double>(buf.size); Y_m.resize(buf.shape);
	} if(divs>1) {
		Y_TT = py::array_t<double>(buf.size); Y_TT.resize(buf.shape);
		Y_Tm = py::array_t<double>(buf.size); Y_Tm.resize(buf.shape);
		Y_mm = py::array_t<double>(buf.size); Y_mm.resize(buf.shape);
	}
}


void Therm::clear_comp_flag(){
	//Set comp_flag to false. This means that all values are
	//recomputed the next time .compute() is called.
	this->comp_flag=false;
}



void Therm::compute_single( double Temp, double mu_B,double * y){

	// computes a value and derivatives at a single point Temp,mu_B
	// y is a 6 element array whose elements are in order are:
	// P, s, n, X_TT, X_Tm, X_mm

	for(int i=0;i<6;i++) y[i]=0.0; //does nothing for base class
} 

//compute the values of Y at the the specified values of Temp and mu_B
//also compute appropriate derivatives if divs>0.
void Therm::compute(){

	if(this->comp_flag) return;
	int N=Temps.request().size;

	double *pTemp= get_arr_pointer( Temps ), *pmu= get_arr_pointer( mu_Bs );

	double *pY= get_arr_pointer( Y );
	double *pY_T= get_arr_pointer( Y_T ), *pY_m= get_arr_pointer( Y_m );
	double *pY_TT= get_arr_pointer( Y_TT ), *pY_Tm= get_arr_pointer( Y_Tm ), *pY_mm= get_arr_pointer( Y_mm );

	#pragma omp parallel for schedule(dynamic) \
	shared(divs) num_threads(N_threads)
	for(int i=0;i<N;i++){
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		if(omp_get_thread_num()%N_threads==0 and display_progress){
			cout<< i<<"/"<<N<<endl;
		}
		double Y[6];
		compute_single(pTemp[i],pmu[i],Y);
		pY[i]=Y[0]; pY_T[i]=Y[1]; pY_m[i]=Y[2];
		pY_TT[i]=Y[3]; pY_Tm[i]=Y[4]; pY_mm[i]=Y[5];
	}
	this->comp_flag=true;
}


//Used by pythons pickling and copying features
py::tuple Therm::getstate(){
	return py::make_tuple(divs, comp_flag, display_progress,
		Temps, mu_Bs, Y, Y_T, Y_m, Y_TT, Y_Tm, Y_mm);
}

void Therm::setstate(py::tuple t){
	if (t.size() != 11)
		throw std::runtime_error("Invalid state!");

	divs=t[0].cast<int>();

	comp_flag=t[1].cast<bool>();
	display_progress=t[2].cast<bool>();
	Temps=t[3].cast<py::array_t<double>>();
	mu_Bs=t[4].cast<py::array_t<double>>();
	Y=t[5].cast<py::array_t<double>>();
	Y_T=t[6].cast<py::array_t<double>>();
	Y_m=t[7].cast<py::array_t<double>>();
	Y_TT=t[8].cast<py::array_t<double>>();
	Y_Tm=t[9].cast<py::array_t<double>>();
	Y_mm=t[10].cast<py::array_t<double>>();
}


//used to include struct in pybind Module
void init_Therm(py::module_ &m){
	py::class_<Therm>(m, "Therm").def(py::init<>(),"Therm Object")//, py::arg("div")=-1)
		.def_readwrite("divs", &Therm::divs)
		.def_readwrite("N_threads", &Therm::N_threads)
		.def_readwrite("comp_flag", &Therm::comp_flag)
		.def_readwrite("display_progress", &Therm::display_progress)
		.def_readwrite("Temps",&Therm::Temps).def_readwrite("mu_Bs",&Therm::mu_Bs)
		.def_readwrite("Y",&Therm::Y).def_readwrite("Y_T",&Therm::Y_T).def_readwrite("Y_m",&Therm::Y_m)
		.def_readwrite("Y_TT",&Therm::Y_TT).def_readwrite("Y_Tm",&Therm::Y_Tm).def_readwrite("Y_mm",&Therm::Y_mm)
		.def("set_Tmu", &Therm::set_Tmu, "test")
		.def("compute", &Therm::compute, "Compute specified Ys")
		.def("__repr__", [](const Therm &a) {
				return "<Therm>";
		}).def(py::pickle(
			[](Therm &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				Therm p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		);
}


