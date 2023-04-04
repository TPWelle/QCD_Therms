#include <iostream>
#include <stdexcept>

#include "Crossover.hpp"
#include <omp.h>
#include <cmath>

using namespace std;

//----------------------------------------------------
//struct for Crossover model (Crossover)

Crossover::Crossover(Therm * pP1, Therm * pP2, Standard * pSF){
	this->pSF=pSF;
	this->pP1=pP1; this->pP2=pP2;
}

void Crossover::set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs){
	Therm::set_Tmu( Temps, mu_Bs);
	clear_comp_flag();
	pSF->divs=divs; pP1->divs=divs; pP2->divs=divs;
	pSF->N_threads=N_threads; pP1->N_threads=N_threads; pP2->N_threads=N_threads;
	pSF->set_Tmu(Temps, mu_Bs); pP1->set_Tmu(Temps, mu_Bs); pP2->set_Tmu(Temps, mu_Bs);
}



void Crossover::clear_comp_flag(){
	this->comp_flag=false;
	pSF->clear_comp_flag();
	pP1->clear_comp_flag();
	pP2->clear_comp_flag();
}


void Crossover::compute_single(double Temp, double mu_B, double * y){

	double yS[6], yP1[6], yP2[6];

	pSF->compute_single(Temp,mu_B,yS);
	pP1->compute_single(Temp,mu_B,yP1); pP2->compute_single(Temp,mu_B,yP2);

	y[0]=(1.0-yS[0])*yP1[0] +yS[0]*yP2[0];
	if(divs>0){
		y[1]=(yP2[0]-yP1[0])*yS[1]+ (1.0-yS[0])*yP1[1] +yS[0]*yP2[1];
		y[2]=(yP2[0]-yP1[0])*yS[2] + (1.0-yS[0])*yP1[2] +yS[0]*yP2[2];
	} if(divs>1) {
		y[3]=(yP2[0]-yP1[0])*yS[3] +2*yS[1]*(yP2[1] - yP1[1])
		+ (1.0-yS[0])*yP1[3] +yS[0]*yP2[3];

		y[4]=(yP2[0]-yP1[0])*yS[4] +yS[2]*(yP2[1] - yP1[1])+
		yS[1]*(yP2[2] -yP1[2]) + (1.0-yS[0])*yP1[4] +yS[0]*yP2[4];

		y[5]=(yP2[0]-yP1[0])*yS[5] +2*yS[2]*(yP2[2] - yP1[2])
		+ (1.0-yS[0])*yP1[5] +yS[0]*yP2[5];
	}
}


void Crossover::compute(){

	if(this->comp_flag) return;

	pSF->compute();
	pP1->compute();	pP2->compute();

	int N=Temps.request().size;

	double *Z= get_arr_pointer( this->Y ), *Z_T= get_arr_pointer( this->Y_T ), *Z_m= get_arr_pointer( this->Y_m );
	double *Z_TT= get_arr_pointer( this->Y_TT ), *Z_Tm= get_arr_pointer( this->Y_Tm ),*Z_mm= get_arr_pointer( this->Y_mm );

	double *S= get_arr_pointer( pSF->Y ), *S_T= get_arr_pointer( pSF->Y_T ), *S_m= get_arr_pointer( pSF->Y_m );
	double *S_TT= get_arr_pointer( pSF->Y_TT ), *S_Tm= get_arr_pointer( pSF->Y_Tm ),*S_mm= get_arr_pointer( pSF->Y_mm );

	double *P1= get_arr_pointer( pP1->Y ), *P1_T= get_arr_pointer( pP1->Y_T ), *P1_m= get_arr_pointer( pP1->Y_m );
	double *P1_TT= get_arr_pointer( pP1->Y_TT ), *P1_Tm= get_arr_pointer( pP1->Y_Tm ),*P1_mm= get_arr_pointer( pP1->Y_mm );

	double *P2= get_arr_pointer( pP2->Y ), *P2_T= get_arr_pointer( pP2->Y_T ), *P2_m= get_arr_pointer( pP2->Y_m );
	double *P2_TT= get_arr_pointer( pP2->Y_TT ), *P2_Tm= get_arr_pointer( pP2->Y_Tm ),*P2_mm= get_arr_pointer( pP2->Y_mm );

	#pragma omp parallel for schedule(dynamic) \
	shared(divs) num_threads(N_threads)
	for(int i=0;i<N;i++){
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();

		if(omp_get_thread_num()%N_threads==0 and display_progress){
			cout<< i<<"/"<<N<<endl;
		}
		Z[i]=(1.0-S[i])*P1[i] +S[i]*P2[i];
		if(divs>0){
			Z_T[i]=(P2[i]-P1[i])*S_T[i] + (1.0-S[i])*(P1_T[i]) +S[i]*(P2_T[i]);
			Z_m[i]=(P2[i]-P1[i])*S_m[i] + (1.0-S[i])*(P1_m[i]) +S[i]*(P2_m[i]);
		} if(divs>1) {
			Z_TT[i]=(P2[i]-P1[i])*S_TT[i] +2*S_T[i]*(P2_T[i] - P1_T[i])
			+ (1.0-S[i])*(P1_TT[i]) +S[i]*(P2_TT[i]);

			Z_Tm[i]=(P2[i]-P1[i])*S_Tm[i] +(S_m[i])*(P2_T[i] - P1_T[i])+
			(S_T[i])*(P2_m[i] - P1_m[i]) + (1.0-S[i])*(P1_Tm[i]) +S[i]*(P2_Tm[i]);

			Z_mm[i]=(P2[i]-P1[i])*S_mm[i] +2*S_m[i]*(P2_m[i] - P1_m[i])
			+ (1.0-S[i])*(P1_mm[i]) +S[i]*(P2_mm[i]);
		}
	}

	this->comp_flag=true;
}



//Used by pythons pickling and copying features
py::tuple Crossover::getstate(){
	return py::make_tuple(10,Therm::getstate(),pSF->getstate()
		,pP1->getstate(),pP2->getstate());
}

void Crossover::setstate(py::tuple t){
	if (t.size() != 5)
		throw std::runtime_error("Invalid state!(Crossover)");

	Therm::setstate(t[1]);
	switch(t[2].cast<py::tuple>()[0].cast<int>()){
		case 0:
			pSF= new Standard; break;
		case 1:
			pSF= new Continuous; break;
	}
	pSF->setstate(t[2]);
	switch(t[3].cast<py::tuple>()[0].cast<int>()){
		case 0:
			pP1= new PT; break;
		case 1:
			pP1= new EXI; break;
		case 2:
			pP1= new EXII; break;
		case 3:
			pP1= new pQCD; break;
		case 4:
			pP1= new LatPar; break;
		case 5:
			pP1= new RMF; break;
		case 10:
			pP1= new Crossover; break;
	}
	pP1->setstate(t[3]);
	switch(t[4].cast<py::tuple>()[0].cast<int>()){
		case 0:
			pP2= new PT; break;
		case 1:
			pP2= new EXI; break;
		case 2:
			pP2= new EXII; break;
		case 3:
			pP2= new pQCD; break;
		case 4:
			pP2= new LatPar; break;
		case 5:
			pP1= new RMF; break;
		case 10:
			pP2= new Crossover; break;
	}
	pP2->setstate(t[4]);
}



//----------------------------------------------------
//pybind things

void init_Crossover(py::module_ &m){
	py::class_<Crossover,Therm>(m, "Crossover").def(py::init<>(),"Crossover Object")
		.def(py::init<Therm *, Therm *, Standard *>(),"Crossover Object")
		.def_readwrite("pSF",&Crossover::pSF)
		.def_readwrite("pP1",&Crossover::pP1).def_readwrite("pP2",&Crossover::pP2)
		.def("set_Tmu", &Crossover::set_Tmu, "Set Ts and mus of contained models")
		.def("__repr__", [](const Crossover &a) {
				return "<Crossover>";
		}).def(py::pickle( [](Crossover &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				Crossover p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		);
}