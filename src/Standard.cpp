/*
	This file computes the Switching Function(Standard) for use
	in crossover equation of state.
*/
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Standard.hpp"
using namespace std;

//----------------------------------------------------
//Standard model struct and pybind things

Standard::Standard(double T0, int r_exponent){
	this->T0=T0; this->mu0=3*pi*T0;
	this->r_exponent=r_exponent;

}

Standard::Standard(double T0, double mu0, int r_exponent){
	this->T0=T0; this->mu0=mu0;
	this->r_exponent=r_exponent;
}

//computes a value and derivatives of the switching function
//at a single point Temp,mu_B
void Standard::compute_single(double Temp, double mu_B, double * y){
	double t=Temp/T0, m=mu_B/mu0;
	int r = r_exponent;

	//Integer determining type of Standard to use. Move this
	// to a model parameter later.
	int Standard_code=0;
	if(Standard_code==0){
		double x=sqrt(t*t + m*m);
		double S=exp(-pow(x,-r));
		y[0]=S;
		if(divs>0){
			y[1]=r*t*pow(x,-(r+2))*S/T0;
			y[2]=r*m*pow(x,-(r+2))*S/mu0;
		} if(divs>1) {
			y[3]=r*S*pow(x,-(r+4))*(r*t*t*pow(x,-r)-r*t*t+ (m*m -t*t))/(T0*T0);
			y[4]=r*S*t*m*pow(x,-(r+4)) *(-(r+2) +r*pow(x,-r))/(T0*mu0);
			y[5]=r*S*pow(x,-(r+4))*(r*m*m*pow(x,-r)-r*m*m+ (t*t -m*m))/(mu0*mu0);
		}
	}else if(Standard_code==1){
		double x=pow(t,r) + pow(m,r);
		double xt=r*pow(t,r-1), xtt=r*(r-1)*pow(t,r-2);
		double xm=r*pow(m,r-1), xmm=r*(r-1)*pow(m,r-2);
		double S=exp(-1.0/x);
		y[0]=S;
		if(divs>0){
			y[1]=xt*pow(x,-2)*S/T0;
			y[2]=xm*pow(x,-2)*S/mu0;
		} if(divs>1) {
			y[3]=S*(xt*xt*(1-2*x) +xtt*x*x)/pow(x,4)/(T0*T0);
			y[4]=S*xt*xm*(1-2*x)/pow(x,4)/(T0*mu0);
			y[5]=S*(xm*xm*(1-2*x) +xmm*x*x)/pow(x,4)/(mu0*mu0);
		}
	}
}

//Used by pythons pickling and copying features
py::tuple Standard::getstate(){
	return py::make_tuple(0,Therm::getstate(),r_exponent,T0,mu0);
}

void Standard::setstate(py::tuple t){
	if (t.size() != 5)
		throw std::runtime_error("Invalid state!");
	Therm::setstate( t[1]);
	r_exponent=t[2].cast<int>();
	T0=t[3].cast<double>(); mu0=t[4].cast<double>();
}

//used to include in pybind Module
void init_Standard(py::module_ &m){

		py::class_<Standard,Therm>(m, "Standard").def(py::init<>(),"Standard Object")
		.def(py::init<double, double>(),"Standard Object")
		.def(py::init<double, double, int>(),"Standard Object")
		.def("__repr__",[](const Standard &a) {
			return "<Standard: T0 = " + std::to_string(a.T0) +
			", mu0 = " + std::to_string(a.mu0) +">";
		}).def(py::pickle(
			[](Standard &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				Standard p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			}))
		.def_readwrite("r_exponent",&Standard::r_exponent)
		.def_readwrite("T0",&Standard::T0).def_readwrite("mu0",&Standard::mu0);
}


