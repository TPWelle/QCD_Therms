/*
	This file computes the Switching Function(Continuous) for use
	in crossover equation of state.
*/
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Continuous.hpp"
using namespace std;

//----------------------------------------------------
//Continuous model struct and pybind things

Continuous::Continuous(double T0, double Tc, double muc,
 double A, int r_exponent){
	this->T0=T0; this->mu0=3*pi*T0;
	this->Tc=Tc; this->muc=muc;
	this->rho_c=sqrt(pow(Tc/T0,2)+pow(muc/mu0,2));
	this->beta2=pow((Tc*mu0)/(muc*T0),2);
	this->A=A;
	this->r_exponent=r_exponent;

}

Continuous::Continuous(double T0, double mu0, double Tc, double muc,
 double A, int r_exponent){
	this->T0=T0; this->mu0=mu0;
	this->Tc=Tc; this->muc=muc;
	this->rho_c=sqrt(pow(Tc/T0,2)+pow(muc/mu0,2));
	this->beta2=pow((Tc*mu0)/(muc*T0),2);
	this->A=A;
	this->r_exponent=r_exponent;
}

double Continuous::f(double x, double y){
	return .25 + (x*y*atan2(y,x) - y*atan2(y,1.0)  - x*atan2(1.0,x))/pi;
}

double Continuous::get_S(double Temp, double mu_B){
	double t2=(Temp/T0)*(Temp/T0), m2=(mu_B/mu0)*(mu_B/mu0);

	double rho=sqrt(t2+m2);
	double x = (t2 - beta2*m2)/(t2 + beta2*m2);
	double y = (pow(rho,r_exponent) - pow(rho_c,r_exponent))
		/(pow(rho,r_exponent) + pow(rho_c,r_exponent));
	return exp( - pow(rho,-r_exponent) -
		2*A*pow(rho_c,-r_exponent)*f(x,y));
}

//computes a value and derivatives of the switching function
//at a single point Temp,mu_B
void Continuous::compute_single(double Temp, double mu_B, double * y){

	double S=get_S(Temp, mu_B);
	y[0]=S;

	if(divs>0){
		double S_p0,S_0m,S_0p,S_m0;
		double dT=.05, dmu=.1;
		S_p0=get_S(Temp+dT, mu_B);
		S_0m=get_S(Temp, mu_B-dmu);
		S_0p=get_S(Temp, mu_B+dmu);
		S_m0=get_S(Temp-dT, mu_B);
		y[1]=(S_p0-S_m0)/(2*dT);
		if (mu_B==0) y[2]= 0.0;
		else y[2]=(S_0p -S_0m)/(2*dmu);

		if(divs>1) {
			double S_pp,S_mm;
			S_pp=get_S(Temp+dT, mu_B+dmu);
			S_mm=get_S(Temp-dT, mu_B-dmu);
			y[3]=(S_p0-2*S+S_m0)/(dT*dT);
			if (mu_B==0) y[4] = 0.0;
			else y[4] = (S_pp+S_mm+2*S-S_0p-S_0m-S_p0-S_m0)/(2*dT*dmu);
			y[5]=(S_0p-2*S+S_0m)/(dmu*dmu);
		}
	}
}

//Used by pythons pickling and copying features
py::tuple Continuous::getstate(){
	return py::make_tuple(1,Standard::getstate(),r_exponent,Tc,muc);
}

void Continuous::setstate(py::tuple t){
	if (t.size() != 5)
		throw std::runtime_error("Invalid state!");
	Standard::setstate(t[1]);
	r_exponent=t[2].cast<int>();
	Tc=t[3].cast<double>(); muc=t[4].cast<double>();
	rho_c=sqrt(pow(Tc/T0,2)+pow(muc/mu0,2));
	beta2=pow((Tc*mu0)/(muc*T0),2);
}

//used to include in pybind Module
void init_Continuous(py::module_ &m){

		py::class_<Continuous,Standard>(m, "Continuous").def(py::init<>(),"Continuous Object")
		.def(py::init<double, double, double, double, int>(),"Continuous Object")
		.def(py::init<double, double, double, double, double, int>(),"Continuous Object")
		.def("__repr__",[](const Continuous &a) {
			return "<Continuous: T0 = " + std::to_string(a.T0) +
			", mu0 = " + std::to_string(a.mu0) +
			", Tc = " + std::to_string(a.Tc) +
			", muc = " + std::to_string(a.muc) +">";
		}).def(py::pickle(
			[](Continuous &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				Continuous p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			}))
		.def_readwrite("r_exponent",&Continuous::r_exponent)
		.def_readwrite("T0",&Continuous::T0).def_readwrite("mu0",&Continuous::mu0)
		.def_readwrite("Tc",&Continuous::Tc).def_readwrite("muc",&Continuous::muc)
		.def_readwrite("A",&Continuous::A);
}


