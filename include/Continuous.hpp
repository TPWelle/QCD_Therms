#ifndef Continuous_H
#define Continuous_H

#include "Therm.hpp"
#include "Standard.hpp"
#include "misc.hpp"

//----------------------------------------------------
//Continuous model struct and pybind things

struct Continuous : Standard {
	Continuous()=default;
	Continuous(double T0, double Tc, double muc,
 		double A, int r_exponent);

	Continuous(double T0, double mu0, double Tc, double muc,
		double A, int r_exponent);

	int r_exponent;

	double T0; double mu0;
	double Tc; double muc;


	double rho_c;
	double beta2;
	double A;

	//Auxilliary function f(x,y)
	double f(double x, double y);

	//computes switching function for a particulat value of Temp and mu_B
	double get_S(double Temp, double mu_B);

	//Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double * y) override;

};

//used to include in pybind Module
void init_Continuous(py::module_ &m);


#endif