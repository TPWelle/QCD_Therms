#ifndef Standard_H
#define Standard_H

#include "Therm.hpp"

//----------------------------------------------------
//Standard model struct and pybind things

struct Standard : Therm{
	Standard()=default;
	Standard(double T0, int r_exponent);
	Standard(double T0, double mu0, int r_exponent);


	int r_exponent;

	double T0; double mu0;

	//Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double * y) override;

};

//used to include in pybind Module
void init_Standard(py::module_ &m);

#endif