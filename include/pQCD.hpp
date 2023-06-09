#ifndef THERMO_PQCD_H
#define THERMO_PQCD_H
#include "Therm.hpp"
#include "misc.hpp"

double alphas_3loop(double Lambda, double C_S = 0.0);

double Total_Pressure_PQCD(double Temp, double mu_B,
 double C_E = 1, double C_S = 0.0, int K = 6);

double Total_Pressure_Free_QCD(double Temp, double mu_B);

//---------------------------------------------------
//These aren't necessary. They perform a finite difference
//which can be done just as well in Therm.cpp
//---------------------------------------------------

double Total_Energy_Dens_PQCD(double Temp, double mu_B,
 double C_E = 1, double C_S = 0.0);

double Total_Baryon_Dens_PQCD(double Temp, double mu_B,
 double C_E = 1, double C_S = 0.0);

double Total_Entropy_Dens_PQCD(double Temp, double mu_B,
 double C_E = 1, double C_S = 0.0);


//----------------------------------------------------
//pQCD model struct and pybind things

struct pQCD : Therm{

	double C_S=0.0; double C_E=1.0;
	int K=6; //Max order in g

	pQCD()=default;
	pQCD(double C_E, double C_S);
	pQCD(double C_E, double C_S, int K);

	//Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double * y) override;

};

//used to include struct in pybind Module
void init_pQCD(py::module_ &m);

#endif