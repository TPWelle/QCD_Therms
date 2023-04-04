#ifndef THERMO_EXII_H
#define THERMO_EXII_H
#include "Therm.hpp"

//-----------------------------------------------------------------------------
//Code to compute thermal properties of a point particle gas of many particle
//species

double Total_Pressure_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess = -1.0);

double Total_Energy_Dens_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess = -1.0);

double Total_Baryon_Dens_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess = -1.0);

double Total_Entropy_Dens_exII(double Temp, double mu_B, double epsilon0 ,
 double Press_exII_guess = -1.0);


//----------------------------------------------------
//EXII model struct and pybind things

struct EXII : Therm{
    double eps0=nan("");

    EXII(double eps0);
    EXII()=default;

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

    void compute_single(double Temp, double mu_B, double * y) override;
};

//used to include struct in pybind Module
void init_EXII(py::module_ &m);

#endif