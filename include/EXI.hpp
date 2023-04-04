#ifndef THERMO_EXI_H
#define THERMO_EXI_H
#include "Therm.hpp"


//-----------------------------------------------------------------------------
//Code to switch between Temp, mu_B and T_star, mu_B_star for exI model

double Temp_exI(double T_star, double mu_B_star, double epsilon0);

double mu_B_exI(double T_star, double mu_B_star, double epsilon0);

void Find_Temp_star_mu_B_star_exI(double Temp, double mu_B, double epsilon0,
 double& T_star, double& mu_B_star);


//-----------------------------------------------------------------------------
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_exI(double Temp, double mu_B, double epsilon0);

double Total_Energy_Dens_exI(double Temp, double mu_B, double epsilon0);

// double Total_Number_Dens_exI(double Temp, double mu_B, double epsilon0);

double Total_Baryon_Dens_exI(double Temp, double mu_B, double epsilon0);

double Total_Entropy_Dens_exI(double Temp, double mu_B, double epsilon0);


//----------------------------------------------------
//EXI model struct and pybind things

struct EXI : Therm{

    double eps0=nan("");

    EXI(double eps0);
    EXI()=default;

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

    void compute_single(double Temp, double mu_B, double * y) override;
};

//used to include struct in pybind Module
void init_EXI(py::module_ &m);

#endif
