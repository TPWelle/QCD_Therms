#ifndef THERMO_PT_H
#define THERMO_PT_H

#include <vector>

#include "particle.hpp"
#include "particlelist.hpp"
#include "misc.hpp"
#include "Therm.hpp"

//Fermi/Bose/Classical stats distribution function
double equilib_distfcn(double x, double y, double z, STATISTICS stats);

// double equilib_distfcn(double Energy, double Temp, double mu_particle,
//  QUANTUM_TYPE quantumtype);


double Thermal_Integral(double y, double z, int a, int b, STATISTICS stats);



//-----------------------------------------------------------------------------
//Code to compute (partial) thermal properties of one species of gas:

double Partial_Pressure_pt(double Temp, double mu_particle, 
 const Particle& particle);

double Partial_Number_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

double Partial_Energy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

double Partial_Entropy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

//-----------------------------------------------------------------------------
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_pt(double Temp, double mu_B);

double Total_Energy_Dens_pt(double Temp, double mu_B);

double Total_Baryon_Dens_pt(double Temp, double mu_B);

double Total_Entropy_Dens_pt(double Temp, double mu_B);

//-----------------------------------------------------------------------------
//PT model struct and pybind things

struct PT : Therm{

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

    void compute_single(double Temp, double mu_B, double * y) override;
};


//used to include struct in pybind Module
void init_PT(pybind11::module_ &m);


#endif
