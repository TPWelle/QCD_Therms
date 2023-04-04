#ifndef LatPar_H
#define LatPar_H
#include "Therm.hpp"
#include "misc.hpp"

//----------------------------------------------------
//LatPar model struct and pybind things

struct LatPar : Therm{
    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

    void compute_single(double Temp, double mu_B, double * y) override;
};

//used to include struct in pybind Module
void init_LatPar(pybind11::module_ &m);

#endif