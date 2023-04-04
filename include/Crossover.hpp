#ifndef CROSSOVER_H
#define CROSSOVER_H

#include <cmath>
#include <iostream>

#include "Standard.hpp"
#include "Continuous.hpp"
#include "misc.hpp"
#include "Therm.hpp"

#include "PT.hpp"
#include "EXI.hpp"
#include "EXII.hpp"
#include "pQCD.hpp"
#include "RMF.hpp"
#include "LatPar.hpp"

//----------------------------------------------------
//struct for Crossover model (Crossover)

struct Crossover : Therm {

	Crossover()=default;
	Crossover(Therm * pP1, Therm * pP2, Standard * pSF);

	Standard * pSF;
	Therm * pP1;
	Therm * pP2;

	void set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) override;

	void clear_comp_flag() override;

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double *y) override;
	void compute() override;

};

void init_Crossover(py::module_ &m);

#endif