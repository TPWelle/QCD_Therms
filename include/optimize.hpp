#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <vector>
#include <stdexcept>
#include "Therm.hpp"
#include "Crossover.hpp"


std::vector<double> Optimize(Crossover* pth);


//Optimize parameters for crossover eos's
double MinimizeChiSquare(Therm *pth, const std::vector<double *> opti);

void init_optimize(py::module_ &m);

#endif