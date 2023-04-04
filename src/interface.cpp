
#include "misc.hpp"
#include "PT.hpp"
#include "EXI.hpp"
#include "EXII.hpp"
#include "pQCD.hpp"
#include "Therm.hpp"
#include "LatPar.hpp"
#include "ColdQCD.hpp"
#include "pQCD.hpp"
#include "optimize.hpp"
#include "Embedded.hpp"

#include "Standard.hpp"
#include "Continuous.hpp"
// #include "Discontinuous.hpp"
#include "Crossover.hpp"

#include "Scaling.hpp"
#include "RMF.hpp"
#include <pybind11/pybind11.h>



//============================================================================
//Pybind11 Binding

//include each structure in the pybind wrapper. 
//this creates a module "therm" with each model
//as a python class.
PYBIND11_MODULE(therm, m) {

	init_Therm(m);
	init_PT(m);
	init_EXI(m);
	init_EXII(m);
	init_pQCD(m);
	init_Standard(m);
	init_Crossover(m);
	init_Continuous(m);
	// init_Discontinuous(m);
	init_LatPar(m);
	init_ColdQCD(m);
	init_Scaling(m);
	init_Embedded(m);
	init_RMF(m);
	init_optimize(m);
}