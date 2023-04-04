#ifndef COLDQCD_H
#define COLDQCD_H
#include "Therm.hpp"
#include "misc.hpp"


struct ColdQCD : Therm{

	ColdQCD()=default;

	int Nf=3; int Nc=3;

	//Used for python's serialization and copying functions
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double * y) override;
};

void init_ColdQCD(py::module_ &m);




#endif
