
#ifndef EMBEDDED_H
#define EMBEDDED_H

#include <cmath>
#include <iostream>
#include<vector>

#include "Standard.hpp"
#include "Crossover.hpp"
#include "misc.hpp"
#include "Therm.hpp"

#include "PT.hpp"
#include "EXI.hpp"
#include "EXII.hpp"
#include "pQCD.hpp"
#include "LatPar.hpp"

//----------------------------------------------------
//struct for Embedded CP EOS model

struct Embedded : Therm {

	Embedded()=default;
	Embedded(Therm * pP);
	Embedded(Therm * pP, double Tc, double muc);

	Therm * pP;

	//3D-Ising
	double k; double p;
	double alpha; double beta; double gamma; double delta;

	double Tc; double muc;
	double nc = 0.0; double Pc =-1.0;
	double a0; double Ta;
	double Td;
	double d_p; double d_m;

	int m_exponent=4;

	std::unordered_map<double, double> mux;
	std::unordered_map<double, double> dmux;
	std::unordered_map<double, double> ddmux;


	py::array_t<double> Rs;
	// py::array_t<double> R_T;
	// py::array_t<double> R_m;

	void set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) override;
	void clear_comp_flag() override;

	void set_critical_exponents(double k, double p);

	void mux_insert(double Temp, double mu, double dmu, double ddmu);
	void get_mux(double Temp);
	void get_mux_ellip(double Temp, double mu0);
	void get_muxs();

	double n_BG(double Temp, double mu_B);

	void get_Q(double x,double y, double * Qs);


	double get_R(double Temp);
	void get_R(double Temp, double mu_B, double * Rs);

	void get_rx(double Temp, double mu, double * rxs);
	void get_a(double Temp, double * as);
	void get_Delta(double Temp, double * Delta);

	void compute_single(double Temp, double mu_B, double * y) override;

	void compute() override;

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

};

void init_Embedded(py::module_ &m);


#endif