#ifndef SCALING_H
#define SCALING_H

#include <cmath>
#include <iostream>

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
//struct for Scaling EOS model

struct Scaling : Therm {

	Scaling()=default;
	Scaling(Therm * pP);
	Scaling(Therm * pP, double Tc, double muc);

	Therm * pP;

	//3D-Ising
	double alpha=.1101; double beta=.3264;
	double gamma=1.2371; double delta=4.7901;

	//Mean-Field
	// double alpha=.0; double beta=.5;
	// double gamma=1.2371; double delta=1.0;

	double Tc; double muc; double nc; double Pc;
	double n0; double P0;

	//Used to select version of EOS
	//0: mux defined by n(T,mux(T)=n_c
	//1: mux=muc
	//2: Rotated Schofield
	int mu_code=0;

	double m0; double h0; double h3; double h5; double theta0;
	double g0; double g1; double g2; double g3;
	double c_star=1.0;

	std::unordered_map<double, double> mux;
	std::unordered_map<double, double> dmux;
	std::unordered_map<double, double> ddmux;

	py::array_t<double> thetas;
	py::array_t<double> Rs;

	//Window Function and derivatives (only used for mu_code=0)
	py::array_t<double> Wfs;
	py::array_t<double> Wf_Ts;
	py::array_t<double> Wf_ms;

	void set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) override;

	void clear_comp_flag() override;

	void set_hs(double h3, double h5);

	void mux_insert(double Temp,  double mu, double dmu, double ddmu);
	void get_mux(double Temp);
	void get_muxs();

	void get_thetas();


	void compute() override;

    //Used by pythons pickling and copying features
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void set_n0P0(double n_rat, double P_rat);
	double n_BG(double Temp, double mu_B);

	double h_func(double theta);
	double g_func(double theta);
	double dh_func(double theta);
	double dg_func(double theta);

	double obj_function(double theta);
	double get_theta(double x, double y);
};

void init_Scaling(py::module_ &m);

#endif