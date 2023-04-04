#ifndef THERM_ARRAY_H
#define THERM_ARRAY_H

#include <cmath>
#include "misc.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//Gets pointer to access elements of PyBind arrays
double * get_arr_pointer(py::array_t<double> Y);

//Stores an array of values running over Temperature
//and chemical potential and optionally its first
//2 orders of derivatives.
struct Therm{

	int divs=2; //order of derivatives to be computed
	int N_threads=7; //Number of threads for parallelization
	bool comp_flag=false; //used to avoid computing values repeatedly
	bool display_progress=false;

	//Arrays of temperatures and chemical potentials
	py::array_t<double> Temps;
	py::array_t<double> mu_Bs;

	//Y stores an array of values (typically Pressures) at the
	//corresponding Temperature and chemical potential. Y_T, Y_m are
	//that values Temp and mu derivatives. Likewise Y_TT, Y_Tm, Y_mm.
	py::array_t<double> Y;
	py::array_t<double> Y_T;
	py::array_t<double> Y_m;

	py::array_t<double> Y_TT;
	py::array_t<double> Y_Tm;
	py::array_t<double> Y_mm;

	//Sets the arrays Temps and mu_Bs and allocates Y, Y_T, etc. 
	//to arrays of the same size.
	virtual void set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs);

	virtual void clear_comp_flag();

	//compute the values of Y at the the specified values of Temp and mu_B
	//also compute appropriate derivatives if divs>0.
	virtual void compute();
	virtual ~Therm()=default;

	//Used by pythons pickling and copying features
	virtual py::tuple getstate();
	virtual void setstate(py::tuple t);

	//computes a value and derivatives at a single point Temp,mu_B
	virtual void compute_single(double Temp, double mu_B,double * y);
};

//used to include struct in pybind Module
void init_Therm(py::module_ &m);

#endif