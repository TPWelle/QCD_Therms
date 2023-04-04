
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_multimin.h>  //for minimization
#include <cstdio>         //to write to file via fprintf
#include <omp.h>          //openmp for parallelization

#include "PT.hpp"
#include "EXI.hpp"
#include "EXII.hpp"
#include "pQCD.hpp"
#include "RMF.hpp"

#include "Standard.hpp"
#include "Crossover.hpp"

#include "optimize.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "latticedata.hpp"


using namespace std;
namespace py = pybind11;

struct chisquare_params{

	//pointer to the model to be optimized. Temps and mu_Bs
	//must already be set to appropriate lattice values.
	Therm *pth;

	//list of pointers to parameters in th
	//which are to be optimized
	vector<double *> opti;

	gsl_vector* pos_init;

	int N_Lats;
	vector<double> Press_Norm_Lat; 
	vector<double> dPress_Norm_Lat;
	vector<double> TraceAnomaly_Lat;
	vector<double> dTraceAnomaly_Lat;
};

/*
	This function is called by the gls minimizer. This function in turn calls
  my function FindChiSquare_TCrossover() to compute the chi-square.

  Also, this function enforces boundary conditions (variables non-negative).
  I put this feature in this caller function because some solvers allow you
  to specify constraints explicitly, so lack of constraints is a gsl issue,
  so keep added constraints specific to the gsl minimizer (ie this function).
 */

double Get_chisquare(const gsl_vector* v, void* p){
	
	chisquare_params par = *(chisquare_params *) p;

	int N=par.N_Lats;
	double chisquare = 0.0;

	Therm *pth = par.pth;

	gsl_vector* pos_init = par.pos_init;

	for(int i=0; i< par.opti.size();i++){
		if(gsl_vector_get(v,i)<0){
			chisquare+=1e8 + 100.0*fabs(gsl_vector_get(v,i));
			cout<<i<<"  "<<gsl_vector_get(v,i)<<" "<<chisquare<<endl;
		} else {
			*(par.opti[i]) = gsl_vector_get(v,i)*
				gsl_vector_get(pos_init,i);
		}
	}

	pth->clear_comp_flag();
	pth->compute();

	double *pTemp= get_arr_pointer(pth->Temps), *pmu= get_arr_pointer(pth->mu_Bs);

	double *pY= get_arr_pointer( pth->Y), *pY_T= get_arr_pointer( pth->Y_T ),
	 *pY_m= get_arr_pointer( pth->Y_m );
	double *pY_TT= get_arr_pointer( pth->Y_TT ), *pY_Tm= get_arr_pointer( pth->Y_Tm ),
	 *pY_mm= get_arr_pointer( pth->Y_mm );

	double Press,TA,Temp;

	for(int i=0;i<N;i++){
		Temp=pTemp[i];
		Press=pY[i];
		TA=Temp*(pY_T[i]) +(pmu[i])*(pY_m[i])-4*Press;
		//Add terms to chi-square from pressure/(Temp^4) and TraceAnomaly
		chisquare += pow( (Press/(pow(Temp, 4.0)) - par.Press_Norm_Lat[i])
		 	/ par.dPress_Norm_Lat[i] , 2.0) +
		  pow( (TA/( pow(Temp, 4.0 ) )- par.TraceAnomaly_Lat[i])
			/ par.dTraceAnomaly_Lat[i] , 2.0);

		if(!gsl_finite(chisquare)) {
			printf("-------------------\n%f %f %f %f\n"
				,chisquare,Temp,Press,TA+3.0*Press);
		}
	}
	return chisquare;
}


//modifies th to be the optimal parameters, returns chi^2
double MinimizeChiSquare(Therm *pth, const vector<double *> opti){

	//gsl parameters
	const double epsabs          = 1.0e-6;//absolute error in minimum position
	const int MAX_NUM_ITERATIONS = 300;


	int N_Lats=Lattice_Temps.size(); //number of lattice points

	pth->Temps = py::array_t<double>(N_Lats);
	pth->mu_Bs = py::array_t<double>(N_Lats);

	double *pTemps= get_arr_pointer(pth->Temps), *pmu= get_arr_pointer(pth->mu_Bs);


	for(int i=0;i<N_Lats;i++){
		pTemps[i]=Lattice_Temps[i];
		pmu[i]=0.0;
	}

	double chisquare;  //chisquares from calculations
	int NumIterations; //store how many iterations used

	int status;

	pth->set_Tmu(pth->Temps, pth->mu_Bs);

	//creating parameter structure for minimization
	chisquare_params par;

	par.pth=pth;
	par.opti=opti;
	par.N_Lats=N_Lats;
	par.Press_Norm_Lat=PressPerT4;
	par.dPress_Norm_Lat=dPressPerT4;
	par.TraceAnomaly_Lat=TraceAnomaly;
	par.dTraceAnomaly_Lat=dTraceAnomaly;

	//Set type of minimizer to Nelder-Mead Simplex, version 2
	const gsl_multimin_fminimizer_type* min_type =
	  gsl_multimin_fminimizer_nmsimplex2;

	//create variable to track the state of the minimizer
	gsl_multimin_fminimizer* state = NULL;
	state = gsl_multimin_fminimizer_alloc(min_type, opti.size()); //initialize

	// position tracks values of fcn arguments as optimizer minimizes fcn
	gsl_vector* position  = gsl_vector_alloc(opti.size());

	//The initial values of the parameters being optimized
	//Used to scale variables
	gsl_vector* pos_init = gsl_vector_alloc(opti.size());

	// Holds gsl step sizes
	gsl_vector* stepsize = gsl_vector_alloc(opti.size());


	for(int k=0; k<opti.size(); k++){
		if( *(opti[k]) != 0.0) gsl_vector_set(pos_init, k, *(opti[k]));
		else gsl_vector_set(pos_init, k, 1.0);
	}

	//set position to initial values. Set each stepsize to .2
	gsl_vector_set_all(position, 1.0);
	gsl_vector_set_all(stepsize, .2);

	par.pos_init=pos_init;

	//gsl structure to store info about the function that gsl will minimize
	gsl_multimin_function gslmin_func;
	gslmin_func.n = opti.size();
	gslmin_func.f = Get_chisquare;
	gslmin_func.params = &par;

	int iter = 0;  //current iteration step
	double size = 0.0;
	gsl_set_error_handler_off();

	//setup problem given minimizer, state, function, position, and stepsize
	gsl_multimin_fminimizer_set(state, &gslmin_func, position, stepsize);
	do {
		iter++;
		//Iterate the solver and return the status.
		status = gsl_multimin_fminimizer_iterate(state);
		//get simplex size given the current state size and
		//check if the size is smaller than the tolerance.
		//if size < epsabs, status == GSL_SUCCESS == 0
		if (status != GSL_SUCCESS){
			printf("Error: gsl_multimin_fminimizer_iterate returned error code"
				" %d in iteration %d.\n", (int)status,iter);
			break;
  		}

		size = gsl_multimin_fminimizer_size(state);
		status = gsl_multimin_test_size(size, epsabs);

		//Print position and status
		cout << iter <<" | ";
		for(int k=0; k<opti.size(); k++){
			cout <<  gsl_vector_get(state->x,k)*gsl_vector_get(pos_init,k)<<" ";
		}
		cout << " | " << (state->fval)/N_Lats<<" "<< size/epsabs <<endl;

	} while (status == GSL_CONTINUE && iter < MAX_NUM_ITERATIONS);
	//get values out of v and put in parameters_new
	chisquare = state->fval;
	NumIterations = iter;

	cout << "=================================="<<endl;
	//Print position and status, set the values in the model
	cout  << NumIterations <<" | ";
	for(int k=0; k<opti.size(); k++){
		double X = gsl_vector_get(state->x, k)*gsl_vector_get(pos_init, k);
		*(opti[k]) = X;
		cout << X <<" ";
	}
	cout << " | " << chisquare/N_Lats<<" "<< size/epsabs <<endl;
	cout << "=================================="<<endl;

	//cleanup gsl
	gsl_vector_free(position);
	gsl_vector_free(stepsize);
	gsl_vector_free(pos_init);
	gsl_multimin_fminimizer_free(state);
	
	return chisquare/N_Lats;
}


vector<double> Optimize(Crossover* pth){

	//Vector to store pointers to the parameters being optimized.
	vector<double *> opti;

	Standard * pSF = dynamic_cast<Standard*>( pth->pSF );
	opti.push_back( &(pSF->T0) );

	if (PT *pP1 = dynamic_cast<PT*> (pth->pP1) ){
		//Do Nothing
	}else if (EXI *pP1 = dynamic_cast<EXI*> (pth->pP1) ){
		opti.push_back( &(pP1->eps0) );
	}else if (EXII *pP1 = dynamic_cast<EXII*> (pth->pP1) ){
		opti.push_back( &(pP1->eps0) );
	}else if ( RMF *pP1 = dynamic_cast<RMF*> (pth->pP1) ){
		// opti.push_back( &(pP1->g_sig) );
		// opti.push_back( &(pP1->g_om) );
		// opti.push_back( &(pP1->g_rho) );
		// opti.push_back( &(pP1->alpha_V) );
		// opti.push_back( &(pP1->z) );
		// opti.push_back( &(pP1->b) );
		// opti.push_back( &(pP1->c) );

		//Do Nothing
	}else throw invalid_argument("Unrecognized Model in Optimize()");


	pQCD* pP2 = dynamic_cast<pQCD*>( pth->pP2 );
	opti.push_back( &(pP2->EScl) );
	opti.push_back( &(pP2->Soft) );

	double chisquare = MinimizeChiSquare(pth, opti);
	vector<double> out(opti.size()+1,0.0);

	out[0] = chisquare;

	for(int k=0; k<opti.size(); k++){
		out[k+1] = *(opti[k]);
	}

	return out;
}





void init_optimize(py::module_ &m){
	m.def("MinimizeChiSquare", &MinimizeChiSquare, "A function for chisquare");
	// m.def("Optimize_PT", &Optimize_PT, "A funcion to optimize PT");
	// m.def("Optimize_EXI", &Optimize_EXI, "A funcion to optimize EXI");
	// m.def("Optimize_EXII", &Optimize_EXII, "A funcion to optimize EXII");
	m.def("Optimize", &Optimize, "A funcion to optimize Crossover");
}



//Here are list of gsl codes to help diagnose gsl errors:
/*
enum  	{
  GSL_SUCCESS = 0, GSL_FAILURE = -1, GSL_CONTINUE = -2, GSL_EDOM = 1,
  GSL_ERANGE = 2, GSL_EFAULT = 3, GSL_EINVAL = 4, GSL_EFAILED = 5,
  GSL_EFACTOR = 6, GSL_ESANITY = 7, GSL_ENOMEM = 8, GSL_EBADFUNC = 9,
  GSL_ERUNAWAY = 10, GSL_EMAXITER = 11, GSL_EZERODIV = 12, GSL_EBADTOL = 13,
  GSL_ETOL = 14, GSL_EUNDRFLW = 15, GSL_EOVRFLW = 16, GSL_ELOSS = 17,
  GSL_EROUND = 18, GSL_EBADLEN = 19, GSL_ENOTSQR = 20, GSL_ESING = 21,
  GSL_EDIVERGE = 22, GSL_EUNSUP = 23, GSL_EUNIMPL = 24, GSL_ECACHE = 25,
  GSL_ETABLE = 26, GSL_ENOPROG = 27, GSL_ENOPROGJ = 28, GSL_ETOLF = 29,
  GSL_ETOLX = 30, GSL_ETOLG = 31, GSL_EOF = 32
}
*/
