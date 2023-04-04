/*
This file computes the EOS for for a gas of point-like (PT) hadron resonances (HRG).
*/
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>

#include <gsl/gsl_errno.h>          //for modifying gsl's error handling
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_result.h>      //contains a struct for returning results 
                                    //of gsl functions
#include <stdio.h>
#include "misc.hpp"
#include "particle.hpp"
#include "particlelist.hpp"
#include "PT.hpp"

using namespace std;


//ignore errors if numerical integration has problems
#define SUPPRESS_INTEGRATION_ERRORS   

//-----------------------------------------------------------------------------
//Code to compute (partial) thermal properties of one species of gas:


double equilib_distfcn(double x, double y, double z, STATISTICS stats){
	//Returns the value of Bose-Einstein, Fermi, or Boltzmann distribution
	//function. Depends on the type of particle, the energy =mass*x,
	//Temp=mass*y, and particle chemical potential mu=mass*z.

	if(y <= 0.0)
		throw domain_error("equilib_distfcn expects positive temperature.");

	const double EXP_ARG_MAX = 700.0, EXP_ARG_MIN = -700.0;
	double argument =(x-z)/y;

	//avoid overflow and underflow
	if(argument > EXP_ARG_MAX) return 0.0;
	else if(argument < EXP_ARG_MIN) return (1.0/stats);   
	//stats==0 can cause problems, avoid by not using classical stats at low Temp
	else return 1.0/(gsl_sf_exp(argument) + stats);
}

//store parameters passed to Thermal_Integrand() via void *p
struct Integration_Params {
	double y; // Temp/mass
	double z; // mu/mass
	STATISTICS stats;
	int a; //power of x in integrand
	int b; //power of sqrt(x^2 -1) in integrand
	// integers a and b specify what is being computed
	//   P: a=0, b=3
	//   e: a=2, b=1
	//   n: a=1, b=1
	// n_s: a=0, b=1
};

double Thermal_Integrand(double x, void *p){
	//Returns dimensionless integrand as a fcn of x = E/mass
	//and other parameters (like Temp, mass, particle chemical potential,...)
	//When integrated over x, this integral will compute the
	//dimensionless component of various quantities.

	Integration_Params * params = (Integration_Params *) p;

	if(x<1.0) throw domain_error("Thermal_Integrand() requires"
		 "x = E/mass >= 1.");

	return pow(x*x - 1.0, .5*params->b)*pow(x,params->a)*
	 equilib_distfcn(x, params->y, params->z, params->stats);
}


double Thermal_Integral(double y, double z, int a, int b, STATISTICS stats){
	//Compute a thermodynamic integral for a particle species in PT model,
	//identified by particle object, at temperature Temp and
	//particle chemical potential mu_particle.

	if(y <=0.0) throw domain_error("Theramal_Integral_pt expects"
		" positive temperature.");

	double lower_limit = 1.0;  //integrate x=E/mass from 1 to infinity

	//do some setup for gsl
	const int MaxNumIntervals = 1000;
	gsl_integration_workspace *work_ptr =
	 gsl_integration_workspace_alloc(MaxNumIntervals);
	double abs_error = 1.0e-8; double rel_error = 1.0e-8;
	double result, error;		//where we will store results

	//create structure holding parameters to pass to integrand function:
	Integration_Params params;
	params.y = y;  params.z = z;
	params.stats = stats;
	params.a=a;  params.b=b;
	void *params_ptr = &params;

	//set up reference to the integrand function to be integrated over:
	gsl_function My_function;
	My_function.function = &Thermal_Integrand;
	My_function.params = params_ptr;

	//Disable gsl's error handling so it won't crash the code if
	//(for instance) an integral converges slowly. We'll check status after.
	gsl_set_error_handler_off();

	//invoke gls's quadrature-adaptive routine, integrate from lower limit
	//to infinity
	int return_code=gsl_integration_qagiu(&My_function, lower_limit,
	 abs_error, rel_error, MaxNumIntervals, work_ptr, &result, &error);
	//Turn gsl's default error handler back on
	//gsl_set_error_handler(NULL);

	gsl_integration_workspace_free(work_ptr);

	// If SUPPRESS_INTEGRATION_ERRORS is defined we
	// allow some non-critical errors:

	// Otherwise, throw an exception to let user know an unacceptable
	// gsl error occured.

	#ifndef SUPPRESS_INTEGRATION_ERRORS
	if(return_code != 0) throw runtime_error("gsl fcn gsl_integration_qagiu()"
		 " experienced unexpected error in fcn Thermal_Integral()."
		 " gsl error code was "+ to_string(return_code)+ ". y = " + to_string(y)+
		 ". z = " + to_string(z) + ". a,b = "+to_string(a)+to_string(b));

	#endif
	//don't allow negative pressure or energy density.
	if(result < 0.0 and not ((a==0 or a==1) and b==1)) {
		result=0;
		// cout<<"negative value warning: result="<< result
		//  <<" a="<<a<<" b="<<b<<endl;
		if(result<-100.0) throw runtime_error("Error: Thermal_Integral() computed"
		 " a negative result. res=" +to_string(result) + " y = " + to_string(y) + " z = "
		 + to_string(z)+ ". a,b = "+to_string(a)+to_string(b)
		 + ". return_code = "+to_string(return_code));
	}
	return result;
}


double Partial_Pressure_pt(double Temp, double mu_particle,
 const Particle& particle){
	//Compute partial pressure of a particle species in PT model,
	//identified by particle object, at temperature Temp and
	//particle chemical potential mu_particle.

	if(Temp <=0.0) throw domain_error("Partial_Pressure_pt expects"
		" positive temperature.");

	double mass = particle.mass;

	STATISTICS stats = particle.stats;

	if(stats == CLASSICAL){
		return 1.0/(2.0*pi*pi)*(mass*mass)*(Temp*Temp)*
		  gsl_sf_bessel_Kn(2, mass/Temp)*gsl_sf_exp(mu_particle/Temp);
	} else {

		double result; //where we will store results
		result = Thermal_Integral( Temp/mass, mu_particle/mass,
			0, 3, stats);

		//multiply integral by phase space factors, and units
		//which convert integral into pressure
		return result * pow(mass, 4.0)/(6.0*pi*pi);
	}
}

//------------------------------------------------------

double Partial_Number_Dens_pt(double Temp, double mu_particle,
 const Particle& particle){
	//Compute partial number density of a particle species per DOF,
	//identified by particle object at temperature Temp and particle
	//chemical potential mu_particle.

	if(Temp <=0.0) throw domain_error("Partial_Number_Dens_pt "
		"expects positive temperature.");

	double mass = particle.mass;
	STATISTICS stats = particle.stats;

	if(stats == CLASSICAL)	{
		//valid for classical stats only:
		return Partial_Pressure_pt(Temp, mu_particle, particle)/Temp;
	} else {
		double result; //where we will store results
		result = Thermal_Integral( Temp/mass, mu_particle/mass,
			1, 1, stats);

		//multiply integral by degeneracy, phase space factors, and units 
		//which convert integral into number density 
		return result * pow(mass, 3.0)/(2.0*pi*pi);
	}
}

//--------------------------------------------------

double Partial_Energy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle){
	//Compute partial energy density of a particle species, identified by 
	//particle object at temperature Temp and particle chemical potential 
	//mu_particle.

	if(Temp <=0.0) throw domain_error("Partial_Energy_Dens_pt"
		 " expects positive temperature.");
  
	double mass = particle.mass;
	STATISTICS stats = particle.stats;

	if(stats == CLASSICAL) {
		return 1.0/(2.0*pi*pi)*pow(mass,4.0)*gsl_sf_exp(mu_particle/Temp)*
			(3.0*gsl_sf_bessel_Kn(2,mass/Temp)/(pow(mass/Temp,2)) +
			gsl_sf_bessel_Kn(1,mass/Temp)/(mass/Temp));
	} else {
		double result; //where we will store results
		result = Thermal_Integral( Temp/mass, mu_particle/mass,
			2, 1, stats);

		//multiply integral by phase space factors, and units
		//which convert integral into energy density
		return result * pow(mass, 4.0)/(2.0*pi*pi) ;
	}
}

double Partial_Entropy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle){
	if(Temp <=0.0) throw domain_error("Partial_Entropy_Dens_pt"
		" expects positive temperature.");

	double energydens = Partial_Energy_Dens_pt(Temp, mu_particle, particle);
	double pressure = Partial_Pressure_pt(Temp, mu_particle, particle);
	double numberdens = Partial_Number_Dens_pt(Temp, mu_particle, particle);

	return (energydens + pressure - mu_particle*numberdens)/Temp;
}


//=============================================================================
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_pt(double Temp, double mu_B) {

	if(Temp <=0.0) throw domain_error("Total_Pressure_pt"
		" expects positive temperature.");

	double mu_particle;
	double Pressure=0.0;
	int degen;

	for(auto const &particle: ParticleList) {

		if(mu_B==0.0){
			Pressure += particle.degen() *
				Partial_Pressure_pt(Temp, 0.0, particle);
		} else {
			// Sum over possible values for baryon number
			for(int B : {-1,0,1}){ //for B=-1,0,1
				degen=particle.degen(B);
				if(degen !=0) Pressure += degen *
					Partial_Pressure_pt(Temp, B*mu_B, particle);
			}
		}

	}
	return Pressure;
}

double Total_Energy_Dens_pt(double Temp, double mu_B){

	if(Temp <=0.0) throw domain_error("Total_Energy_Dens_pt"
		" expects positive temperature.");

	double mu_particle;
	double EnergyDens=0.0;
	int degen;

	for(auto const &particle: ParticleList) {

		if(mu_B==0.0){
			EnergyDens += particle.degen() *
				Partial_Energy_Dens_pt(Temp, 0.0, particle);
		} else {
			// Sum over possible values for baryon number
			for(int B : {-1,0,1}){ //for B=-1,0,1
				degen=particle.degen(B);
				if(degen !=0) EnergyDens += degen *
					Partial_Energy_Dens_pt(Temp, B*mu_B, particle);
			}
		}
	}
	return EnergyDens;
}

double Total_Baryon_Dens_pt(double Temp, double mu_B){

	if(Temp <=0.0) throw domain_error("Total_Baryon_Dens_pt"
		" expects positive temperature.");

	if(mu_B==0.0) return 0.0;

	double mu_particle;
	double baryonDens=0.0;
	int degen;

	for(auto const &particle: ParticleList) {
		// Sum over possible values for baryon number
		for(int B : {-1,1}){
			degen=particle.degen(B);
			if(degen != 0) baryonDens += degen*B*
				Partial_Number_Dens_pt(Temp, B*mu_B, particle);
		}
	}
	return baryonDens;
}

double Total_Entropy_Dens_pt(double Temp, double mu_B){

	if(Temp <=0.0) throw domain_error("Total_Entropy_Dens_pt"
		" expects positive temperature.");

	double energyDens = Total_Energy_Dens_pt(Temp, mu_B);
 	double pressure = Total_Pressure_pt(Temp, mu_B);
	double baryonDens = Total_Baryon_Dens_pt(Temp, mu_B);
	//use standard thermodynamic relation 
	return (energyDens + pressure - mu_B*baryonDens)/Temp;
}

//----------------------------------------------------
//PT model struct and pybind things

void PT::compute_single(double Temp, double mu_B, double * y){
	y[0]=Total_Pressure_pt(Temp, mu_B);
	if(divs>0){
		double E = Total_Energy_Dens_pt(Temp, mu_B);
		if (mu_B==0) y[2]= 0.0;
		else y[2]=Total_Baryon_Dens_pt(Temp, mu_B);
		y[1]=(E+y[0] -mu_B*y[2])/Temp;
	} if(divs>1) {
		double dT=.1, dmu=.1;
		if (mu_B==0) y[4]= 0.0;
		else y[4]=.5*(Total_Baryon_Dens_pt(Temp+dT, mu_B) -Total_Baryon_Dens_pt(Temp-dT, mu_B))/dT;
		y[5]=.5*(Total_Baryon_Dens_pt(Temp, mu_B+dmu) -Total_Baryon_Dens_pt(Temp, mu_B-dmu))/dmu;
		double E_T=.5*(Total_Energy_Dens_pt(Temp+dT, mu_B)-Total_Energy_Dens_pt(Temp-dT, mu_B))/dT;
		y[3]=(E_T-mu_B*y[4])/Temp;
	}
}


namespace py = pybind11;

//Used by pythons pickling and copying features
py::tuple PT::getstate(){
	return py::make_tuple(0,Therm::getstate());
}
void PT::setstate(py::tuple t){
	if (t.size() != 2)
		throw runtime_error("Invalid state!(PT)");
	Therm::setstate(t[1]);
}

//used to include struct in pybind Module
void init_PT(py::module_ &m){
	py::class_<PT,Therm>(m, "PT").def(py::init<>(),"PT model")
		.def("__repr__", [](const PT &a) {
				return "<PT>";
		}).def(py::pickle(
		[](PT &p){
			return p.getstate();
		}, [](py::tuple t){ // __setstate__
			PT p; /* Create a new C++ instance */
			p.setstate(t); /* Assign additional state */
			return p;
		})
	);
}

