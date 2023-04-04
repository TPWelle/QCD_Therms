/*
	This file computes the EOS for for a gas of hadron resonances (HRG) with
	an excluded volume proportional to mass (EXII).
*/

#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdexcept>
#include "misc.hpp"
#include "particle.hpp"
#include "PT.hpp"
#include "EXI.hpp"
#include "EXII.hpp"

using namespace std;
//-----------------------------------------------------------------------------
//Code to find exII pressure by solving a non-linear equation.


double PexII_fcn_to_zero(double Press_exII, double Temp, double mu_B,
 double epsilon0 ){
 	//When we find a correct value Press_exII, this function will equal zero.

	double mu_effective;
	double Sum = 0;
	int degen;
	for(auto const &particle: ParticleList)	{
		//particle's excluded volume:
		double volume = particle.mass / epsilon0;

		if(mu_B==0.0){
			Sum += particle.degen() *Partial_Pressure_pt(Temp,
				-Press_exII*volume, particle);
		} else {
		// Sum over possible values for baryon number
			for(int B : {-1,0,1}){
				degen=particle.degen(B);
				if(degen !=0) {
					mu_effective =B*mu_B - Press_exII*volume;
					Sum += degen *Partial_Pressure_pt(Temp,
						mu_effective, particle);
				}
			}
		}
	}
	//this = 0 when solution is found
	return (Press_exII - Sum);
}


double Total_Pressure_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess){
	//Compute the total pressure in the EXII model.

	const double Tol = 0.1;  //solver tolerance
	
	if(Temp <=0.0)	throw domain_error("Total_Pressure_exII "
		"expects positive temperature.");

	//default value for Press_exII_guess is -1.0 (means ignore)
	if(Press_exII_guess < 0.0) //no exII pressure guess provided, so generate one
		Press_exII_guess = Total_Pressure_exI(Temp, mu_B, epsilon0);

	//intialize guesses for rootfinding:
	//P is for pressure, f is for fcn to zero:  f(P) == 0 when find root P
	double P0 = Press_exII_guess; //1st guess
	double dP = 100.0;
	double P1 = P0 + dP;          //2nd guess: pick a nearby starting point

	double f0 = PexII_fcn_to_zero(P0, Temp, mu_B, epsilon0);
	double f1 = PexII_fcn_to_zero(P1, Temp, mu_B, epsilon0);

	double P2, f2, divisor, error;

	//loop and perform search, updating guesses, until find solution.
	//Uses secant method.
	for(int i = 0; i < 1*MAX_NUMBER_ITER; i++){
		divisor = (f1 - f0);
		if(divisor != 0.0)
			P2 = P1 - f1*(P1 - P0)/divisor;  //secant method
		else
			P2 = Press_exII_guess*(1.0 + 0.2*i);  //fail, so shift starting point

		f2 = PexII_fcn_to_zero(P2, Temp, mu_B, epsilon0);
		error = std::abs(f2);
		if(error <= Tol) break;   //found solution within tolerance
		//update variables for next iteration
		P0 = P1; f0 = f1;
		P1 = P2; f1 = f2;
		//printf("%f %f %f \n",P0/pow(Temp,4),Press_exII_guess/pow(Temp,4),error);
	}if(error > Tol){
		printf("\n %8.8e \n",f2);
		printf("error: %f %f %f %f %e\n",P0,error,Temp,mu_B,epsilon0);
		throw runtime_error("Function Total_Pressure_exII failed to find"
			" a solution.");
    }
	return P2;  //this is the exII pressure which is a solution
}

//-----------------------------------------------------------------------------

double Total_Energy_Dens_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess){
	//Compute the total energy density in the EXII model.

	if(Temp <=0.0) throw domain_error("Total_Energy_Dens_exII"
		" expects positive temperature.");

	double mu_effective, volume;
	double numerator = 0.0, denominator = 1.0;
	int degen;
	//need exII pressure to get effective chem potential.Use guess to speed up.
	double Press_exII = Total_Pressure_exII(Temp, mu_B, epsilon0,
	 Press_exII_guess);

	for(auto const &particle: ParticleList)	{
		//particle's excluded volume:
		volume = particle.mass / epsilon0;

		if(mu_B==0.0){
			degen=particle.degen();

			numerator += degen*Partial_Energy_Dens_pt(Temp,
				- Press_exII*volume, particle);

			denominator += volume * degen * Partial_Number_Dens_pt(
				Temp, - Press_exII*volume, particle);

		} else {
			// Sum over possible values for baryon number
			// for(int B=-1; B<2; B++){ //for B=-1,0,1
			for(int B : {-1,0,1}){
				degen=particle.degen(B);
				if(degen !=0) {
					mu_effective =B*mu_B - Press_exII*volume;

					numerator += degen*Partial_Energy_Dens_pt(Temp,
						mu_effective, particle);

					denominator += volume * degen * Partial_Number_Dens_pt(
						Temp, mu_effective, particle);
				}
			}
		}
	}

	return numerator/denominator;
}

//----------------------------------------------------------------------------

double Total_Baryon_Dens_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess){
	//Compute the total baryon density in the EXII model.

	if(Temp <=0.0) throw domain_error("Total_Baryon_Dens_exII "
		"expects positive temperature.");

	if(mu_B==0) return 0.0;

	double mu_effective, volume;
	double numerator = 0.0, denominator = 1.0;
	int degen;
	//need exII pressure to get effective chem potential.Use guess to speed up.
	double Press_exII = Total_Pressure_exII(Temp, mu_B, epsilon0,
	 Press_exII_guess);

	for(auto const &particle: ParticleList){
		//particle's excluded volume:
		volume = particle.mass / epsilon0;

		for(int B : {-1,0,1}){
			degen=particle.degen(B);
			if(degen !=0) {
				mu_effective =B*mu_B - Press_exII*volume;

				if(B !=0 ) numerator += degen*B*Partial_Number_Dens_pt(
					Temp, mu_effective, particle);

				denominator += volume * degen * Partial_Number_Dens_pt(
					Temp, mu_effective, particle);
			}
		}
	}
	return numerator/denominator;
}


//----------------------------------------------------------------------------

double Total_Entropy_Dens_exII(double Temp, double mu_B, double epsilon0,
 double Press_exII_guess){
	//Compute the total entropy density in the EXII model.

	if(Temp <=0.0)
		throw domain_error("Total_Entropy_Dens_exII expects positive"
		 " temperature.");

	//find pressure, using guess to speed it up
	double Pressure = Total_Pressure_exII(Temp, mu_B, epsilon0,
 	 Press_exII_guess);

	//find other quantities, using pressure as a pressure guess to speed
	//the solution
	double EnergyDens = Total_Energy_Dens_exII(Temp, mu_B, epsilon0, Pressure);
	double BaryonDens = Total_Baryon_Dens_exII(Temp, mu_B, epsilon0, Pressure);
	double EntropyDens = (Pressure + EnergyDens - mu_B*BaryonDens)/Temp;
	return EntropyDens;
}


//----------------------------------------------------
//EXII model struct and pybind things

EXII::EXII(double eps0){
	this->eps0=eps0;
}

void EXII::compute_single(double Temp, double mu_B, double * y){

	if(isnan(eps0))
		throw invalid_argument("epsilon0 is nan in EXII.compute");
	y[0]=Total_Pressure_exII(Temp, mu_B,eps0);
	if(divs>0){
		y[2]=Total_Baryon_Dens_exII(Temp, mu_B,eps0,y[0]);

		double E=Total_Energy_Dens_exII(Temp, mu_B,eps0,y[0]);
		y[1]=(y[0] + E - mu_B*y[2])/Temp;
	} if(divs>1) {
		double dT=.1, dmu=.1;
		y[4]=.5*(Total_Baryon_Dens_exII(Temp+dT, mu_B,eps0) -Total_Baryon_Dens_exII(Temp-dT, mu_B,eps0))/dT;
		y[5]=.5*(Total_Baryon_Dens_exII(Temp, mu_B+dmu,eps0) -Total_Baryon_Dens_exII(Temp, mu_B-dmu,eps0))/dmu;
		double E_T=.5*(Total_Energy_Dens_exII(Temp+dT, mu_B,eps0)
			-Total_Energy_Dens_exII(Temp-dT, mu_B,eps0)) /dT;

		y[3]=(E_T-y[4]*mu_B)/Temp;
	}
}

//Used by pythons pickling and copying features
py::tuple EXII::getstate(){
	return py::make_tuple(2,Therm::getstate(),eps0);
}

void EXII::setstate(py::tuple t){
	if (t.size() != 3)
		throw std::runtime_error("Invalid state!(EXII)");

	Therm::setstate(t[1]);
	eps0=t[2].cast<double>();
}

//used to include struct in pybind Module
void init_EXII(py::module_ &m){
	py::class_<EXII,Therm>(m, "EXII").def(py::init<>(),"EXII Object")
		.def(py::init<double>(),"EXII Object")
		.def("__repr__",
	        [](const EXII &a) {
	            return "<EXII: eps0 = " + std::to_string(a.eps0) + ">";
	        }
        ).def(py::pickle(
			[](EXII &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				EXII p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		).def_readwrite("eps0", &EXII::eps0);
}
