/*
	This file computes the EOS for for a gas of hadron resonances (HRG) with
	an excluded volume proportional to energy (EXI).
*/

#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>

#include "misc.hpp"
#include <stdexcept>
#include "PT.hpp"
#include "EXI.hpp"


using namespace std;

//-----------------------------------------------------------------------------
//Functions for computing Temp, mu_B in model exI given Temp* = T_star and
//mu_B* = mu_B_star.

double Temp_exI(double T_star, double mu_B_star, double epsilon0) {
	if(T_star <=0.0)
		throw domain_error("Temp_exI expects positive T_star.");

	if(epsilon0 <= 0.0)
		throw domain_error("Temp_exI expects positive epsilon0.");

	double Temp = T_star/(1.0 - 
	 Total_Pressure_pt(T_star, mu_B_star)/epsilon0);
	return Temp;
}

double mu_B_exI(double T_star, double mu_B_star, double epsilon0){

	if(T_star <=0.0)
		throw domain_error("mu_B_exI expects positive T_star.");

	if(epsilon0 <= 0.0)
		throw domain_error("mu_B_exI expects positive epsilon0.");

	double mu_B = mu_B_star/(1.0 - 
	 Total_Pressure_pt(T_star, mu_B_star)/epsilon0);
	return mu_B;
}


//-----------------------------------------------------------------------------

void Find_Temp_star_mu_B_star_exI(double Temp, double mu_B, double epsilon0, double& Tstar, double& mu_Bstar){
	//Compute Temp* (Tstar) and mu_B* (mu_Bstar) in model exI given Temp, mu_B.
	//returns values in reference variables T_star, mu_B_star.
	//Based on 3rd generation solver developed in python code.

	if(Temp <= 0.0) throw domain_error(
		"Find_Temp_star_mu_B_star_exI expects positive temperature.");

	if(epsilon0 <= 0.0)	throw domain_error(
		"Find_Temp_star_mu_B_star_exI expects positive epsilon0.");


	double Tol_Temp = 0.00000001;       //Tolerance to which we try to know Temp.
	//Bounds on Tstar search window. TstarMax found empirically, seems to work.
	double TstarMin = 0.0, TstarMax = 400.0;
	double TempLinearCutoff = 200.0; //Below here, Temp~Tstar. Found empirically.


  	double Tol_Tstar = 1e-10;   //If (TstarMax-TstarMin) < Tol_star, we say that
								//TstarMax == TstarMin == Tstar2, and we're done

	//-----------------------------------------------------------------
	//Simple initialization.  Could improve in future.

	double Tstar0, Tstar1;     //Temp* guesses
	double f0, f1;             //fcn values at Tstar0, Tstar1

	//generate 1st guess = Tstar0:

	if(Temp < TempLinearCutoff){  //for small Temps, Tstar ~ Temp, but Tstar <Temp
		Tstar0 = Temp - 1.0;   //so guess Tstar <~ Temp
		if(Tstar0 < TstarMin)   //revise guess to avoid values below TstarMin
			Tstar0 = (TstarMin + TstarMax)/2.0;
	}
	else    //Simply start in the middle of the Tstar search window
		Tstar0 = (TstarMin + TstarMax)/2.0;

	//evaluate fcn at Tstar0 guess.  fcn==0 when you find sln:
	f0 = Temp_exI(Tstar0, (mu_B/Temp)*Tstar0, epsilon0) - Temp;
	if(f0 > 0.0 || f0 < -Temp){ //Tstar0 is too big
		TstarMax = Tstar0;
		Tstar1 = Tstar0 - 1.0;   //get a second nearby guess
	}else{                          //Tstar0 is too small
		TstarMin = Tstar0;
		Tstar1 = Tstar0 + 1.0;
	}

	if(Tstar1 < TstarMin || Tstar1 > TstarMax)  //don't go out of bounds
		Tstar1 = (TstarMin + TstarMax)/2.0;

	//eval fcn at Tstar1.  fcn==0 when you find sln
	f1 = Temp_exI(Tstar1, (mu_B/Temp)*Tstar1, epsilon0) - Temp;

	//-----------------------------------------------------------------
	//Now, iterate, using previous 2 guesses in each step to pick next guess.
	//If a smart secant-method based guess fails, fall back to bisection method.

	bool UseSecantMethod = true; //use secant method to accelerate solving

	double Tstar2, divisor, f2, error;
	for(int iteration = 0; iteration < MAX_NUMBER_ITER; iteration++){

		//Secant method is not finding the solution quickly,
		// turn it off and use bisection method which is guaranteed to converge.
		if(iteration == MAX_NUMBER_ITER / 2) UseSecantMethod = false;

		//make new guess Tstar2
		if(UseSecantMethod == true){
			divisor = f1 - f0;
			Tstar2 = (Tstar0*f1 - Tstar1*f0)/divisor; //use secant method

			if(!(Tstar2 >= TstarMin && Tstar2 <= TstarMax)) //if new guess invalid,
				Tstar2 = (TstarMin + TstarMax)/2.0;   //fall back to bisection method
		}
		else Tstar2 = (TstarMin + TstarMax)/2.0; //use bisection method

		//find new f value at Tstar2:
		f2 = Temp_exI(Tstar2, (mu_B/Temp)*Tstar2, epsilon0) - Temp;
		error = fabs(f2);

		//Note: if fabs(TstarMax-TstarMin) < Tol_Tstar, you cannot
		//improve any further because to available accuracy
		//TstarMax == TstarMin, so we should quit.

		if( (error < Tol_Temp) || ( fabs(TstarMax-TstarMin) < Tol_Tstar)){
			Tstar = Tstar2;                  //found solution
			mu_Bstar = (mu_B/Temp)*Tstar;
			break;
		}

		//update variables for next iteration:
		if(f2 > 0.0 || f2 < -Temp) //Tstar2 is too big
			TstarMax = Tstar2;
		else  		//Tstar2 is too small
			TstarMin = Tstar2;
		Tstar0 = Tstar1;
		Tstar1 = Tstar2;
		f0 = f1;
		f1 = f2;
	}

	if( (error >= Tol_Temp) && (fabs(TstarMax-TstarMin) >= Tol_Tstar) )
		throw runtime_error("Function Find_Temp_star_mu_B_star_exI() failed to "
			" find a solution.");
}



//-----------------------------------------------------------------------------
//Code to compute thermal properties of a exI gas of many particle
//species


double Total_Pressure_exI(double Temp, double mu_B, double epsilon0){
	//Computes the total pressure in the EXI model.

	if(Temp <=0.0)
		throw domain_error("Total_Pressure_exI expects positive temperature.");

	if(epsilon0 <= 0.0)
		throw domain_error("Total_Pressure_exI expects positive epsilon0.");

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, T_star,
	 mu_B_star);

	double Press_pt = Total_Pressure_pt(T_star, mu_B_star);

	double Press_exI = Press_pt/(1.0 - Press_pt/epsilon0);

	return Press_exI;
}

double Total_Energy_Dens_exI(double Temp, double mu_B, double epsilon0){
	//Computes the total energy density in the EXI model.

	if(Temp <=0.0)
		throw domain_error("Total_Energy_Dens_exI expects positive"
		 " temperature.");

	if(epsilon0 <= 0.0)
		throw domain_error("Total_Energy_Dens_exI expects positive epsilon0.");

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, T_star,
	 mu_B_star);

	double EnergyDens_pt = Total_Energy_Dens_pt(T_star, mu_B_star);
	double EnergyDens_exI = EnergyDens_pt/(1.0 + EnergyDens_pt/epsilon0);
	return EnergyDens_exI;
}

double Total_Baryon_Dens_exI(double Temp, double mu_B, double epsilon0){
	//Computes the total baryon density in the EXI model.

	if(Temp <=0.0)
		throw domain_error("Total_Baryon_Dens_exI expects positive"
		 " temperature.");

	if(epsilon0 <= 0.0)
		throw domain_error("Total_Baryon_Dens_exI expects positive epsilon0.");


	if(mu_B==0) return 0.0;

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, T_star,
	 mu_B_star);

	double EnergyDens_pt = Total_Energy_Dens_pt(T_star, mu_B_star);
	double BaryonDens_pt = Total_Baryon_Dens_pt(T_star, mu_B_star);
	double BaryonDens_exI = BaryonDens_pt/(1.0 + EnergyDens_pt/epsilon0);
	return BaryonDens_exI;
}

double Total_Entropy_Dens_exI(double Temp, double mu_B, double epsilon0){
	//Computes the total entropy density in the EXI model.

	if(Temp <=0.0)
		throw domain_error("Total_Entropy_Dens_exI expects positive"
		 "temperature.");

	if(epsilon0 <= 0.0)
		throw domain_error("Total_Entropy_Dens_exI expects positive epsilon0.");

	double EnergyDens_exI = Total_Energy_Dens_exI(Temp, mu_B, epsilon0);
	double Press_exI = Total_Pressure_exI(Temp, mu_B, epsilon0);
	double BaryonDens_exI = Total_Baryon_Dens_exI(Temp, mu_B, epsilon0);
	//use standard thermodynamic relation
	return (EnergyDens_exI + Press_exI - mu_B*BaryonDens_exI)/Temp;
}


//----------------------------------------------------
//EXI model struct and pybind things

EXI::EXI(double eps0){
	this->eps0=eps0;
}

void EXI::compute_single(double Temp, double mu_B, double * y){

	if(isnan(eps0))
		throw invalid_argument("epsilon0 is nan in EXI.compute");
	y[0]=Total_Pressure_exI(Temp, mu_B,eps0);
	if(divs>0){
		y[2]=Total_Baryon_Dens_exI(Temp, mu_B,eps0);
		double E=Total_Energy_Dens_exI(Temp, mu_B,eps0);
		y[1]=(E+y[0]-mu_B*y[2])/Temp;
	} if(divs>1) {
		double dT=.1, dmu=.1;
		y[4]=.5*(Total_Baryon_Dens_exI(Temp+dT, mu_B,eps0)
			-Total_Baryon_Dens_exI(Temp-dT, mu_B,eps0))/dT;
		y[5]=.5*(Total_Baryon_Dens_exI(Temp, mu_B+dmu,eps0)
			-Total_Baryon_Dens_exI(Temp, mu_B-dmu,eps0))/dmu;
		double E_T=.5*(Total_Energy_Dens_exI(Temp+dT, mu_B,eps0)
			-Total_Energy_Dens_exI(Temp-dT, mu_B,eps0)) /dT;
		y[3]=(E_T -mu_B*y[4])/Temp;
	}
}

//Used by pythons pickling and copying features
py::tuple EXI::getstate(){
	return py::make_tuple(1,Therm::getstate(),eps0);
}

void EXI::setstate(py::tuple t){
	if (t.size() != 3)
		throw std::runtime_error("Invalid state!(EXI)");
	Therm::setstate(t[1]);
	eps0=t[2].cast<double>();
}

//used to include struct in pybind Module
void init_EXI(py::module_ &m){

	py::class_<EXI,Therm>(m, "EXI").def(py::init<>(),"EXI Object")
		.def(py::init<double>(),"EXI Object")
		.def("__repr__",
			[](const EXI &a) {
				return "<EXI: eps0 = " + to_string(a.eps0) + ">";
		}).def(py::pickle(
			[](EXI &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				EXI p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		).def_readwrite("eps0", &EXI::eps0);
}
