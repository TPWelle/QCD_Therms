/*
	This file computes the EOS for QCD with mass-less quarks
	and gluons including interactions.
*/


#include <cmath>
#include <stdexcept>

#include "pQCD.hpp"
#include "misc.hpp"

using namespace std;

/*Lambda is the energy scale.  This is formula 9.5 from the pdg review
	of QCD.  (My Lambda is pdg's mu_R.)  3-loop, so set b3=0.
	SoftenConstant == 0 is standard QCD--non-zero values soften divergences
	at small temperature.  */
double alphas_3loop(double Lambda, double SoftenConstant)
{
	double LambdaMS = 290.0;   //MeV
	double Nf = 3.0;

	double t = log(SoftenConstant + pow(Lambda/LambdaMS, 2.0) );

	double b0 = (33.0 - 2.0*Nf)/(12.0*pi);
	double b1 = (153.0 - 19.0*Nf)/(24.0*pi*pi);
	double b2 = (2857.0 - 5033.0/9.0*Nf + 325.0/27.0*Nf*Nf)/(128.0*pi*pi*pi);
	//double b2 = 0.0;

//To improve:  check for t == 0.0, because 1/t == infinity == trouble
//I could throw an exception (which would break execution), or return 0,
//or return a large fixed value.  Think further
	double alphas = 1.0/(b0*t)*(1.0 - b1*log(t)/(b0*b0*t)
	 + (b1*b1*(log(t)*log(t) - log(t) - 1.0 ) + b0*b2)/( pow(b0, 4.0)*t*t )
	 - ( b1*b1*b1*( pow(log(t),3.0) - 5.0/2.0*log(t)*log(t) - 2.0*log(t)
	 + 1.0/2.0  ) + 3.0*b0*b1*b2*log(t) )/( pow(b0*b0*t, 3.0) )    );

	return alphas;
}


//-----------------------------------------------------------------------------

//QCD pressure for free quarks and gluons
double Total_Pressure_Free_QCD(double Temp, double mu_B){
	int Nf=3;
	return pi*pi*(8.0/45.0 + 7.0*Nf/60.0)*pow(Temp,4)
		+.5*Nf*pow(Temp*mu_B,2)+.25*pow(mu_B,4)*Nf/(pi*pi);
}


/*Computes the pressure of QCD at finite Temp, mu_B using order
	(alpha_s)^3*log(alpha_s) perturbative results.  Formula is from eqs 1-7
	from arXiv:1212.1797v2 by Strickland which is simply quoting the result
	calulated by Vuorinen.  For Vuorinen's derivation, see ref 17 of
	arXiv:1212.1797v2.  (That's arXiv:hep-ph/0212283)
	I used the result quoted by Strickland, instead of the
	original work by Vuorinen, because Vuorinen did not finish dimensional
	regularization in his expression.*/
double Total_Pressure_PQCD(double Temp, double mu_B,
 double EnergyScaleConstant, double SoftenConstant, int K ){
	if(Temp <=0.0)
		throw domain_error("Total_Pressure_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw domain_error("Total_Pressure_PQCD expects positive"
		 " EnergyScaleConstant.");

	if(K<0 or K>6)
		throw invalid_argument("invalid K in Total_Pressure_PQCD");

	if(K==0 or K==1)
		return Total_Pressure_Free_QCD(Temp, mu_B);
	else{
		double Nf = 3.0, Nc = 3.0;

		double mu_q = mu_B / 3.0;           //quark chemical potential
		double mu_qh = mu_q/(2*pi*Temp);    //chem potential hat
		double mu_qh2 = pow(mu_qh, 2.0), mu_qh4 = pow(mu_qh, 4.0);

		double Lambda = EnergyScaleConstant*sqrt( Temp*Temp + pow(mu_q/pi,2.0) );

		// double CU_CM = .20;
		// double Lambda = EnergyScaleConstant*
		// 	pow( pow(Temp,4)+ pow(CU_CM*mu_q/pi,4.0) ,.25);

		double Lambdah = Lambda/(2.0*pi*Temp);  //hat

		double Pressure;

		double alphas = alphas_3loop(Lambda, SoftenConstant);

		Pressure = 1.0 + 21.0/32.0*Nf*(1.0 + 120.0/7.0*mu_qh2 + 240.0/7.0*mu_qh4);
		if(K>=2)
			Pressure += (-15.0/4.0*(1.0 + 5.0/12.0*Nf*(1.0 + 72.0/5.0*mu_qh2
				 + 144.0/5.0*mu_qh4)))*alphas/pi;

		if(K>=3)
			Pressure += (30.0 * pow((1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf), 1.5))*pow(alphas/pi, 1.5);

		if(K>=4)
			Pressure += ( 237.223 + (15.963 + 124.773*mu_qh2 - 319.849*mu_qh4)*Nf
				 -(0.415 + 15.926*mu_qh2 + 106.719*mu_qh4)*Nf*Nf
				 +135.0/2.0*(1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf )*
				 log(alphas/pi*(1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf))
				 -165.0/8.0*(1.0 + 5.0/12.0*
					(1.0 + 72.0/5.0*mu_qh2 + 144.0/5.0*mu_qh4)*Nf)*(1.0 - 2.0/33.0*Nf)
				 *log(Lambdah))*pow(alphas/pi, 2.0);

		if(K>=5)
			Pressure += ( -sqrt(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*(799.149 +
				 (21.963 - 136.33*mu_qh2 + 482.171*mu_qh4)*Nf
				 + (1.926 + 2.0749*mu_qh2 - 172.07*mu_qh4)*Nf*Nf )
				 + 495.0/2.0*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*(1.0 - 2.0/33.0*Nf)*
				 log(Lambdah))*pow(alphas/pi, 2.5);

		if(K>=6)
			Pressure += ( -(659.175 + (65.888 - 341.489*mu_qh2 + 1446.514*mu_qh4)*Nf
				 + (7.653 + 16.225*mu_qh2 - 516.210*mu_qh4)*Nf*Nf
				 - 1485.0/2.0*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*
				 (1.0 - 2.0/33.0*Nf)*log(Lambdah)
				 )*log(alphas/pi*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*4.0*pi*pi )
				 - 475.587*log(alphas*4.0*pi*Nc))*pow(alphas/pi, 3.0);


		Pressure=8.0*pi*pi/45.0*pow(Temp, 4.0)*Pressure;

		//inclusion of m_s
		// double m_s=100.0;
		// Pressure=Pressure-Nc*m_s*m_s*pow(Temp, 2.0)*(1/12.0 +mu_qh2);

		return Pressure;
	}
}


//finite difference:  q_B = dP/dmu_B
double Total_Baryon_Dens_PQCD(double Temp, double mu_B,
 double EnergyScaleConstant, double SoftenConstant){
	if(Temp <=0.0)
		throw domain_error("Total_Baryon_Dens_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw domain_error("Total_Baryon_Dens_PQCD expects positive"
		 " EnergyScaleConstant.");

	double dmu_B = 0.01;  //reasonable step size

	double Baryon = ( Total_Pressure_PQCD(Temp, mu_B + dmu_B,
		EnergyScaleConstant, SoftenConstant) -
	 Total_Pressure_PQCD(Temp, mu_B - dmu_B,
		EnergyScaleConstant, SoftenConstant) ) / (2.0*dmu_B);

	 return Baryon;
}

//Finite difference: s = dP/dT
double Total_Entropy_Dens_PQCD(double Temp, double mu_B,
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw domain_error("Total_Entropy_Dens_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw domain_error("Total_Entropy_Dens_PQCD expects positive"
		 " EnergyScaleConstant.");

	double T = Temp;
	double dT = 0.01;   //reasonable step size

	double Entropy = ( Total_Pressure_PQCD(T + dT, mu_B,
		EnergyScaleConstant, SoftenConstant) -
	 Total_Pressure_PQCD(T - dT, mu_B,
		EnergyScaleConstant, SoftenConstant) ) / (2.0*dT);

	 return Entropy;
}


double Total_Energy_Dens_PQCD(double Temp, double mu_B,
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw domain_error("Total_Energy_Dens_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw domain_error("Total_Energy_Dens_PQCD expects positive"
		 " EnergyScaleConstant.");

	double Press = Total_Pressure_PQCD(Temp, mu_B,
							EnergyScaleConstant, SoftenConstant);

	double Baryon = Total_Baryon_Dens_PQCD(Temp, mu_B,
							EnergyScaleConstant, SoftenConstant);

	double Entropy = Total_Entropy_Dens_PQCD(Temp, mu_B,
							EnergyScaleConstant, SoftenConstant);

	return (Temp*Entropy - Press + mu_B*Baryon); //from thermo identity
}

//----------------------------------------------------
//pQCD model struct and pybind things

pQCD::pQCD(double EScl, double Soft){
	this->EScl=EScl; this->Soft=Soft;
}

pQCD::pQCD(double EScl, double Soft, int K){
	this->EScl=EScl; this->Soft=Soft;
	this->K=K;
}

void pQCD::compute_single( double Temp, double mu_B,double * y){

	double Z=Total_Pressure_PQCD(Temp, mu_B,EScl,Soft,K);
	y[0]=Z;

	if(divs>0){
		double Z_p0,Z_0m,Z_0p,Z_m0;
		double dT=.1, dmu=.1;
		Z_p0=Total_Pressure_PQCD(Temp+dT, mu_B,EScl,Soft,K);
		Z_0m=Total_Pressure_PQCD(Temp, mu_B-dmu,EScl,Soft,K);
		Z_0p=Total_Pressure_PQCD(Temp, mu_B+dmu,EScl,Soft,K);
		Z_m0=Total_Pressure_PQCD(Temp-dT, mu_B,EScl,Soft,K);
		y[1]=(Z_p0-Z_m0)/(2*dT);
		if (mu_B==0) y[2]= 0.0;
		else y[2]=(Z_0p -Z_0m)/(2*dmu);

		if(divs>1) {
			double Z_pp,Z_mm;
			Z_pp=Total_Pressure_PQCD(Temp+dT, mu_B+dmu,EScl,Soft,K);
			Z_mm=Total_Pressure_PQCD(Temp-dT, mu_B-dmu,EScl,Soft,K);
			y[3]=(Z_p0-2*Z+Z_m0)/(dT*dT);
			if (mu_B==0) y[4]= 0.0;
			else y[4]=(Z_pp+Z_mm+2*Z-Z_0p-Z_0m-Z_p0-Z_m0)/(2*dT*dmu);
			y[5]=(Z_0p-2*Z+Z_0m)/(dmu*dmu);
		}
	}
}

//Used by pythons pickling and copying features
py::tuple pQCD::getstate(){
	return py::make_tuple(3,Therm::getstate(),Soft,EScl,K);
}

void pQCD::setstate(py::tuple t){
	if (t.size() != 5)
		throw std::runtime_error("Invalid state!(pQCD)");

	Therm::setstate(t[1]);
	Soft=t[2].cast<double>();
	EScl=t[3].cast<double>();
	K=t[4].cast<int>();
}

//used to include struct in pybind Module
void init_pQCD(py::module_ &m){
	py::class_<pQCD,Therm>(m, "pQCD").def(py::init<>(),"pQCD Object")
		.def(py::init<double, double>(),"pQCD Object")
		.def(py::init<double, double, int>(),"pQCD Object")
		.def("__repr__",
			[](const pQCD &a) {
				return "<pQCD: EScl = " + std::to_string(a.EScl) +
				", Soft = " + std::to_string(a.Soft) +
				", K = " + std::to_string(a.K) + ">";
			})
		.def(py::pickle(
			[](pQCD &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				pQCD p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		).def_readwrite("Soft", &pQCD::Soft).def_readwrite("EScl", &pQCD::EScl)
		.def_readwrite("K", &pQCD::K);
}