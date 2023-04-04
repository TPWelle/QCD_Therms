
#include <cmath>
#include "ColdQCD.hpp"

using namespace std;


//Low-T quark gas from Gerhold, Ipp, Rebhan 2004 in PRD

double alphas_3loop(double Lambda, double SoftenConstant, int Nf)
{
	double LambdaMS = 290.0;   //MeV

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

//Ancillary function used by Pressure
double F_bar(double mu2, int Nf){
	if(Nf==1) return 0;
	double P=2*Nf*Nf*mu2*mu2*log(Nf);
	P=P+(8.0/3.0)*Nf*(Nf-1)*mu2*mu2*L2;
	return P;
}


double Press4_T0(double mu2, int Nf, int Nc){

	int Ng= Nc*Nc-1;

	double P=0.0;

	double M2=mu2*Nf;

	double soft=2.5;

	//alpha_s / (4*pi)
	double al_4pi=alphas_3loop(sqrt(M2),soft,Nf)/(4*pi);

	double mux=(1/3.0)*(11*Nc -2*Nf)*log(mu2/M2)-2.25*Nc +.409*Nf
		-3.697 -4.236/Nc;

	mux=Nf*mu2*mu2*(Nc/3.0 +Ng*al_4pi*(-1.0 + al_4pi*mux));
	P=mux-Ng*pow(al_4pi,2) *(mu2*mu2*(2*log(al_4pi)-.476)+F_bar(mu2, Nf));

	return P/(4*pi*pi);

}


double Pressure(double Temp, double mu_B, int Nf, int Nc){

	int Ng= Nc*Nc-1;
	double T2=Temp*Temp, mu2= mu_B*mu_B/9.0;
	double EnergyScaleConstant=2*pi*1.5;
	double Lambda=EnergyScaleConstant*sqrt( T2 + mu2/(pi*pi) );

	double P_gas=Ng*pi*pi*T2*T2/45.0 +Nc*Nf*(
		7*pi*pi*T2*T2/180.0 +T2*mu2/6.0 +mu2*mu2/(12.0*pi*pi));

	double soft=2.5;

	double g2=2*pi*Nf*alphas_3loop(Lambda,soft,Nf);
	double Lg2=log(g2*mu2/T2);
	double P=0.0;
	P=mu2*T2*g2/pow((12*pi),2)*(Lg2-2.5119776);
	P=P-.00887493874*pow(g2*mu2*T2*T2,2.0/3.0);
	P=P+.01444672361792*pow(g2*mu2,1.0/3.0)*pow(T2,5.0/3.0);
	P=Ng*(P-.025827014176*T2*T2*(Lg2+8.198697024078));
	double P4=Press4_T0(mu2,Nf,Nc);

	// if(!isfinite(P+P_gas+P4)){
	// 	cout<<"!"<<Temp<<" | "<<mu_B<<endl;
	// 	cout<<sqrt(g2)<<endl;
	// }
	return P+P_gas+P4;
}


void ColdQCD::compute_single( double Temp, double mu_B, double * y){

	double Z=Pressure(Temp, mu_B, Nf, Nc);
	y[0]=Z;

	double Z_p0,Z_0m,Z_0p,Z_m0, Z_pp,Z_mm;
	double dT=.1, dmu=.1;
	Z_p0=Pressure(Temp+dT, mu_B, Nf, Nc); Z_m0=Pressure(Temp-dT, mu_B, Nf, Nc);
	Z_0p=Pressure(Temp, mu_B+dmu, Nf, Nc); Z_0m=Pressure(Temp, mu_B-dmu, Nf, Nc);
	Z_pp=Pressure(Temp+dT, mu_B+dmu, Nf, Nc); Z_mm=Pressure(Temp-dT, mu_B-dmu, Nf, Nc);

	y[1]=(Z_p0-Z_m0)/(2*dT);
	if (mu_B==0) {
		y[2]= 0.0; y[4]= 0.0;
	} else {
		y[2]=(Z_0p -Z_0m)/(2*dmu);
		y[4]=(Z_pp+Z_mm+2*Z-Z_0p-Z_0m-Z_p0-Z_m0)/(2*dT*dmu);
	}
	y[3]=(Z_p0-2*Z+Z_m0)/(dT*dT);
	y[5]=(Z_0p-2*Z+Z_0m)/(dmu*dmu);

}


py::tuple ColdQCD::getstate(){
	return py::make_tuple(5,Nf,Nc,Therm::getstate());
}

void ColdQCD::setstate(py::tuple t){
	if (t.size() != 4)
		throw std::runtime_error("Invalid state!(ColdQCD)");

	Therm::setstate(t[1]);
	Nf=t[2].cast<double>();
	Nc=t[3].cast<double>();
}

void init_ColdQCD(py::module_ &m){
	py::class_<ColdQCD,Therm>(m, "ColdQCD").def(py::init<>(),"ColdQCD model")
		.def_readwrite("Nf",&ColdQCD::Nf)
		.def_readwrite("Nc",&ColdQCD::Nc)
		.def("__repr__", [](const ColdQCD &a) {
				return "<ColdQCD>";
		}).def(py::pickle([](ColdQCD &p){
			return p.getstate();
		}, [](py::tuple t){ // __setstate__
			ColdQCD p; /* Create a new C++ instance */
			p.setstate(t); /* Assign additional state */
			return p;
		})
	);
}

