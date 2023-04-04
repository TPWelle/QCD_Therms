/*
	This file computes the EOS for the lattice parameterization used in:
	"Lattice-based equation of state at finite baryon number, electric charge and
	strangeness chemical potential" P. Parotto et. al., Jan 2020
*/
#include <cmath>
#include <algorithm>
#include <iostream>
// #include <vector>
#include <array>
#include "LatPar.hpp"

using namespace std;

//Specifies which susceptibility (X_BQS) is being obtained
enum Coef_code {c000, c200,c020,c002, c110,c101,c011, c400,c040,c004,
	c310, c301,c031, c130,c103,c013, c220,c202,c022, c211,c121,c112};

// Computes the numerator(A) and the denominator(B), as a function of u=1/t (t=154/Temp), 
// for the parameterization of the susceptibility specified by 'code'.
// Also gets derivatives dA, ddA, dB, ddB.
void poly_sum(double u, Coef_code code, double &A, double &dA, double &ddA, double &B, double &dB, double &ddB){

	array<double, 10> as; array<double, 10> bs;

	switch(code){

		case c000:
			as={7.53891,-6.18858,-5.37961,7.0875,-0.97797,0.0302636,0,0,0,0};
			bs={ 2.2453,-6.02568,15.3737,-19.6331,10.24,0.799479,0,0,0,0};
			break;
		case c200:
			throw invalid_argument("In poly_sum, c200 requires exceptional function.");
			break;
		case c020:
			as={-1.254,13.7781,-20.8361,11.4637,-1.52145,0.0563044,0,0,0,0};
			bs={-2.08695,22.3712,-33.4035,19.9497,-6.67937,4.1127,0,0,0,0};
			break;
		case c002:
			as={0.728917,-1.73212,1.61219,-0.706361,0.192223,-0.0164219,-0.0040308,0.00044212,0,0};
			bs={0.634185,-0.484646,-3.02879,7.29526,-5.94029,0.954829,0.782178,0.0848009,0,0};
			break;
		case c110:
			as={0.611997,0.260951,0.439882,4.04624,-0.492197,-0.479177,0.108023,-0.000271088,0,0};
			bs={506.969,2.07112,-1310.73,-47.3907,1855.62,207.417,-2635.03,1616.16,0,0};
			break;
		case c101:
			as={-3.42744,0.0807472,0.155933,1.76331,-0.350538,-0.547143,0.0641196,-0.000271926,0,0};
			bs={8.81578,4.53879,70.1272,-212.977,-287.925,1688.03,-2130.95,901.004,0,0};
			break;
		case c011:
			as={0.975914,-2.2118,1.99441,-0.710665,0.100002,-0.0043751,0,0,0,0};
			bs={ 2.94254,-5.97226,4.37484,0.723152,-4.3139,3.70245,0,0,0,0};
			break;
		case c400:
			as={0.0697892,-0.0759267,0.0270699,-0.00183789,-0.00102026,0.000248834,-0.0000205803,5.78113e-7,0,0};
			bs={3.3139,-2.34182,-3.05239,0.281088,3.36387,-1.47861,0.232943,-0.00920141,0,0};
			break;
		case c040:
			as={0.519384,-2.61484,6.99796,-9.37407,5.50677,-0.933273,0.0628049,-0.00149075,0,0};
			bs={2.78757,-7.70015,7.34828,5.60254,-16.5647,10.1847,-1.46422,0.258243,0,0};
			break;
		case c004:
			as={3.99178,-10.8564,11.4807,-5.56961,1.43254,-0.204083,0.0152834,-0.00047076,0,0};
			bs={7.26105,-25.0961,41.2002,-31.1539,-3.87268,19.7369,-9.31673,2.02404,0,0};
			break;
		case c310:
			as={0.000214078,-0.00277202,0.0107602,-0.0189801,0.0163346,-0.00649086,0.00102683,-0.0000118454,0,0};
			bs={0.628355,-1.27107,-0.0555062,0.801392,0.649844,-0.248501,-1.16057,0.662302,0,0};
			break;
		case c301:
			as={-0.606637,0.940635,-0.609091,0.211817,-0.0423212,0.0048043,-0.000283315,6.59604e-6,0,0};
			bs={22.8266,-19.1507,-33.6479,25.4636,17.3853,-0.671223,-19.7378,9.96533,0,0};
			break;
		case c031:
			as={1.39052,-2.95215,2.99901,-1.3976,0.337495,-0.0441243,0.0029685,-0.0000804859,0,0};
			bs={52.129,-92.6007,24.1788,32.9419,-12.5404,-1.67767,1.02439,0.502227,0,0};
			break;
		case c130:
			as={1.33817,-0.36966,-7.73766,12.6268,-7.54688,2.27058,-0.380023,0.0357606,-0.00175991,0.0000349795};
			bs={32.3922,-36.2407,-44.2609,31.2543,50.794,17.5211,-7.80941,-13.3867,-118.309,93.7845};
			break;
		case c103:
			as={-0.0853497,0.09878,-0.0477156,0.0124373,-0.00188339,0.000165099,-7.72499e-6,1.47927e-7,0,0};
			bs={0.285383,0.769297,-3.15803,1.59797,3.54785,-0.652119,-6.48277,4.28691,0,0};
			break;
		case c013:
			as={0.23137,-0.607108,0.574083,-0.232842,0.0476026,-0.00514917,0.000279883,-5.97476e-6,0,0};
			bs={1.12154,-2.86563,2.35378,-0.14257,-0.827056,0.35061,-0.0544297,0.125906,0,0};
			break;
		case c220:
			as={0.131897,-0.151923,0.0728375,-0.0188047,0.00281673,-0.000244096,0.0000112936,-2.14344e-7,0,0};
			bs={2.46229,-1.78965,3.86743,-3.007,-4.28013,0.190242,3.36159,-0.215634,0,0};
			break;
		case c202: //!!!!!!!!
			as={0.0481773,-0.0633491,0.034631,-0.0101557,0.00171648,-0.000166247,8.49007e-6,-1.75132e-7,0,0};
			bs={.505109,0.555159,-2.50987,0.346874,2.47285,0.611415,-3.84829,2.02716,0,0};
			break;
		case c022:
			as={1.03006,-2.50946,2.44698,-1.00851,0.207566,-0.0225078,0.00122422,-0.000026132,0,0};
			bs={15.1999,-40.1845,44.1416,-19.6254,-13.5991,25.2683,-12.6079,2.72985,0,0};
			break;
		case c211:
			as={.146608,-0.533936,0.834892,-0.645642,0.260112,-0.0540238,0.00527256,-0.000182156,0,0};
			bs={5.80204,-15.5399,5.25306,18.444,-1.81185,-20.1787,-4.61059,13.9429,0,0};
			break;
		case c121:
			as={-1.27191,2.11351,-1.4598,0.538497,-0.113414,0.0134801,-0.000826188,0.0000197941,0,0};
			bs={56.5761,-106.452,123.146,-162.408,94.5282,51.273,-77.7255,29.3669,0,0};
			break;
		case c112:
			as={-2.61752,4.37997,-3.00258,1.08514,-0.221478,0.0253053,-0.00148533,0.0000342702,0,0};
			bs={43.2755,-108.526,180.836,-134.256,-38.6051,46.669,6.94258,17.7581,0,0};
			break;
	}

	A=as[0] + u*as[1]; B=bs[0]+ u*bs[1];
	dA=as[1]; dB=bs[1];

	for(int k=2; k<10; k++){

		if(as[k]==0) break;

		A=A+as[k]*pow(u,k);
		dA=dA+k*as[k]*pow(u,k-1);
		ddA=ddA+k*(k-1)*as[k]*pow(u,k-2);

		B=B+bs[k]*pow(u,k);
		dB=dB+k*bs[k]*pow(u,k-1);
		ddB=ddB+k*(k-1)*bs[k]*pow(u,k-2);
	}
}

// Computes susceptibility X_B^2 and first 2 T derivatives. 
// Uses a different functional form than other susceptabilities.
void XB2(double u, double &X, double &dX, double &ddX){
	// double h1=-0.325372, h2=0.497729, f3=0.148987, f4=6.66388, f5=-5.07725;
	// double t=Temp/200.0;

	double h1=-0.422561039, h2=0.839482206, f3=0.148987, f4=5.1311876, f5=-5.07725;
	// double u=154.0/Temp;
	double th=tanh(f4/u + f5);

	X= f3*exp(-h1*u -h2*u*u)*(1.0 +th);
	dX=-( h1+2*h2*u + (f4/(u*u))*(1-th))*X;
	ddX=(pow( h1+2*h2*u + (f4/(u*u))*(1-th), 2) -2*h2
		+2*f4*(1-th)/pow(u,3)+f4*f4*(1-th*th)/pow(u,4))*X;

}

//Computes susceptibility (X_BQS) specified by 'code' and first 2 T derivatives.
void X_func(double u, Coef_code code , double &X, double &dX, double &ddX){

	if(code==c200){
		XB2(u,X,dX,ddX);
		return;
	}

	double A,dA,ddA,B,dB,ddB;
	poly_sum(u,code,A,dA,ddA,B,dB,ddB);

	X=A/B;
	dX=dA/B - A*dB/(B*B);
	ddX=ddA/B - ( 2*dA*dB + A*(ddB - 2*dB*dB/B))/(B*B);

	if(code==c002)X=X+.00083;
	else if(code==c011)X=X+.00012;
	else if(code==c002)X=X-00007;

}

int factorial(int n){
	if(n<0) throw invalid_argument("argument of factorial can't be negative");
	if(n==0 or n==1) return 1;
	return n*factorial(n-1);
}

//----------------------------------------------------
//LatPar model struct and pybind things

void LatPar::compute_single( double Temp, double mu_B, double * y){

	double P=0.0,s=0.0,n=0.0;
	double X_TT=0.0, X_Tm=0.0, X_mm=0.0;

	double u=154.0/Temp;
	// double vB=mu_B/154.0;
	double vB=mu_B/Temp;

	Coef_code codes[3]={c000,c200,c400};

	for(int k=0; k<3; k++){
		double X,dX,ddX;
		X_func(u,codes[k],X,dX,ddX );
		P=P+X*pow(vB,2*k)/factorial(2*k);
		if(k>0) n=n+X*pow(vB,2*k-1)/factorial(2*k-1);
		s=s+(pow(vB,2*k)/factorial(2*k))*((4-2*k)*X -u*dX);
		if( divs>1 ) {
			if(k>0) X_mm=X_mm+X*pow(vB,2*k-2)/factorial(2*k-2);
			if(k>0) X_Tm=X_Tm+(pow(vB,2*k-1)/factorial(2*k-1))*((4-2*k)*X -u*dX);
			X_TT=X_TT+(pow(vB,2*k)/factorial(2*k))*( (k-4)*(k-3)*X +(2*k-6)*u*dX +u*u*ddX);
		}
	}

	y[0]=P*pow(Temp,4);
	y[1]=s*pow(Temp,3);
	y[2]=n*pow(Temp,3);
	y[3]=X_TT*pow(Temp,2);
	y[4]=X_Tm*pow(Temp,2);
	y[5]=X_mm*pow(Temp,2);

}

//Used by pythons pickling and copying features
py::tuple LatPar::getstate(){
	return py::make_tuple(4,Therm::getstate());
}

void LatPar::setstate(py::tuple t){
	if (t.size() != 2)
		throw std::runtime_error("Invalid state!(LatPar)");
	Therm::setstate(t[1]);
}

//used to include struct in pybind Module
void init_LatPar(py::module_ &m){
	py::class_<LatPar,Therm>(m, "LatPar").def(py::init<>(),"LatPar model")
		.def("__repr__", [](const LatPar &a) {
				return "<LatPar>";
		}).def(py::pickle(
		[](LatPar &p){
			return p.getstate();
		}, [](py::tuple t){ // __setstate__
			LatPar p; /* Create a new C++ instance */
			p.setstate(t); /* Assign additional state */
			return p;
		})
	);
}

