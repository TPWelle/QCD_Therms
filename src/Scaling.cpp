#include <iostream>
#include <stdexcept>

#include "Scaling.hpp"
#include <omp.h>
#include <cmath>

using namespace std;




//----------------------------------------------------
//struct for Scaling EOS model

Scaling::Scaling(Therm * pP):
pP(pP){}

Scaling::Scaling(Therm * pP, double Tc, double muc):
pP(pP){
	this->Tc=Tc; this->muc=muc;
}

void Scaling::set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs){

	Therm::set_Tmu( Temps, mu_Bs);
	clear_comp_flag();
	pP->divs=2;//Constituent model needs two derivatives
	pP->N_threads=N_threads;
	pP->set_Tmu( Temps, mu_Bs);
	//Allocate Rs and thetas
	py::buffer_info buf=Temps.request();
	thetas = py::array_t<double>(buf.size); thetas.resize(buf.shape);
	Rs = py::array_t<double>(buf.size); Rs.resize(buf.shape);
	Wfs = py::array_t<double>(buf.size); Wfs.resize(buf.shape);
	Wf_Ts = py::array_t<double>(buf.size); Wf_Ts.resize(buf.shape);
	Wf_ms = py::array_t<double>(buf.size); Wf_ms.resize(buf.shape);

}

void Scaling::clear_comp_flag(){
	this->comp_flag=false;
	pP->clear_comp_flag();
}

// set values for parameters h3 and h5, and compute associated
// values theta0, and g0,...,g3
void Scaling::set_hs(double h3, double h5){
	this->h3=h3; this->h5=h5;
	if(h5>=0.0) theta0=sqrt((-h3 -sqrt(h3*h3 -4*h5))/(2*h5));
	else theta0=sqrt((-h3 +sqrt(h3*h3 -4*h5))/(2*h5));
	g0=(1+h3+h5)/(delta+1.0);
	g1=(1+h3+h5 -2*beta*(1+2*h3+3*h5))/(2.0-2*alpha);
	g2=(.5*h3+h5 -beta*(h3+3*h5))/alpha;
	g3=h5*(beta-.5)/((1+alpha));
}

double Scaling::n_BG(double Temp, double mu_B){
	double yS[6];
	pP->compute_single(Temp,mu_B,yS);
	return yS[2];
}

//h and g functions used in Scaling EOS and derivatives.
double Scaling::h_func(double theta){
	double X=theta*theta;
	return theta*(1+h3*X +h5*X*X);
}
double Scaling::dh_func(double theta){
	double X=theta*theta;
	return 1+3.0*h3*X +5.0*h5*X*X;
}
double Scaling::g_func(double theta){
	double X=(1.0-pow(theta,2));
	return g0+g1*X +g2*X*X +g3*pow(X,3);
}
double Scaling::dg_func(double theta){
	double X=(1.0-pow(theta,2));
	return -2*theta*(g1 +2.0*g2*X +3*g3*X*X);
}
//Used for root-finding to obtain theta
double Scaling::obj_function(double theta){
	return h_func(theta)/pow(abs(1-theta*theta),beta*delta);
}


// compute theta numerically from x=(mu-mux)/muc and y=T/Tc -1
double Scaling::get_theta(double x, double y){

	double theta1, theta2, F1, F2;
	double Z=abs(x)/h0 /pow(abs(y),beta*delta);

	// Set up left and right bounds for root finding
	if(y==0.0){
		if(x>0.0) return 1.0;
		else if(x<0.0) return -1.0;
		else return 0.0;
	}else if(y>0.0){
		if(x==0) return 0.0;
		else {
			theta1=0.0; theta2=.999; F1=-Z;//slightly less
			do {
				theta2=.9+theta2/10.0;
				F2=obj_function(theta2)-Z;
			} while (F2<0.0);
		}
	}else if(y<0.0){
		if(x==0) return theta0;
		else {
			theta2=1.0001; theta1=theta0; F1=-Z;//slightly less
			do {
				theta2=(theta2-1.0)/10.0 +1.0;
				F2=obj_function(theta2)-Z;
			} while (F2<0.0);
		}
	}
	//Use root finding to get theta, using both secant
	//and bisection steps (works the best)
	double theta, F=1.0,tol=1e-8;
	int iter=0;
	while(abs(F)>tol and iter<100){
		iter++;
		if(iter%3==0){ //bisection step
			theta=.5*(theta1+theta2);
			F=obj_function(theta)-Z;
		}else{ //secant step
			theta=(F2*theta1-F1*theta2)/(F2-F1);
			F=obj_function(theta)-Z;
		}
		if(F <= 0.0){
			theta1=theta;
			F1=F;
		} else {
			theta2=theta;
			F2=F;
		}
	}
	return copysign(theta,x);
}

void Scaling::get_thetas(){

	get_muxs();
	int N=Temps.request().size;

	double *pthetas = get_arr_pointer(thetas);
	double *pRs = get_arr_pointer(Rs);
	double *pTemp= get_arr_pointer(Temps);
	double *pmu= get_arr_pointer(mu_Bs);

	// int N_threads=7;
	if(mu_code==0 or mu_code==1){
		#pragma omp parallel for schedule(dynamic) \
		num_threads(N_threads)
		for(int i=0;i<N;i++){

			double Temp=pTemp[i],mu_B=pmu[i];

			pthetas[i]=get_theta((mu_B-mux.at(Temp))/muc,Temp/Tc -1.0);
			if(Temp==Tc)
				pRs[i]=pow(abs( (mu_B/muc -1.0)/(h0*h_func(1.0))) ,1.0/(beta*delta));
			else
				pRs[i]=abs((Temp/Tc -1.0)/(1.0-pthetas[i]*pthetas[i]));

		}
	} else if (mu_code==2){
		double rho=2.0, w=1.0, alp1=.0672;
		double ca1=cos(alp1), sa1=sin(alp1);
		double ca2=-sa1, sa2=ca1;//alp2== pi/2 + alp1
		#pragma omp parallel for schedule(dynamic) \
		shared(ca1,ca2,sa1,sa2,rho,w) num_threads(N_threads)
		for(int i=0;i<N;i++){

			double tt=pTemp[i]/Tc -1.0;
			// double mm=pmu[i]/muc - 1.0;
			double mm=(pmu[i]-muc)/Tc ;
			double tau=-(ca2*tt +sa2*mm)/(rho*w);
			double zeta=(ca1*tt +sa1*mm)/w;

			pthetas[i]=get_theta(zeta,tau);
			if(tau==0.0)
				pRs[i]=pow(abs( zeta/(h0*h_func(1.0))) ,1.0/(beta*delta));
			else
				pRs[i]=abs(tau/(1.0-pthetas[i]*pthetas[i]));

		}
	}
}

inline void Scaling::mux_insert(double Temp,
  double mu, double dmu, double ddmu){
    mux.insert(pair<double, double>(Temp,mu));
    dmux.insert(pair<double, double>(Temp,dmu));
    ddmux.insert(pair<double, double>(Temp,ddmu));
}

//Find mu_x(T), the chemical potential where n_B=nc
void Scaling::get_mux(double Temp){
	if(mux.count(Temp)) return;

    if(mu_code==1){
        mux_insert(Temp, muc, 0.0, 0.0);
        return;
    }else if(mu_code==2){
        double alp1=.0672; double ca1=1.0/tan(alp1);
        double mu=muc-ca1*(Temp-Tc);
        double dmu=-ca1;
        if(mu<=0) mux_insert(Temp, muc, 0.0, 0.0);
        else mux_insert(Temp, mu, dmu, 0.0);
        return;
    }

	//target pressure for n_BG
	double n_t=nc-n0;
	double mu1,mu2,mu,mup,mum,n1,n2,dT=.01;
	//initializing
	if (Temp>Tc){
		double y[6];
		pP->compute_single(Temp,0,y);
		mu1=n_t/y[5]; n1=n_BG(Temp,mu1);
		mu2=muc; n2=n_t;
		if(n1>=n_t){
			mu2=mu1; n2=n1;
			mu1=0.0; n1=0.0;
		}
	} else {
		mu1=muc; mu2=1.2*muc;
		n1=n_BG(Temp,mu1); n2=n_BG(Temp,mu2);
		while(n2<n_t){
			mu1=mu2; n1=n2;
			mu2=1.5*mu2; n2=n_BG(Temp,mu2);
		}
	}
	double nB=n_t+1.0,tol=1e-4;
	int iter=0;
	//compute mux
	while(abs(nB-n_t)>tol and abs(mu1-mu2)>tol and iter<80){
		iter++;
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		//Swapping between secant and bisection method gives best results
		if(iter%2==0){ //bisection step
			mu=.5*(mu1+mu2);
		} else { // Secant Step
			mu=((n2-n_t)*mu1-(n1-n_t)*mu2)/(n2-n1);
		}
		nB=n_BG(Temp,mu);
		if(nB <= n_t){
			mu1=mu; n1=nB;
		} else {
			mu2=mu; n2=nB;
		}
	}
	//compute mup
	iter=0; nB=n_t+1.0;
	mu2=mu; n2=n_BG(Temp+dT,mu2);
	mu1=.95*mu2; n1=n_BG(Temp+dT,mu1);

	while(n1>n_t){
		mu2=mu1; n2=n1;
		mu1=.9*mu2;
		n1=n_BG(Temp+dT,mu1);
	}
	while(abs(nB-n_t)>tol  and abs(mu1-mu2)>tol and iter<80){
		iter++;
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		if(iter%2==0){ //bisection step
			mup=.5*(mu1+mu2);
		} else { // Secant Step
			mup=((n2-n_t)*mu1-(n1-n_t)*mu2)/(n2-n1);
		}
		nB=n_BG(Temp+dT,mup);
		if(nB <= n_t){
			mu1=mup; n1=nB;
		} else {
			mu2=mup; n2=nB;
		}
	}
	//compute mum
	iter=0; nB=n_t+1.0;
	mu1=mu; n1=n_BG(Temp-dT,mu1);
	mu2=1.1*mu1; n2=n_BG(Temp-dT,mu2);
	while(n2<n_t){
		mu1=mu2; n1=n2;
		mu2=1.05*mu1;
		n2=n_BG(Temp-dT,mu2);
	}
	while(abs(nB-n_t)>tol and abs(mu1-mu2)>tol and iter<80){
		iter++;
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		if(iter%2==0){ //bisection step
			mum=.5*(mu1+mu2);
		} else { // Secant Step
			mum=((n2-n_t)*mu1-(n1-n_t)*mu2)/(n2-n1);
		}
		nB=n_BG(Temp-dT,mum);
		if(nB <= n_t){
			mu1=mum; n1=nB;
		} else {
			mu2=mum; n2=nB;
		}
	}
    mux_insert(Temp, mu, (mup-mum)/(2*dT), (mup-2*mu +mum)/(dT*dT));

}


void Scaling::get_muxs(){

	int N=Temps.request().size;
	double * pTemp = get_arr_pointer( Temps );

	double yS[6];
	pP->compute_single(Tc,muc,yS);
	nc= yS[2]+n0; Pc= yS[0]+P0;//Critical Baryon Density and Pressure

	get_mux(Tc); //get mu_x for Tc

	#pragma omp parallel for schedule(dynamic) \
	num_threads(N_threads)
	for(int i=0;i<N;i++){
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		bool prev_T=false;

		//Avoid calculating mu_x multiple times if Temps are the same
		for(int j=0; j<i;j++){
			prev_T = prev_T or (pTemp[j]==pTemp[i]);
		}
		if(not prev_T)
			get_mux(pTemp[i]);

	}
}

void Scaling::compute(){

	if(comp_flag) return;

	get_thetas();

	pP->compute();

	int N=Temps.request().size;

	double *pTemp= get_arr_pointer(Temps), *pmu= get_arr_pointer(mu_Bs);
	double *pthetas= get_arr_pointer(thetas),*pRs= get_arr_pointer(Rs);

	double *Z= get_arr_pointer( Y ), *Z_T= get_arr_pointer( Y_T ), *Z_m= get_arr_pointer( Y_m );
	double *Z_TT= get_arr_pointer( Y_TT ), *Z_Tm= get_arr_pointer( Y_Tm ),*Z_mm= get_arr_pointer( Y_mm );

	double *P= get_arr_pointer( pP->Y ), *P_T= get_arr_pointer( pP->Y_T ), *P_m= get_arr_pointer( pP->Y_m );
	double *P_TT= get_arr_pointer( pP->Y_TT ), *P_Tm= get_arr_pointer( pP->Y_Tm ),*P_mm= get_arr_pointer( pP->Y_mm );

	double *pWf=get_arr_pointer(Wfs), *pWf_T=get_arr_pointer(Wf_Ts),*pWf_m=get_arr_pointer(Wf_ms);


	#pragma omp parallel for schedule(dynamic) \
	shared(divs) num_threads(N_threads)
	for(int i=0;i<N;i++){

		//Display Progress
		if(omp_get_thread_num()%N_threads==0 and display_progress){
			cout<< i<<"/"<<N<<endl;
		}
		//Allow Ctrl-C abort
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();

		double theta=pthetas[i],th2=theta*theta, R=pRs[i];

		double nn=P_m[i]/(nc-n0); //Background density normalized
		double mux1=mux.at(pTemp[i]), dmux1=dmux.at(pTemp[i]);

		//
		int W_code=1;
		double Wf,dWfdmu,dWfdT,dWfdmumu,dWfdTmu,dWfdTT;
		if(W_code==0){//Window Function///////////////
			Wf=exp(-pow( (nn*nn-1.0)/((c_star/(1.0-n0/nc))*nn) ,2));
			double dWf=(2.0/pow(c_star/(1.0-n0/nc),2))*(1.0-pow(nn,4))*Wf/(pow(nn,3)*nc);
			//Will be used when div=2 is included. Currently incomplete.
			// double ddWf=(2.0/pow(c_star/(1.0-n0/nc),2))*(1.0-pow(nn,4))*dWf/(pow(nn,3)*nc)
			// 	-(2.0/pow(nc*c_star/(1.0-n0/nc),2))*(3.0*pow(nn,-4)-1.0)*Wf;
			dWfdT=dWf*P_Tm[i]; dWfdmu=dWf*P_mm[i];

		}else if(W_code==1){ //Window Function with mus///////////////


			// //j=1
			// Wf=exp(-pow( (pmu[i]*pmu[i]-mux1*mux1)/(c_star*muc*pmu[i]) ,2));
			// dWfdmu=(2.0/pow(c_star*muc,2))*(pow(mux1,4)-pow(pmu[i],4))*Wf/pow(pmu[i],3);
			// dWfdT=(4.0*mux1*dmux1/pow(c_star*muc*pmu[i],2))*(pow(pmu[i],2)-pow(mux1,2))*Wf;

			//j>1
			int J=2;
			Wf=exp(-pow( (pow(pmu[i],2*J)-pow(mux1,2*J))/(c_star*pow(muc*pmu[i],J)) ,2));
			dWfdmu=(2.0/(pow(muc,2*J))*c_star*c_star)*(pow(mux1,4)-pow(pmu[i],4))/pow(pmu[i],2*J+1);
			dWfdT= (2*J*dmux1*Wf/(c_star*c_star))*(pow(mux1*pmu[i],2*J)-pow(mux1,4*J))
				/(pow(pmu[i],2*J+1)*pow(muc,2*J));


			//Will be used when div=2 is included. Currently incomplete.
			// dWfdmumu=(-2.0/pow(c_star*muc,2))*( (3*pow(mux1,4)+pow(pmu[i],4))*Wf/pow(pmu[i],4)
			// 	(pow(pmu[i],4)-pow(mux1,4))*dWfdmu/pow(pmu[i],3));
			// dWfdTmu=(-2.0/pow(c_star*muc,2))*( (pow(mux1,4)-pow(pmu[i],4))*dWfdT/pow(pmu[i],3)
			// 	-(4*dmux1*pow(mux1,3))*Wf/pow(pmu[i],3));
			// dWfdTT=(4.0/pow(c_star*muc*pmu[i],,2))*(
			// 	(ddmux1*mux1*(pmu[i]*pmu[i]-mux1*mux1) +dmux1*dmux1*(pmu[i]*pmu[i]-3*mux1*mux1) )*Wf
			// 	+mux1*dmux1*(pmu[i]*pmu[i]-mux1*mux1)*dWfdT );
		}

		double P1;

		if(mu_code==0 or mu_code==1){
			P1=P0 + h0*muc*n0*pow(R,beta*delta)*h_func(theta)
				+m0*h0*muc*n0*pow(R,beta*(delta+1))*(theta*h_func(theta) -g_func(theta));
		}else if(mu_code==2){
			P1=m0*h0*muc*n0*pow(R,beta*(delta+1))*(theta*h_func(theta) -g_func(theta));
		}

		if (nn==1.0 or mu_code>0){
			Wf=1.0; dWfdT=0.0; dWfdmu=0.0;
		}else if(pmu[i]==0.0){
			Wf=0.0; dWfdT=0.0; dWfdmu=0.0;
		}

		pWf[i]=Wf; pWf_T[i]=dWfdT; pWf_m[i]=dWfdmu;
		if(mu_code==0) Z[i]=P[i]+Wf*P1;
		else if(mu_code==1) Z[i]=P1;
		else if(mu_code==2) Z[i]=P1;

		if(divs>0){
			double s1,n1;
				double g_tilde=(beta*theta*h_func(theta)-(2-alpha)*g_func(theta))/(1-th2);

			if(mu_code==0){
				n1=n0+m0*n0*theta*pow(R,beta);
				s1=(m0*h0*muc*n0/Tc)*pow(R,1-alpha)*g_tilde - n1*dmux1;
				if (pmu[i]==0) Z_m[i]= 0.0;
				else Z_m[i]=P_m[i] + n1*Wf + P1*dWfdmu;
				if(theta==0) Z_T[i]=P_T[i] -s1;
				else Z_T[i]=P_T[i] - s1*Wf + P1*dWfdT;
			}else if(mu_code==1){
				n1=n0+m0*n0*theta*pow(R,beta);
				s1=(m0*h0*muc*n0/Tc)*pow(R,1-alpha)*g_tilde;
				if (pmu[i]==0) Z_m[i]= 0.0;
				else Z_m[i]=n1;
				Z_T[i]=s1;
			}else if(mu_code==2){

				double rho=2.0, w=1.0, alp1=.0672;
				double ca1=cos(alp1), sa1=sin(alp1);
				double ca2=-sa1, sa2=ca1;//alp2== pi/2 + alp1

				n1=n0*m0*(muc/Tc)*(sa1*rho*theta*pow(R,beta) -sa2*h0*g_tilde*pow(R,1-alpha))/w;
				s1=n0*m0*(muc/Tc)*(ca1*rho*theta*pow(R,beta) -ca2*h0*g_tilde*pow(R,1-alpha))/w;
				if (pmu[i]==0) Z_m[i]= 0.0;
				else Z_m[i]=n1;
				Z_T[i]=s1;
			}
			//if(divs>1) see code at the bottom for a starting point on this
		}
	}
}

//Used to set n_0 and P_0 to be a particular fraction of n_c and P_c
void Scaling::set_n0P0(double n_rat, double P_rat){

	double y[6];
	pP->compute_single(Tc,muc,y);
	Pc=y[0]/(1.0-P_rat); P0=P_rat*Pc;
	nc=y[2]/(1.0-n_rat); n0=n_rat*nc;
}

//Used by pythons pickling and copying features
py::tuple Scaling::getstate(){
	return py::make_tuple(11,Therm::getstate(),pP->getstate()
		,alpha,beta,gamma,delta,Tc,muc,nc,Pc,n0,P0,m0,h0,h3,h5
		,c_star,mux,dmux,ddmux,thetas,Rs,Wfs,Wf_Ts,Wf_ms);
}

void Scaling::setstate(py::tuple t){
	if (t.size() != 26)
		throw std::runtime_error("Invalid state! (Scaling)");

	Therm::setstate(t[1]);
	switch(t[2].cast<py::tuple>()[0].cast<int>()){
		case 0:
			pP= new PT; break;
		case 1:
			pP= new EXI; break;
		case 2:
			pP= new EXII; break;
		case 3:
			pP= new pQCD; break;
		case 4:
			pP= new LatPar; break;
		case 10:
			pP= new Crossover; break;
	}
	pP->setstate(t[2]);
	// pWF->pP = this;
	alpha=t[3].cast<double>(); beta=t[4].cast<double>();
	gamma=t[5].cast<double>(); delta=t[6].cast<double>();

	Tc=t[7].cast<double>(); muc=t[8].cast<double>();
	nc=t[9].cast<double>(); Pc=t[10].cast<double>();
	n0=t[11].cast<double>(); P0=t[12].cast<double>();

	m0=t[13].cast<double>(); h0=t[14].cast<double>();
	set_hs(t[15].cast<double>(),t[16].cast<double>());
	c_star=t[17].cast<double>();

	mux=t[18].cast<std::unordered_map<double, double>>();
	dmux=t[19].cast<std::unordered_map<double, double>>();
	ddmux=t[20].cast<std::unordered_map<double, double>>();
	thetas=t[21].cast<py::array_t<double>>();
	Rs=t[22].cast<py::array_t<double>>();

	Wfs=t[23].cast<py::array_t<double>>();
	Wf_Ts=t[24].cast<py::array_t<double>>();
	Wf_ms=t[25].cast<py::array_t<double>>();

}

//----------------------------------------------------
//pybind things

void init_Scaling(py::module_ &m){
	py::class_<Scaling,Therm>(m, "Scaling").def(py::init<>(),"Scaling Object")
		.def(py::init<Therm *>(),"Scaling Object")
		.def(py::init<Therm *, double, double>(),"Scaling Object")
		.def_readwrite("alpha",&Scaling::alpha).def_readwrite("beta",&Scaling::beta)
		.def_readwrite("gamma",&Scaling::gamma).def_readwrite("delta",&Scaling::delta)
		.def_readwrite("Tc",&Scaling::Tc).def_readwrite("muc",&Scaling::muc)
		.def_readwrite("Pc",&Scaling::Pc).def_readwrite("nc",&Scaling::nc)
		.def_readwrite("P0",&Scaling::P0).def_readwrite("n0",&Scaling::n0)
		.def_readwrite("m0",&Scaling::m0).def_readwrite("h0",&Scaling::h0)
		.def_readwrite("theta0",&Scaling::theta0).def_readwrite("mu_code",&Scaling::mu_code)
		.def_readwrite("h3",&Scaling::h3).def_readwrite("h5",&Scaling::h5)
		.def_readwrite("c_star",&Scaling::c_star).def_readwrite("Rs",&Scaling::Rs)
		.def_readwrite("mux",&Scaling::mux).def_readwrite("dmux",&Scaling::dmux).def_readwrite("ddmux",&Scaling::ddmux)
		.def_readwrite("thetas",&Scaling::thetas)
		.def_readwrite("Wfs",&Scaling::Wfs)
		.def_readwrite("Wf_Ts",&Scaling::Wf_Ts).def_readwrite("Wf_ms",&Scaling::Wf_ms)
		.def_readwrite("pP",&Scaling::pP)//.def_readwrite("pWF",&Scaling::pWF)
		.def("get_thetas", &Scaling::get_thetas, "get thetas")
		.def("set_n0P0", &Scaling::set_n0P0, "set n0 and P0")
		.def("g_func", &Scaling::g_func, "g_func").def("h_func", &Scaling::h_func, "h_func")
		.def("get_mux", &Scaling::get_mux, "Get the value of mux at a specific")
		.def("get_muxs", &Scaling::get_muxs, "Get the values of mux at each Temp")
		.def("set_hs", &Scaling::set_hs, "Get the values of mux at each Temp")
		.def("__repr__", [](const Scaling &a) {
				return "<Scaling:T_c="+to_string(a.Tc)+"  mu_c="
				+to_string(a.muc) +">";
		}).def(py::pickle([](Scaling &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				Scaling p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		);
}


