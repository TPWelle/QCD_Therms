#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <omp.h>

#include "RMF.hpp"
#include <gsl/gsl_errno.h>          //for modifying gsl's error handling
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_result.h>      //contains a struct for returning results 
                                    //of gsl functions
#include "misc.hpp"
#include "PT.hpp"

using namespace std;


vector<double> RMF::get_m_star(double sig){
	double *pgs= get_arr_pointer(gs); //!
	//Compute m_stars for baryons given mean scalar field sig
	vector<double> m_stars={};
	for(int k=0;k<N_Bar;k++){
		double m_star=m_baryons[k]-pgs[0*N_Bar+k] * sig;
		m_stars.push_back(m_star);
	}
	return m_stars;
}

vector<double> RMF::get_mu_star( vector<double> mus, vector<double> RMFs){
	//Compute mu_star for baryons given mean field values
	double *pgs= get_arr_pointer(gs);
	vector<double> mu_stars={};
	for(int k=0;k<N_Bar;k++){
		double mu_star=mus[k];
		for(int j=1; j<N_Mes; j++){
			mu_star=mu_star - pgs[j*N_Bar+k]*RMFs[j];
		}
		mu_stars.push_back(mu_star);
	}
	return mu_stars;
}


//calculate the chemical potential for each baryon given mu_B (mu_Q)
//Currently just returns mu_B for each
vector<double> RMF::get_mus(double mu_B){
	vector<double> mus (N_Bar,mu_B);
	return mus;
}
vector<double> RMF::get_n_scl(double Temp, vector<double> m_star, vector<double> mu_star){
	vector<double> n_scl;
	for(int j=0;j<N_Bar;j++){

		//This may cause issues when |m_star| is very small

		double y = Temp/abs(m_star[j]), z = mu_star[j]/abs(m_star[j]);

		n_scl.push_back((2*pow(m_star[j],3)/(pi*pi))
			*(Thermal_Integral(y, z, 0, 1, FERMION) + Thermal_Integral(y, -z, 0, 1, FERMION)));

	}
	return n_scl;
}
vector<double> RMF::get_n(double Temp, vector<double> m_star, vector<double> mu_star){
	vector<double> ns;
	for(int j=0;j<N_Bar;j++){

		//This may cause issues when |m_star| is very small

		double y = Temp/abs(m_star[j]), z = mu_star[j]/abs(m_star[j]);
		ns.push_back((2*pow(abs(m_star[j]),3)/(pi*pi))
			*(Thermal_Integral(y, z, 1, 1, FERMION) - Thermal_Integral(y, -z, 1, 1, FERMION)) );

	}
	return ns;
}


vector<double> RMF::RMF_obj(vector<double> RMFs,
	vector<double> ns, vector<double> n_scl){

	vector<double> result(N_Mes,0.0);
	double *pgs= get_arr_pointer(gs);
	double gN = pgs[0];

	result[0] =(b*939*gN*pow(RMFs[0]*gN,2)+c*gN*pow(RMFs[0]*gN,3));
	result[0] += m_mesons[0]*m_mesons[0] * RMFs[0];
	for(int j=0;j<N_Bar;j++){
		result[0]-=pgs[0*N_Bar+j]*   n_scl[j];
	}

	for(int k=1;k<N_Mes;k++){
		result[k] += m_mesons[k]*m_mesons[k] * RMFs[k];
		for(int j=0;j<N_Bar;j++){
			result[k]-=pgs[k*N_Bar+j]*   ns[j];
		}
	}
	return result;
}


vector<double> RMF::RMF_obj_test(vector<double> RMFs, double Temp, double mu_B){

	vector<double> m_star = get_m_star(RMFs[0]);
	vector<double> mu_star = get_mu_star(get_mus(mu_B), RMFs);

	vector<double> n_scl = get_n_scl(Temp, m_star, mu_star);
	vector<double> ns = get_n(Temp, m_star, mu_star);

	return RMF_obj(RMFs, ns, n_scl);
}



//Parameters for the root finder
struct rf_params {double Temp; double mu_B; RMF * pRMF; };

//Objective function for root finding
int rf_obj (const gsl_vector * x, void * p, gsl_vector * f) {
    struct rf_params * params= (rf_params *) p;
    const double Temp = (params->Temp);
    const double mu_B = (params->mu_B);
    RMF *pRMF = (params->pRMF);

    const int N = (*x).size;
    vector<double> RMFs (N, 0.0);

    for(int i=0;i<N;i++){
        RMFs[i]=gsl_vector_get(x,i);
    }
    vector<double> m_stars=pRMF->get_m_star(RMFs[0]);
    vector<double> mu_stars=pRMF->get_mu_star(pRMF->get_mus(mu_B), RMFs);
	vector<double> ns = pRMF->get_n(Temp, m_stars, mu_stars);
	vector<double> n_scl = pRMF->get_n_scl(Temp, m_stars, mu_stars);
	vector<double> obj=pRMF->RMF_obj(RMFs, ns, n_scl);
    for(int i=0;i<N;i++){
    	// cout<<RMFs[i]<<":"<<obj[i]/1.0e6<<"|";
        gsl_vector_set(f, i, obj[i]);
    }
    // cout<<endl;

	if(RMFs[0]<0.0) obj[0] = 1e10 * RMFs[0];

    return GSL_SUCCESS;
}




double sig_obj (double  x, void * p) {
    struct rf_params * params= (rf_params *) p;
    const double Temp = (params->Temp);
    const double mu_B = (params->mu_B);
    RMF *pRMF = (params->pRMF);
    double N_Bar = pRMF->N_Bar;
    vector<double> m_stars=pRMF->get_m_star(x);
    vector<double> mu_stars (N_Bar, mu_B) ;

	vector<double> n_scl = pRMF->get_n_scl(Temp, m_stars, mu_stars);


	double result = 0;
	double *pgs= get_arr_pointer(pRMF->gs);
	double gN = pgs[0];

	result =(pRMF->b*939*gN*pow(x*gN,2)+pRMF->c*gN*pow(x*gN,3));
	result += pow(pRMF->m_mesons[0],2)* x;
	for(int j=0;j<N_Bar;j++){
		result -= pgs[j]*(n_scl[j]);
		// result -= 8.0*(n_scl[j]);
	}

	if(x<0.0) result = 1e10 * x;
    return result;
}

//Calculates the Mean field value of sigma when mu=0
double RMF::get_sig(double Temp){

	gsl_root_fsolver *s;

	int status;
	size_t NN, iter = 0;
	struct rf_params params = { Temp, 0.0, this};
	double x_lo, x_hi;

	if(Temp<100){
		x_lo = 0.0; x_hi = 5.0;
	} else if(Temp<200){
		x_lo = 0.0; x_hi = 100.0;
	} else {
		x_lo = 50.0; x_hi = 200.0;
	}
	double r = .5*(x_lo + x_hi);

	gsl_function F;
	F.function = &sig_obj; F.params = &params;
	const gsl_root_fsolver_type *T;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);


	double epsabs=1.0e-18,epsrel=1.0e-7;

	do {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi,  epsabs,epsrel);

		// printf ("%5d [%.7f, %.7f] %.7f %.7f %.7f\n",
		//       iter, x_lo, x_hi, r, x_hi - x_lo, sig_obj(r, &params));
	} while (status == GSL_CONTINUE && iter < 400);

	return r;
}


//Calculates the Mean field values of the various meson species
//at the given temperature and baryon chemical potential
vector<double> RMF::get_RMFs(double Temp, double mu_B){
	vector<double> RMFs;
	if (mu_B == 0){
		// cout<<"a"<<endl;
		RMFs.push_back(get_sig(Temp));
		return RMFs;
	}
	gsl_multiroot_fsolver *s;

	int status;
	size_t iter = 0;
	struct rf_params params = { Temp, mu_B, this};

	gsl_multiroot_function F = {&rf_obj,N_Mes,&params};

    vector<double> mus=get_mus(mu_B);
    vector<double> ns;
    if(mu_B==0){
    	for(int i=0; i<N_Bar; i++) ns.push_back(0.0);
    } else ns=get_n(Temp, m_baryons, mus);

	//allocate and initialize vector
	gsl_vector *x =RMFs_initialize( Temp, mu_B);

	// cout<<" T:"<<Temp<<"  mu:"<<mu_B<<endl;
	// cout<<"==============================="<<endl;
	const gsl_multiroot_fsolver_type *T;
	// T = gsl_multiroot_fsolver_hybrids;
	T = gsl_multiroot_fsolver_broyden;
	// T = gsl_multiroot_fsolver_dnewton;
	s = gsl_multiroot_fsolver_alloc(T, N_Mes);
	gsl_multiroot_fsolver_set (s, &F, x);

	double res=1.0e-7;
	double epsabs=1.0e-8, epsrel=1.0e-6;
	do{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);

		if (status) break; /* check if solver is stuck */

		// check for convergence using function residual
		status =gsl_multiroot_test_residual(s->f, res);

		// check for convergence using change in root
		// status =gsl_multiroot_test_delta(gsl_multiroot_fsolver_dx(s),
		// 	gsl_multiroot_fsolver_root(s), epsabs, epsrel);
	}while (status == GSL_CONTINUE && iter < 200);

	for(int k = 0; k<N_Mes; k++){
		RMFs.push_back(gsl_vector_get(s->x,k));
	}
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);

	return RMFs;
}



// gsl_vector * RMF::sig_n_initialize( double Temp, double mu_B){
gsl_vector * RMF::RMFs_initialize( double Temp, double mu_B){

	int NN; if(mu_B==0) NN=1;
	else NN=N_Mes;
	gsl_vector *x=gsl_vector_alloc(NN);
	//Initial Guess for mu_Bs
	double sig_init=1.0;
	if(Temp>150) sig_init+=5;
	if(Temp>170) sig_init+=20;
	if(Temp>200) sig_init+=55;
	if(Temp>300) sig_init+=40;
	if(Temp>400) sig_init+=55;

	if(mu_B>700) sig_init+=20+3*g_sig; 
	if(mu_B>1500) sig_init+=20+3*g_sig;

	gsl_vector_set(x,0, sig_init);
	double n_init=mu_B/200.0;
	if(mu_B >700) n_init+=(mu_B -600)/10;
	if(Temp > 180) n_init=mu_B/10;

	if(mu_B != 0.0){
		gsl_vector_set(x,1,n_init);
		gsl_vector_set(x,2,-.4*n_init);
		gsl_vector_set(x,3,.1);
	}
	return x;
}


void RMF::compute_single(double Temp, double mu_B, double * y) {

	if(Temp <=0.0)
		throw domain_error("RMF::compute_single expects positive temperature.");

	vector<double> RMFs= get_RMFs(Temp, mu_B);
	vector<double> mus=get_mus(mu_B);
	vector<double> m_star = get_m_star(RMFs[0]);
	vector<double> mu_star = get_mu_star(mus, RMFs);
	vector<double> ns = get_n(Temp, m_star, mu_star);

	double P=0.0, eps=0.0, n_B=0.0 , s= 0.0;


	for(int j=0;j<N_Bar;j++){
		double y = Temp/abs(m_star[j]), z = mu_star[j]/abs(m_star[j]);
		P += (pow(m_star[j],4)/(3*pi*pi))
			*(Thermal_Integral(y, z, 0, 3, FERMION) + Thermal_Integral(y, -z, 0, 3, FERMION));
		eps += (pow(m_star[j],4)/(pi*pi))
			*(Thermal_Integral(y, z, 2, 1, FERMION) + Thermal_Integral(y, -z, 2, 1, FERMION));
		n_B += ns[j];
		s -= ns[j]*mus[j];
	}

	double A=pow(m_mesons[0]*RMFs[0],2)/2.0 + b*pow(RMFs[0],3)/3.0 + c*pow(RMFs[0],4)/4.0;
	double B=0;
	for(int i=1; i<N_Mes; i++) B+=pow(m_mesons[i]*RMFs[i],2)/2;

	P+= (B-A);
	eps+=(B+A);
	s+= eps+P;
	s=s/Temp;

	y[0]=P; y[1]=s; y[2]=n_B;
	// if(divs>0) 	if(divs>1)
			//IOU 2nd derivatives
}

void RMF::set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) {

	Therm::set_Tmu( Temps, mu_Bs);
	py::buffer_info buf=Temps.request();
	vector<py::ssize_t> M_shape=buf.shape;
	vector<py::ssize_t> B_shape=buf.shape;

	//Allocate Array for mean field values
	M_shape.push_back(N_Mes);
	Meson_RMFs = py::array_t<double>(buf.size*N_Mes); Meson_RMFs.resize(M_shape);

	//Allocate array for Baryon densities
	B_shape.push_back(N_Bar);
	ns_Bar = py::array_t<double>(buf.size*N_Bar); ns_Bar.resize(B_shape);
	ns_scl_Bar = py::array_t<double>(buf.size*N_Bar); ns_scl_Bar.resize(B_shape);

	// Use with Particle List

	// //Allocate Array for mean field values
	// M_shape.push_back(Mesons.size());
	// Meson_RMFs = py::array_t<double>(buf.size*Mesons.size()); Meson_RMFs.resize(M_shape);
	// //Allocate array for Baryon densities
	// B_shape.push_back(Baryons.size());
	// ns_Bar = py::array_t<double>(buf.size*Baryons.size()); ns_Bar.resize(B_shape);
	// ns_scl_Bar = py::array_t<double>(buf.size*Baryons.size()); ns_scl_Bar.resize(B_shape);
}

//Setting Coupling Constants

void RMF::set_gs() {
	//z =   g8/(g1 sqrt(6))
	vector<py::ssize_t> gs_shape={N_Mes,N_Bar};
	// vector<py::ssize_t> gs_shape={Mesons.size(),Baryons.size()};
	// p,n,Lambda, Sigma,Sigma,Sigma, Xi,Xi
	double denom = ( (4*alpha_V-1)*z +6 );


	//sigma-hyperon coupling ratios
	double R_L= 0.5033, R_S=0.5033, R_X = 0.2365;


	gs = py::array_t<double>(N_Bar*N_Mes); gs.resize(gs_shape);
	double *pgs= get_arr_pointer(gs);
	// cout<<N_Bar<<"  "<<N_Mes<<endl;
	// cout<<"g: "<<g_sig<<"  "<<g_om<<"  "<<g_rho<<endl;
	//=====sigma======
	pgs[0*N_Bar+0]=g_sig; pgs[0*N_Bar+1]=g_sig;
	pgs[0*N_Bar+2]=R_L*g_sig;
	pgs[0*N_Bar+3]=R_S*g_sig; pgs[0*N_Bar+4]=R_S*g_sig; pgs[0*N_Bar+5]=R_S*g_sig;
	pgs[0*N_Bar+6]=R_X*g_sig; pgs[0*N_Bar+7]=R_X*g_sig;

	//=====omega======
	double X;
	pgs[1*N_Bar+0]=g_om; pgs[1*N_Bar+1]=g_om; //N

	pgs[1*N_Bar+2]=g_om*(-2*(1-alpha_V)*z +6)/denom; //Lambda
	X=g_om*(2*(1-alpha_V)*z +6)/denom; //Sigma
	for(int i=3;i<6;i++) pgs[1*N_Bar+i]=X;
	X=g_om*(-(2*alpha_V+1)*z +6)/denom; //Xi
	pgs[1*N_Bar+6]=X; pgs[1*N_Bar+7]=X;

	//======phi=======

	X=g_om*((4*alpha_V-1)*z -3)*sqrt(2)/denom; //N

	pgs[2*N_Bar+0]=X; pgs[2*N_Bar+1]=X;
	pgs[2*N_Bar+2]=g_om*(-2*(1-alpha_V)*z -3)*sqrt(2)/denom; //Lambda
	X=g_om*(2*(1-alpha_V)*z -3)*sqrt(2)/denom; //Sigma
	for(int i=3;i<6;i++) pgs[2*N_Bar+i]=X;
	X=g_om*(-(2*alpha_V +1)*z -3)*sqrt(2)/denom; //Xi
	pgs[2*N_Bar+6]=X; pgs[2*N_Bar+7]=X;

	//======rho=======

	pgs[3*N_Bar+0]=(.5*g_rho); pgs[3*N_Bar+1]=(-.5*g_rho); 
	pgs[3*N_Bar+2]=(0.0);
	// pgs[3*N_Bar+3]=(.5*g_rho); pgs[3*N_Bar+4]=(0.0); pgs[3*N_Bar+5]=(-.5*g_rho);
	pgs[3*N_Bar+3]=g_rho; pgs[3*N_Bar+4]=0.0; pgs[3*N_Bar+5]=-g_rho;
	pgs[3*N_Bar+6]=(.5*g_rho); pgs[3*N_Bar+7]=(-.5*g_rho);
}


void RMF::compute() {

	if(comp_flag) return;

	set_gs();

	int N=Temps.request().size;
	double *pTemp= get_arr_pointer(Temps), *pmu= get_arr_pointer(mu_Bs);
	double *pY= get_arr_pointer(Y), *pY_T= get_arr_pointer(Y_T), *pY_m= get_arr_pointer(Y_m);
	double *pY_TT= get_arr_pointer(Y_TT), *pY_Tm= get_arr_pointer(Y_Tm), *pY_mm= get_arr_pointer(Y_mm);

	double *pns_Bar= get_arr_pointer( ns_Bar );
	double *pns_scl_Bar= get_arr_pointer( ns_scl_Bar );
	double *pMeson_RMFs= get_arr_pointer( Meson_RMFs );

	#pragma omp parallel for schedule(dynamic) \
	shared(divs) num_threads(N_threads)
	for(int i=0;i<N;i++){
		if (PyErr_CheckSignals() != 0) throw py::error_already_set();
		if(omp_get_thread_num()%N_threads==0 and display_progress)
			cout<< i<<"/"<<N<<endl;

		vector<double> RMFs= get_RMFs(pTemp[i], pmu[i]);
		vector<double> mus=get_mus(pmu[i]);

		vector<double> m_star = get_m_star(RMFs[0]);
		vector<double> mu_star = get_mu_star(mus, RMFs);
		vector<double> ns = get_n(pTemp[i], m_star, mu_star);
		vector<double> n_scl=get_n_scl(pTemp[i], m_star,  mu_star);



		// for computing X_mm
		double dmu =.2;

		vector<double> m_light_mesons
		 ={139.57,134.98,547.86,475.0,775.3,782.7,957.8,990.0,980.0, 1019.5};
		vector<double> degen_light_mesons
		 ={2,1,1,1,3,3,1,1,1,3};

		double P=0.0, eps=0.0, n_B=0.0 , s= 0.0;
		double X_2=0.0;

		for(int j=0;j<m_light_mesons.size();j++){
			double mm=m_light_mesons[j];
			double deg=degen_light_mesons[j];

			double y = pTemp[i]/abs(mm);
			P+=deg*(pow(mm,4)/(6*pi*pi))*Thermal_Integral(y,0,0,3, BOSON);
			eps+=deg*(pow(mm,4)/(2*pi*pi))*Thermal_Integral(y,0,2,1, BOSON);
		}

		for(int j=0;j<N_Mes;j++){
			pMeson_RMFs[i*N_Mes + j] = RMFs[j];
		}

		for(int j=0;j<N_Bar;j++){
			double y = pTemp[i]/abs(m_star[j]), z = mu_star[j]/abs(m_star[j]);

			P+=(pow(m_star[j],4)/(3*pi*pi))*
			(Thermal_Integral(y, z, 0, 3, FERMION) + Thermal_Integral(y, -z, 0, 3, FERMION));
			eps += (pow(m_star[j],4)/(pi*pi))*
			(Thermal_Integral(y, z, 2, 1, FERMION) + Thermal_Integral(y, -z, 2, 1, FERMION));
			n_B+=ns[j];
			// X_2+=(ns_p[j] - ns_n[j])/(2*dmu);
			pns_Bar[i*N_Bar + j] = ns[j];
			pns_scl_Bar[i*N_Bar + j] = n_scl[j];
			s-= ns[j]*mus[j];
		}

		double A=pow(m_mesons[0]*RMFs[0],2)/2.0 + b*pow(RMFs[0],3)/3.0 + c*pow(RMFs[0],4)/4.0;

		double B=0;
		for(int i=1; i<N_Mes; i++) B+=pow(m_mesons[i]*RMFs[i],2)/2;

		P+= (B-A);
		eps+=(B+A);
		s+= eps+P;
		s=s/pTemp[i];

		pY[i]=P; pY_T[i]=s; pY_m[i]=n_B;
		//2nd derivs are 0 for now
		pY_TT[i]=0.0; pY_Tm[i]=0.0; pY_mm[i]=X_2;
	}
	comp_flag=true;
}

//Used by pythons pickling and copying features
py::tuple RMF::getstate(){
	return py::make_tuple(5,Therm::getstate(),N_Mes,N_Bar,b,c,
		m_baryons,m_mesons,gs,Meson_RMFs);
}

void RMF::setstate(py::tuple t){
	if (t.size() != 9)
		throw std::runtime_error("Invalid state!(RMF)");

	Therm::setstate(t[1]);

	N_Mes=t[2].cast<int>(); N_Bar=t[3].cast<int>();
	b=t[4].cast<double>(); c=t[5].cast<double>();

	m_baryons=t[5].cast<vector<double>>();
	m_mesons=t[6].cast<vector<double>>();
	gs = t[7].cast<py::array_t<double>>();
	Meson_RMFs = t[8].cast<py::array_t<double>>();
}

void init_RMF(py::module_ &m){
	py::class_<RMF,Therm>(m, "RMF").def(py::init<>(),"RMF model Object")
		.def(py::init<>(),"RMF model Object")
		.def_readwrite("N_Mes",&RMF::N_Mes).def_readwrite("N_Bar",&RMF::N_Bar)
		.def_readwrite("b",&RMF::b).def_readwrite("c",&RMF::c)
		.def_readwrite("m_baryons",&RMF::m_baryons).def_readwrite("m_mesons",&RMF::m_mesons)
		.def_readwrite("g_sig",&RMF::g_sig).def_readwrite("g_om",&RMF::g_om)
		.def_readwrite("g_rho",&RMF::g_rho).def_readwrite("gs",&RMF::gs)
		.def_readwrite("Meson_RMFs",&RMF::Meson_RMFs)
		.def_readwrite("ns_Bar",&RMF::ns_Bar).def_readwrite("ns_scl_Bar",&RMF::ns_scl_Bar)
		.def("get_n_scl",&RMF::get_n_scl).def("get_n",&RMF::get_n)
		.def("get_m_star",&RMF::get_m_star).def("get_mu_star",&RMF::get_mu_star)
		.def("get_mus",&RMF::get_mus).def("RMF_obj",&RMF::RMF_obj)
		.def("RMF_obj_test",&RMF::RMF_obj_test)
		.def("set_gs",&RMF::set_gs)
		.def("__repr__", [](const RMF &a) {
				return "<RMF>";
		}).def(py::pickle( [](RMF &p){
				return p.getstate();
			}, [](py::tuple t){ // __setstate__
				RMF p; /* Create a new C++ instance */
				p.setstate(t); /* Assign additional state */
				return p;
			})
		);
}
