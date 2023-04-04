
#include <cmath>
#include <typeinfo>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_multimin.h>  //for minimization
#include <omp.h>

#include "Embedded.hpp"
#include "misc.hpp"

using namespace std;

Embedded::Embedded(Therm * pP):
pP(pP){
    //Mean-Field
    // set_critical_exponents(4.0/3.0, 3.0/2.0);
    //3D-Ising
    set_critical_exponents(1.209, 1.564);
}

Embedded::Embedded(Therm * pP, double Tc, double muc):
pP(pP){
    this->Tc=Tc; this->muc=muc;
    //Mean-Field
    // set_critical_exponents(4.0/3.0, 3.0/2.0);
    //3D-Ising
    set_critical_exponents(1.209, 1.564);
}

void Embedded::set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs){
    Therm::set_Tmu( Temps, mu_Bs);
    clear_comp_flag();

    py::buffer_info buf = Temps.request();
    Rs = py::array_t<double>(buf.size); Rs.resize(buf.shape);
    // Rs_T = py::array_t<double>(buf.size); Rs_T.resize(buf.shape);
    // Rs_m = py::array_t<double>(buf.size); Rs_m.resize(buf.shape);

    pP->divs=2;//Constituent model needs two derivatives
    pP->N_threads=N_threads;
    pP->set_Tmu( Temps, mu_Bs);
}

void Embedded::clear_comp_flag(){
    this->comp_flag=false;
    pP->clear_comp_flag();
}


void Embedded::set_critical_exponents(double k, double p){
    this->k = k; this->p = p;
    this->alpha = 2.0 - k*p;
    this->beta  = p*(k - 1.0);
    this->gamma = p*(2.0 - k);
    this->delta = 1.0/(k - 1.0);
}

inline void Embedded::mux_insert(double Temp,
  double mu, double dmu, double ddmu){
    mux.insert(pair<double, double>(Temp,mu));
    dmux.insert(pair<double, double>(Temp,dmu));
    ddmux.insert(pair<double, double>(Temp,ddmu));
}

//use an elliptical mux
void Embedded::get_mux_ellip(double Temp, double mu0){
    if(mux.count(Temp)) return;

    double T0 = Tc / sqrt(1 - pow(muc/mu0,2));
    double zz = sqrt(1 - pow(Temp/T0,2));
    mux_insert(Temp, mu0*zz, -Temp*mu0/(T0*T0*zz) , mu0/(T0*T0*pow(zz,3) ));
}


//Find mu_x(T), the chemical potential where n_B=nc
void Embedded::get_mux(double Temp){
    if(mux.count(Temp)) return;

    int mu_code=3;

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

    if (Pc==-1) {
        double yS[6];
        pP->compute_single(Tc,muc,yS);
        nc = yS[2]*get_R(Tc); Pc = yS[0]*get_R(Tc);
    }


    //target density for n_BG
    double n_t=nc/get_R(Temp);
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
            if (PyErr_CheckSignals() != 0) throw py::error_already_set();
            mu1=mu2; n1=n2;
            mu2=1.5*mu2; n2=n_BG(Temp,mu2);

        }
    }
    mu=.5*(mu1+mu2);

    double nB=n_t*2,tol=1e-4;
    int iter=0;
    //compute mux
    while(abs(nB/n_t-1.0)>tol and abs(mu1-mu2)>tol and iter<80){
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
    //compute mup (mu_x(T+dT))
    iter=0; nB=n_t*2;
    n_t=nc/get_R(Temp+dT);
    mu2=mu; n2=n_BG(Temp+dT,mu2);
    mu1=.95*mu2; n1=n_BG(Temp+dT,mu1);

    while(n1>n_t){
        mu2=mu1; n2=n1;
        mu1=.9*mu2;
        n1=n_BG(Temp+dT,mu1);
    }
    while(n2<n_t){
        mu1=mu2; n1=n2;
        mu2=1.05*mu1;
        n2=n_BG(Temp+dT,mu2);
    }
    while(abs(nB/n_t-1.0)>tol  and abs(mu1-mu2)>tol and iter<80){
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
    iter=0; nB=n_t*2;
    n_t=nc/get_R(Temp-dT);
    mu1=mu; n1=n_BG(Temp-dT,mu1);
    mu2=1.1*mu1; n2=n_BG(Temp-dT,mu2);

    while(n1>n_t){
        mu2=mu1; n2=n1;
        mu1=.9*mu1;
        n1=n_BG(Temp-dT,mu1);
    }
    while(n2<n_t){
        mu1=mu2; n1=n2;
        mu2=1.05*mu1;
        n2=n_BG(Temp-dT,mu2);
    }
    while(abs(nB/n_t-1.0)>tol and abs(mu1-mu2)>tol and iter<80){
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

double Embedded::n_BG(double Temp, double mu_B){
    double yS[6];
    pP->compute_single(Temp,mu_B,yS);
    return yS[2];
}

void Embedded::get_Q(double x,double y, double * Qs){
    // Qs={Q,Q1,Q2,Q11,Q12,Q22};
    double Z = sqrt(x*x +y*y);

    Qs[0]=pow(Z + y,k);
    Qs[1]=0;Qs[2]=0;Qs[3]=0;Qs[4]=0;Qs[5]=0;
    if(y==0.0 and divs>0){
        Qs[1]=k*pow(x,k-1); Qs[2]=Qs[1];
        if(divs>1){
            Qs[3]=k*(k-1.0)*pow(x,k-2);
            Qs[4]=Qs[3]; Qs[5]=Qs[3]/(k-1.0);
        }
    }else if(divs>0){
        Qs[2]=k*Qs[0]/Z;
        Qs[1]=k*x*pow(Z + y,k-1)/Z;
        if(divs>1){
            if(y==-1.0) Qs[3]=(k/(Z*Z))*
                pow(Z+y,k-1)*(1.0/Z +(k-1)*(Z+1));
            else Qs[3]=(k/(Z*Z))*
                pow(Z+y,k-2)*(y*y*(1.0+y/Z) +(k-1.0)*x*x);

            Qs[4]=k*(Qs[1]/Z -x*Qs[0]/pow(Z,3));
            Qs[5]=(k*Qs[0]/(Z*Z))*(k-y/Z);
        }
    }
}

void Embedded::get_rx(double Temp, double mu, double * rxs){
    // rs={r,r_T,r_mu,r_TT,r_Tmu,r_mumu};

    if (PyErr_CheckSignals() != 0) throw py::error_already_set();

    rxs[1]=0;rxs[2]=0;rxs[3]=0;rxs[4]=0;rxs[5]=0;
    if(mu == 0){
        rxs[0] = -1.0;
        return;
    }
    int m = m_exponent;
    get_mux(Temp);

    double mx = mux[Temp], dmx = dmux[Temp], ddmx = ddmux[Temp];
    double M = pow(mu,m), Mx= pow(mx,m);
    double MM = M + Mx;

    rxs[0] = (M-Mx)/MM;
    if(divs>0){
        rxs[1] = -2*m*dmx*M*Mx/(MM*MM*mx);
        rxs[2] = 2*m*M*Mx/(MM*MM*mu);
        if(divs>1){
            rxs[5] =2*m*M*Mx*((m-1)*Mx - (m+1)*M)/(pow(MM,3)*mu*mu);
            rxs[4] = m*rxs[0]*rxs[2]*dmx/mx;
            rxs[3] = -2*m*M*Mx*(  (ddmx/mx)/pow(MM,2)+
                pow(dmx/mx,2)*((m-1)*M-(m+1)*Mx)/pow(MM,3));
        }
    }
}

void Embedded::get_a(double Temp, double * as){
    // as={a,a_T,a_TT};
    as[0]=a0*exp(-Temp/Ta);
    as[1]= -as[0]/Ta; as[2]= -as[1]/Ta;
}

void Embedded::get_Delta(double Temp, double * Delta){
    // Delta={Delta,Delta_T,Delta_TT};
    Delta[1] = 0.0;
    Delta[2] = 0.0;
    if(Temp <= Tc){
        Delta[0] = d_m*pow(1.0-Temp/Tc,p)*exp(-Temp/Td);
        if(divs>0) Delta[1] = -Delta[0]/Td
            -d_m*(p/Tc)*pow(1.0-Temp/Tc,p-1)*exp(-Temp/Td);
        if(divs>1) Delta[2] = -2*Delta[1]/Td - Delta[0]/(Td*Td)
            +d_m*(p*(p-1)/(Tc*Tc))*pow(1.0-Temp/Tc,p-2)*exp(-Temp/Td);
    } else {
        Delta[0] =d_p*pow(Temp/Tc-1.0,p)*exp(-Temp/Td);
        if(divs>0) Delta[1] = -Delta[0]/Td
            +d_p*(p/Tc)*pow(Temp/Tc-1.0,p-1)*exp(-Temp/Td);
        if(divs>1) Delta[2] = -2*Delta[1]/Td - Delta[0]/(Td*Td)
            +d_p*(p*(p-1)/(Tc*Tc))*pow(Temp/Tc-1.0,p-1)*exp(-Temp/Td);
    }
}

double Embedded::get_R(double Temp){
    // R={R,R_T,R_TT};
    double Q_p, Q_m, Delta;
    double as[3];

    get_a(Temp, as);

    if(Temp <= Tc)
        Delta = d_m*pow(1.0-Temp/Tc,p)*exp(-Temp/Td);
    else
        Delta = d_p*pow(Temp/Tc-1.0,p)*exp(-Temp/Td);

    Q_p=pow(sqrt(Delta*Delta +1.0) + 1.0,k);
    Q_m=pow(sqrt(Delta*Delta +1.0) - 1.0,k);

    if(Temp <= Tc)
        return 1 - as[0]*( Q_p - pow(Delta,k) );
    else
        return 1 - as[0]*( Q_p + Q_m - 2*pow(Delta,k) );
}

void Embedded::get_R(double Temp, double mu_B, double * R){
    // R={R,R1,R2,R11,R12,R22};

    double Q_pr[6], Q_mr[6], Q_p1[6], Q_m1[6];
    double Delta[3], rxs[6], as[3];

    get_Delta(Temp, Delta);

    get_a(Temp, as);

    get_rx( Temp, mu_B, rxs);

    get_Q(Delta[0],rxs[0], Q_pr);
    get_Q(Delta[0],-rxs[0], Q_mr);
    get_Q(Delta[0],1.0, Q_p1);
    get_Q(Delta[0],-1.0, Q_m1);

    R[1]=0;R[2]=0;R[3]=0;R[4]=0;R[5]=0;

    if(Temp <= Tc and rxs[0]>=0.0){

        R[0] = 1 + as[0]*( Q_pr[0] - Q_p1[0] );
        if(divs<1) return;
        R[1] = as[1]*( Q_pr[0] - Q_p1[0] ) + as[0]*(Delta[1]*( Q_pr[1] - Q_p1[1] )
             + rxs[1]* Q_pr[2]);
        R[2] =  as[0]*rxs[2]*Q_pr[2];
        if(divs<2) return;
        R[3] = as[2]*( Q_pr[0] - Q_p1[0] ) + 2*as[1]*(Delta[1]*( Q_pr[1] - Q_p1[1] )
            + rxs[1]* Q_pr[2])+ as[0]*( Delta[3]*(Q_pr[1] - Q_p1[1])
            + Delta[1]*Delta[1]*( Q_pr[3] - Q_p1[3] )+ rxs[3]* Q_pr[2]
            +rxs[1]*( 2*Delta[1]*Q_pr[4] + rxs[1]*Q_pr[5] ));

        R[4] =  as[1]*rxs[1]*Q_pr[2] + as[0]*Delta[1]*rxs[1]*Q_pr[4]
             + as[0]*(rxs[1]* Q_pr[5] + rxs[4]* Q_pr[2]);

        R[5] =  as[0]*(rxs[5]*Q_pr[2] + rxs[2]*rxs[2]*Q_pr[5]);

    }else if(Temp <= Tc and rxs[0]<0.0){


        R[0] = 1 + as[0]*( Q_mr[0] - Q_p1[0] );

        if(divs<1) return;

        R[1] = as[1]*( Q_mr[0] - Q_p1[0] ) + as[0]*(Delta[1]*( Q_mr[1] - Q_p1[1] )
             - rxs[1]* Q_mr[2]);
        R[2] =  -as[0]*rxs[2]*Q_mr[2];

        if(divs<2) return;

        R[3] = as[2]*( Q_mr[0] - Q_p1[0] ) + 2*as[1]*(Delta[1]*( Q_mr[1] - Q_p1[1] )
            - rxs[1]* Q_mr[2])+ as[0]*( Delta[3]*(Q_mr[1] - Q_p1[1])
            + Delta[1]*Delta[1]*( Q_mr[3] - Q_p1[3] )- rxs[3]* Q_mr[2]
            +rxs[1]*( -2*Delta[1]*Q_mr[4] + rxs[1]*Q_mr[5] ));

        R[4] =  -as[1]*rxs[1]*Q_mr[2] - as[0]*Delta[1]*rxs[1]*Q_mr[4]
             + as[0]*(rxs[1]*rxs[1]* Q_mr[5] - rxs[4]* Q_mr[2]);

        R[5] =  as[0]*(-rxs[5]*Q_mr[2] + rxs[2]*rxs[2]*Q_mr[5]);

    }else{
        R[0]= 1 + as[0]*( Q_pr[0] + Q_mr[0] - Q_p1[0] - Q_m1[0] );

        if(divs<1) return;

        R[1]= as[1]*( Q_pr[0] + Q_mr[0] - Q_p1[0] - Q_m1[0] )
            + as[0]*Delta[1]*( Q_pr[1] + Q_mr[1] - Q_p1[1] - Q_m1[1] )
            + as[0]*rxs[1]* ( Q_pr[2] - Q_mr[2] );
        R[2]=  as[0]*rxs[2]*( Q_pr[2] - Q_mr[2] );

        if(divs<2) return;

        R[3] = as[2]*( Q_pr[0] + Q_mr[0] - Q_p1[0]- Q_m1[0] )
            + 2*as[1]*(Delta[1]*( Q_pr[1] + Q_mr[1] - Q_m1[1] - Q_p1[1] )
            + rxs[1]* (Q_pr[2] - Q_mr[2]))+ as[0]*(Delta[3]*(Q_pr[1] + Q_mr[1] - Q_p1[1] - Q_m1[1])
            + Delta[1]*Delta[1]*( Q_pr[3] + Q_mr[3] - Q_m1[3] - Q_p1[3] )+ rxs[3]* (Q_pr[2]-  Q_mr[2])
            +rxs[1]*( 2*Delta[1]*(Q_pr[4]-Q_mr[4]) + rxs[1]*(Q_pr[5] +Q_mr[5])) );

        R[4] =  as[1]*rxs[1]*(Q_pr[2] - Q_mr[2]) + as[0]*Delta[1]*rxs[1]*(Q_pr[4]- Q_mr[4])
             + as[0]*(rxs[1]* (Q_pr[5]+ Q_mr[5]) + rxs[4]* (Q_pr[2] -  Q_mr[2]));

        R[5] =  as[0]*(rxs[5]*(Q_pr[2]-Q_mr[2]) + rxs[2]*rxs[2]*(Q_pr[5]+Q_mr[5]));
    }

}

void Embedded::get_muxs(){


    int N=Temps.request().size;
    double * pTemp = get_arr_pointer( Temps );

    double yS[6];
    pP->compute_single(Tc,muc,yS);
    nc = yS[2] * get_R(Tc); //Critical Baryon Density
    Pc = yS[0] * get_R(Tc); //Critical Pressure
    get_mux(Tc); //get mu_x for Tc

    #pragma omp parallel for schedule(dynamic) \
    num_threads(N_threads)
    for(int i=0;i<N;i++){
        if (PyErr_CheckSignals() != 0) throw py::error_already_set();
        bool prev_T=false;

        //Avoid calculating mu_x multiple times if Temps are the same
        for(int j=0; j<i;j++){
            prev_T=prev_T or (pTemp[j]==pTemp[i]);
        }
        if(not prev_T)
            // get_mux_ellip(pTemp[i], 970.0); //for  squiggly version
            get_mux(pTemp[i]);

    }
}


void Embedded::compute_single( double Temp, double mu_B, double * y){


    double z[6], R[6];
    pP->compute_single(Temp, mu_B, z);

    if(nc==0){
        double yS[6];
        pP->compute_single(Tc,muc,yS);
        nc = yS[2] * get_R(Tc); //Critical Baryon Density
    }

    get_R(Temp, mu_B, R);

    y[0] = z[0]*R[0];

    y[1]=0.0; y[2]=0.0; y[3]=0.0; y[4]=0.0; y[5]=0.0;

    if (divs<1) return;
    y[1] = z[1]*R[0] + z[0]*R[1];
    y[2] = z[2]*R[0] + z[0]*R[2];
    if (divs<2) return;
    y[3] = z[3]*R[0] + 2*z[1]*R[1] + z[0]*R[3];
    y[4] = z[4]*R[0] + z[1]*R[2] + z[2]*R[1] + z[0]*R[4];
    y[5] = z[5]*R[0] + 2*z[2]*R[2] + z[0]*R[5];

}


void Embedded::compute(){
    if(comp_flag) return;
    comp_flag=true;
    pP->compute();
    get_muxs();
    int N=Temps.request().size;

    double *pTemp= get_arr_pointer(Temps), *pmu= get_arr_pointer(mu_Bs);

    double *Z= get_arr_pointer( Y ), *Z_T= get_arr_pointer( Y_T ), *Z_m= get_arr_pointer( Y_m );
    double *Z_TT= get_arr_pointer( Y_TT ), *Z_Tm= get_arr_pointer( Y_Tm ),*Z_mm= get_arr_pointer( Y_mm );

    double *P= get_arr_pointer( pP->Y ), *P_T= get_arr_pointer( pP->Y_T ), *P_m= get_arr_pointer( pP->Y_m );
    double *P_TT= get_arr_pointer( pP->Y_TT ), *P_Tm= get_arr_pointer( pP->Y_Tm ),*P_mm= get_arr_pointer( pP->Y_mm );

    double *pRs= get_arr_pointer( Rs );

    #pragma omp parallel for schedule(dynamic) \
    shared(divs) num_threads(N_threads)
    for(int i=0;i<N;i++){

        //Display Progress
        if(omp_get_thread_num()%N_threads==0 and display_progress){
            cout<< i<<"/"<<N<<endl;
        }
        //Allow Ctrl-C abort
        if (PyErr_CheckSignals() != 0) throw py::error_already_set();

        double R[6];

        get_R(pTemp[i], pmu[i], R);

        pRs[i] = R[0];

        Z[i] = P[i]*R[0];
        Z_T[i]=0.0; Z_m[i]=0.0; Z_TT[i]=0.0; Z_Tm[i]=0.0; Z_mm[i]=0.0;

        if (divs<1) continue;
        Z_T[i] = P_T[i]*R[0] + P[i]*R[1];
        Z_m[i] = P_m[i]*R[0] + P[i]*R[2];
        if (divs<2) continue;
        Z_TT[i] = P_TT[i]*R[0] + 2*P_T[i]*R[1] + P[i]*R[3];
        Z_Tm[i] = P_Tm[i]*R[0] + P_T[i]*R[2]
            + P_m[i]*R[1] + P[i]*R[4];
        Z_mm[i] = P_mm[i]*R[0] + 2*P_m[i]*R[2] + P[i]*R[5];

        // if(mu_code==0) Z[i]=P[i]
        // else if(mu_code==1) Z[i]=P1;
        // else if(mu_code==2) Z[i]=P1;
    }
}

//Used by pythons pickling and copying features
py::tuple Embedded::getstate(){
    return py::make_tuple(12,Therm::getstate(),pP->getstate()
        ,k,p,Tc,muc,Pc,nc,a0,Ta,Td,d_p,d_m,m_exponent,mux,dmux,ddmux);
}

void Embedded::setstate(py::tuple t){
    if (t.size() != 17)
        throw std::runtime_error("Invalid state! (Embedded)");

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

    set_critical_exponents(t[3].cast<double>(), t[4].cast<double>());

    Tc=t[5].cast<double>(); muc=t[6].cast<double>();
    nc=t[7].cast<double>(); Pc=t[8].cast<double>();
    a0=t[9].cast<double>(); Ta=t[10].cast<double>();
    Td=t[11].cast<double>();
    d_p=t[12].cast<double>(); d_m=t[13].cast<double>();
    m_exponent=t[14].cast<int>();

    mux=t[15].cast<std::unordered_map<double, double>>();
    dmux=t[16].cast<std::unordered_map<double, double>>();
    ddmux=t[17].cast<std::unordered_map<double, double>>();

}





void init_Embedded(py::module_ &m){
    py::class_<Embedded,Therm>(m, "Embedded").def(py::init<>(),"Embedded Object")
        .def(py::init<Therm *>(),"Embedded Object")
        .def(py::init<Therm *, double, double>(),"Embedded Object")
        .def_readwrite("k",&Embedded::k).def_readwrite("p",&Embedded::p)
        .def_readwrite("alpha",&Embedded::alpha).def_readwrite("beta",&Embedded::beta)
        .def_readwrite("gamma",&Embedded::gamma).def_readwrite("delta",&Embedded::delta)
        .def_readwrite("Tc",&Embedded::Tc).def_readwrite("muc",&Embedded::muc)
        .def_readwrite("a0",&Embedded::a0).def_readwrite("Ta",&Embedded::Ta)
        .def_readwrite("Td",&Embedded::Td)
        .def_readwrite("m_exponent",&Embedded::m_exponent)
        .def_readwrite("d_p",&Embedded::d_p).def_readwrite("d_m",&Embedded::d_m)
        .def_readwrite("nc",&Embedded::nc).def_readwrite("Pc",&Embedded::Pc)
        .def_readwrite("mux",&Embedded::mux).def_readwrite("dmux",&Embedded::dmux).def_readwrite("ddmux",&Embedded::ddmux)
        .def_readwrite("pP",&Embedded::pP)
        .def_readwrite("Rs",&Embedded::Rs)
        .def("set_critical_exponents", &Embedded::set_critical_exponents,
         "Sets the values of critical exponents using k and p.")
        .def("get_mux", &Embedded::get_mux, "Get the value of mux at a specific")
        .def("get_muxs", &Embedded::get_muxs, "Get the values of mux at each Temp")
        .def("__repr__", [](const Embedded &a) {
                return "<Embedded:T_c="+to_string(a.Tc)+"  mu_c="
                +to_string(a.muc) +">";
        }).def(py::pickle([](Embedded &p){
                return p.getstate();
            }, [](py::tuple t){ // __setstate__
                Embedded p; /* Create a new C++ instance */
                p.setstate(t); /* Assign additional state */
                return p;
            })
        );
}


