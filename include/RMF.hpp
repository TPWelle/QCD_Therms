#ifndef RMF_H
#define RMF_H

#include "misc.hpp"
#include "Therm.hpp"
#include <vector>
#include <gsl/gsl_multiroots.h>


struct RMF : Therm {

	RMF()=default;

	void set_Tmu(py::array_t<double> Temps, py::array_t<double> mu_Bs) override;

	int N_Mes=4; int N_Bar=8;
	// int N_Mes=4; int N_Bar=2;

	double b = 0.0; double c = 0.0;

	std::vector<double> m_baryons={938.27,939.57,1115.7,1189.4,1192.6,1197.4,1314.9,1321.7};
	std::vector<double> m_mesons={500.0, 783.0, 1019.5, 775.5};


	double g_sig;
	double g_om;
	double g_rho;
	double alpha_V = 1.0; 
	double z = 1.0;//z =   g8/(g1 sqrt(6))


	//p,n,Lambda,Sigma,Sigma,Sigma,Xi,Xi
	py::array_t<double> gs;


	void set_gs();


	//get chemical potential for each baryon
	std::vector<double> get_mus(double mu_B); //mu_Q, mu_S
	py::array_t<double> Meson_RMFs;

	//vector and scalar density of each baryon species

	py::array_t<double> ns_Bar;
	py::array_t<double> ns_scl_Bar;

	//For Root Finding prev
	std::vector<double> get_m_star(double sig);
	std::vector<double> get_mu_star( std::vector<double> mus, std::vector<double> RMFs);
	std::vector<double> get_n(double Temp, std::vector<double> m_star, std::vector<double> mu_star);
	std::vector<double> get_n_scl(double Temp, std::vector<double> m_star, std::vector<double> mu_star);

	std::vector<double> RMF_obj(std::vector<double> RMFs,
		std::vector<double> ns, std::vector<double> n_scl);
	std::vector<double> RMF_obj_test(std::vector<double> RMFs, double Temp, double mu_B);
	gsl_vector * RMFs_initialize( double Temp, double mu_B);

	std::vector<double> get_RMFs(double Temp, double mu_B);
	double get_sig(double Temp);

	//Used for python's serialization and copying functions
    py::tuple getstate() override;
    void setstate(py::tuple t) override;

	void compute_single(double Temp, double mu_B, double * y) override;
	void compute() override;

};

void init_RMF(py::module_ &m);

#endif