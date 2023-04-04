#ifndef PARTICLE_H
#define PARTICLE_H

// Contains Particle struct which stores the properties of a particle species.

#include <string>
#include "misc.hpp"
#include <vector>


//Gives a particular flavor state in terms of baryon number,
//charge, and strangeness (charm)
struct Flavor {
	int B;
	int Q;
	int S;
	// int C = 0;
};

// Frequently used lists of Quantum numbers for particle species
//Baryons
const std::vector<Flavor> F_N = { //N Baryons
	{1,1,0},{1,0,0},
	{-1,-1,0},{-1,0,0}
};
const std::vector<Flavor> F_S = { //Sigma Barons
	{1,1,-1},{1,0,-1},{1,-1,-1},
	{-1,-1,1},{-1,0,1},{-1,1,1},
};
const std::vector<Flavor> F_X = { //Xi Baryons
	{1,0,-2},{1,-1,-2},
	{-1,0,2},{-1,1,2}
};
const std::vector<Flavor> F_Del = { //Delta Baryons
	{1,2,0},{1,1,0},{1,0,0},{1,-1,0},
	{-1,-2,0},{-1,-1,0},{-1,0,0},{-1,1,0}
};
const std::vector<Flavor> F_O = { //Omega Baryons
	{1,-1,-3},	{-1,1,3}
};
const std::vector<Flavor> F_L = { //Lambda Baryons
	{1,0,-1},	{-1,0,1}
};
//Mesons
const std::vector<Flavor> F_Iscl = { //Isoscalar Mesons (eta/eta', omega, phi, h, f)
	{0,0,0}
};
const std::vector<Flavor> F_Ivect = { //Unflavored Isovector Mesons (pi, rho, a, b)
	{0,1,0},{0,0,0},{0,-1,0}
};
const std::vector<Flavor> F_K = { //Kaons
	{0,1,1},{0,0,1},{0,0,-1},{0,-1,-1}
};

//Charmed Baryons
// const std::vector<Flavor> F_Lc = { //Charmed Lambda Baryons
// 	{1,1,0,1},	{-1,-1,0,-1}
// };
// const std::vector<Flavor> F_Sc = { //Charmed Sigma Barons
// 	{1,2,0,1},{1,1,0,1},{1,0,0,1},
// 	{-1,-2,0,-1},{-1,-1,0,-1},{-1,0,0,-1},
// };

// const std::vector<Flavor> F_Xc1 = { // Charged Charmed Xi Baryons (and Xi')
// 	{1,1,-1,1},     {-1,-1,1,-1}
// };
// const std::vector<Flavor> F_Xc0 = { //Neutral Charmed Xi Baryons (and Xi')
// 	{1,0,-1,1}, 	{-1,0,1,-1}
// };
// //const std::vector<Flavor> F_Xcc = { //Double Charmed Xi Baryons
// //	{1,2,0,2},{1,1,0,2},
// //	{-1,-2,0,-2},{-1,-1,0,-2}
// //};

// const std::vector<Flavor> F_Oc = { //Charmed Omega Baryons
// 	{1,0,-2,1},	{-1,0,2,-1}
// };


//Charmed Mesons
// const std::vector<Flavor> F_D1 = { //Charged D Mesons
// 	{0,1,0,1},{0,-1,0,-1}
// };
// const std::vector<Flavor> F_D0 = { //Neutral D Mesons
// 	{0,0,0,1},{0,0,0,-1}
// };
// const std::vector<Flavor> F_Ds = { //Strange D Mesons
// 	{0,1,1,1},{0,-1,-1,-1}
// };

//-----------------------------------------------------------------------------

struct Particle{
	//stores properties of one species of particle

	Particle(std::string name, double mass, int J_degen, std::vector<Flavor> Fs,
	 STATISTICS stats);

	//Gives number of states (simply size of Fs)
	int degen();

	//Gives number of states with baryon number B
	int degen(int B);

	// //Gives number of states with a given baryon
	// //number B and charge Q.
	// int degen(int B, int Q);

	std::string name;
	double mass;
	int J_degen; //Spin Degeneracy (2J+1)

	//List of
	std::vector<Flavor> Fs;

	STATISTICS stats;
};

#endif