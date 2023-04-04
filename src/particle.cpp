#include <string>
#include "misc.hpp"
#include "particle.hpp"

#include <vector>

// Contains Particle struct which stores the properties of a particle species.

using namespace std;

Particle::Particle(string name, double mass, int J_degen,
 vector<Flavor> Fs,	STATISTICS stats){
	this->name = name;
	this->mass = mass;
	this->J_degen = J_degen;
	this->Fs = Fs;
	this->stats = stats;
}


//Gives number of states (simply size of Fs times spin degeneracy)
int Particle::degen() {
	return J_degen*Fs.size();
}

//degeneracy of states with a given baryon number B
int Particle::degen(int B){
	int X=0;
	for(auto F : Fs){
		if(F.B==B) X++;
	}
	return J_degen*X;
}

// //degeneracy of states with a given baryon number B
// //and charge Q.
// int Particle::degen(int B, int Q){
// 	int X=0;
// 	for(auto F : Fs){
// 		if(F.B==B and F.Q=Q) X++;
// 	}
// 	return J_degen*X;
// }

