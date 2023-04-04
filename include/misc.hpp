#ifndef MISC_H
#define MISC_H

const int MAX_NUMBER_ITER = 200;  //max number of iterations for solvers I coded

//Useful numbers

const double pi = 3.14159265358979323846;
const double gammaE=0.57721566490153286060; //Euler's gamma

//Zeta function at 3 and 5
const double Zeta3=1.2020569031595942854;
const double Zeta5=1.0369277551433699263;

//log(2) computed for convenience and speed
const double L2=0.69314718055994530942;

//logarithmic derivative of Zeta function
//at -1 and -3
const double Z1=1.98505372440541115057;
const double Z3=0.64542916293291613733;


//specifies the quantum statistics
enum STATISTICS : int {FERMION = 1, BOSON = -1, CLASSICAL = 0};

#endif