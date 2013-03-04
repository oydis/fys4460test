#include <iostream>
#include <armadillo>
#include "Lattice.h"
#include "Atom.h"
#include "Integrator.h"

using namespace std;
using namespace arma;
int main()
{
    int Nc = 8;
    int dim = 3;
    double b = 1.545; // = 5.260;
    int numpart = 4*Nc*Nc*Nc;

//    double kB = 1.3806488e-23; // m²kg/(s²K)
//    double T = 100; // K
//    double m = 39.948*1.660538921e-27; // kg
//    double sigma = sqrt(kB*T/m); // m/s


    double T_bath = 300/119.74; // K
    double eps = 1.0318e-2; // eV
    double m = 39.948; // u
    double v0 = sqrt(T_bath);

    Lattice Lat;
    Lat.createLattice(Nc,b,dim,v0);
    Atom* atoms = Lat.getAtoms();

    double dt = 0.006;
    double sig = 3.405; //Å
    //double eps = 119.8e20*kB; // Å²kg/s²



    Integrator Int(dt,m,numpart,dim, atoms,b,Nc,sig,eps,T_bath);
    Int.integrate();

}

