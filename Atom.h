#ifndef ATOM_H
#define ATOM_H
#include <armadillo>

using namespace arma;

class Atom
{
private:
    vec position;
    vec velocity;
    vec force;
//    vec potential;
    vec position_0;

public:
    Atom();
    void setPosition(const vec& position);
    vec getPosition();
    void setVelocity(vec velocity);
    vec getVelocity();
    void setForce(vec force);
    vec getForce();
    void setPosInit(vec position_0);
    vec getPosInit();
};

#endif // ATOM_H
