#pragma once
#include "Vecteur2D.h"
#include "Montagne.h"
#include <vector>


class CubedAir{

};

class Ciel {
private:
    std::vector<std::vector<std::vector<CubedAir>>> collection3D;
    unsigned int Nx,Ny,Nz;
    double lambda;
};


