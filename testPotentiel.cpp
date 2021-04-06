#include <iostream>
#include "Vecteur2D.h"
#include "Montagne.h"
#include "Champ_potentiel.h"
#include <vector>
using namespace std;


int main() {
    Montagne Himalaya(3, 3, 0.4, 3, 3);
    ChampPotentiels CP(30,30,30,0.69);
    CP.set_xyz();
    CP.initialise(15,Himalaya);
    CP.calcule_laplaciens();
    CP.affichage();

    return 0;
}
