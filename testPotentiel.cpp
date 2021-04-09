#include <iostream>
#include "Vecteur2D.h"
#include "Montagne.h"
#include "Champ_potentiel.h"
#include <vector>
using namespace std;


int main() {
    Montagne Himalaya(15, 15, 15, 5, 5);
    ChampPotentiels CP(30,30,30,0.689655);
    CP.initialise(15,Himalaya);
    //CP.calcule_laplaciens();
    //CP.affichage();
    //CP.lapla_affichage();
    Montagne m(15, 15, 15, 5, 5);
    ChampPotentiels champ(30,30,30,0.689655);
    champ.initialise(15,m);
    //champ.calcule_laplaciens();
    champ.resolution(2.2621843e-5,2000,false);
    champ.affiche_total();
    return 0;
}
