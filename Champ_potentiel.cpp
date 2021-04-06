#include "Champ_potentiel.h"
#include "Vecteur2D.h"
#include "Montagne.h"
#include <vector>
using namespace std;

Vecteur2D Potentiel::get_potentiel() const {
    return potentiel;
}

Vecteur2D Potentiel::get_laplacien() const {
    return laplacien;
}

void Potentiel::set_potentiel(double x, double y) {
    potentiel.set_coord(x,y);
}

void Potentiel::set_laplacien(double x, double y) {
    laplacien.set_coord(x,y);
}


// méthode initialisant le tableau contenant les coordonnées xi yj et zk -- fonctionne
void ChampPotentiels:: set_xyz(){
    double r;
    for(double i(0); i < Nx; ++i){
        X.push_back(i*lambda);
    }
    for(double j(0); j < Ny; ++j){
        Y.push_back(j*lambda - (lambda*(Ny-1))/2);
    }
    for(double k(0); k < Nz; ++k){
        Z.push_back(k*lambda);
    }


}

// méthode initialisant le potentiel vecteur en tout point (xi,yj,zk) -- fonctionne
void ChampPotentiels:: initialise(double const& v_infini, Montagne const& M) {
    //initialisation du potentiel vecteur en tout point et du laplacien à (0,0)
    for (int i(0); i < Nx; ++i){
        for (int j(0); j < Ny; ++j){
            for (int k(0); k < Nz; ++k){
                collection3D[i][j][k].set_potentiel((-v_infini/2) * Z[k],(-v_infini/2) * Y[j]);
                collection3D[i][j][k].set_laplacien(0,0);
                if((Z[k] <= M.altitude(X[i], Y[j]) and (i != 0) and (i != Nx - 1) and (j != 0) and (Ny - 1 != j))){
                    collection3D[i][j][k].set_potentiel(0,0);

                }
            }
        }
    }

}




void ChampPotentiels:: calcule_laplaciens(){ // fonctionne :) mais le compilateur considère la boucle comme infini => segmentation fault
    for (int i(1); i < Nx-1; ++i){
        for (int j(1); j < Ny-1; ++j){
            for (int k(1); k < Nz-1; ++k){
                collection3D[i][j][k].set_laplacien((collection3D[i-1][j][k].get_potentiel()
                                                     + collection3D[i][j-1][k].get_potentiel()+ collection3D[i][j][k-1].get_potentiel()
                                                     - (6 * collection3D[i][j][k].get_potentiel()) + collection3D[i+1][j][k].get_potentiel()
                                                     + collection3D[i][j+1][k].get_potentiel() + collection3D[i][j][k+1].get_potentiel()).get_x(),
                                                    (collection3D[i-1][j][k].get_potentiel()
                                                     + collection3D[i][j-1][k].get_potentiel()
                                                     + collection3D[i][j][k-1].get_potentiel() - (6 * collection3D[i][j][k].get_potentiel())
                                                     + collection3D[i+1][j][k].get_potentiel() + collection3D[i][j+1][k].get_potentiel()
                                                     + collection3D[i][j][k+1].get_potentiel()).get_y());
               // cout << i <<" "<< j <<" "<< k <<": "<< collection3D[i][j][k].get_laplacien() << endl;
                // j'affiche pour savoir jusqu'où va la boucle, pour l'instant:
                // max i =1 max j= 28 et max k = 29(normal puisque Nz = 30)

            }

        }
       // cout << "ok" << endl; // ne s'affiche pas donc le pb vient de cette boucle
    }

}

// méthode affichage des vecteurs potentiels et des laplaciens
void ChampPotentiels:: affichage(){// fonctionne par contre ça ne correspond pas à leur affichage
    for (int i(0); i < Nx; ++i){
        for (int j(0); j < Ny; ++j){
            for (int k(0); k < Nz; ++k){
                cout << i<< " " << j << " " << k << " ";
                cout << collection3D[i][j][k].get_potentiel() <<endl;
            }
        }
    }
}

void ChampPotentiels::set_lambda(double x) {
    lambda=x;
    cout<<lambda<<endl;
    cout<<x<<endl;

}