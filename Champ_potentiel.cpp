#include "Champ_potentiel.h"
#include "Vecteur2D.h"
#include "Montagne.h"
#include <vector>
using namespace std;
//MODIF
double ChampPotentiels::epsilon(0.1);

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
        Y.push_back(j*lambda - (lambda*(Ny-1)/2));
    }
    for(double k(0); k < Nz; ++k){
        Z.push_back(k*lambda);
    }


}

// méthode initialisant le potentiel vecteur en tout point (xi,yj,zk) -- fonctionne
//les valeurs affichées sont correctes !
void ChampPotentiels:: initialise(double const& v_infini, Montagne const& M) {
    //initialisation du potentiel vecteur en tout point et du laplacien à (0,0)
    for (int i(0); i < Nx; ++i){
        for (int j(0); j < Ny; ++j){
            for (int k(0); k < Nz; ++k){
                collection3D[i][j][k].set_potentiel((-v_infini/2)*Z[k],Y[j]*(v_infini/2));
                collection3D[i][j][k].set_laplacien(0,0);
                if((Z[k] <= M.altitude(X[i], Y[j])) and (i != 0) and (i != Nx - 1) and (j != 0) and (Ny - 1 != j)){
                    collection3D[i][j][k].set_potentiel(0,0);
                }
            }
        }
    }

}




void ChampPotentiels:: calcule_laplaciens(){//n'affiche pas les bonnes valeurs
  Vecteur2D v(0,0);

    for (int i(1); i < Nx-1; ++i){
        for (int j(1); j < Ny-1; ++j){
            for (int k(1); k < Nz-1; ++k){
            v=(collection3D[i-1][j][k].get_potentiel()
               + collection3D[i][j-1][k].get_potentiel()
               + collection3D[i][j][k-1].get_potentiel() - (6 * collection3D[i][j][k].get_potentiel())
               + collection3D[i+1][j][k].get_potentiel() + collection3D[i][j+1][k].get_potentiel()
               + collection3D[i][j][k+1].get_potentiel());

                collection3D[i][j][k].set_laplacien(v.get_x(),v.get_y());

                cout << i <<" "<< j <<" "<< k <<": "<< collection3D[i][j][k].get_laplacien() << endl;
                // j'affiche pour savoir jusqu'où va la boucle, pour l'instant:
                // max i =1 max j= 28 et max k = 29(normal puisque Nz = 30)

            }

        }
       // cout << "ok" << endl; // ne s'affiche pas donc le pb vient de cette boucle
    }

}

// méthode affichage des vecteurs potentiels et des laplaciens
void ChampPotentiels:: affichage(){// ceci fonctionne !! mais gnuplot n'affiche rien !
    for (int i(0); i < Nx; ++i){
        for (int j(0); j < Ny; ++j){
            for (int k(0); k < Nz; ++k){
                cout << i<< " " << j << " " << k << " ";
                cout << collection3D[i][j][k].get_potentiel() <<endl;
            }
        }
    }
}

void ChampPotentiels::lapla_affichage() {
    for (int i(0); i < Nx-1; ++i){
        for (int j(0); j < Ny-1; ++j){
            for (int k(0); k < Nz-1; ++k){
                cout << i<< " " << j << " " << k << " ";
                cout << collection3D[i][j][k].get_laplacien() <<endl;
            }
        }
    }
}


// méthode erreur() qui renvoie la somme des carrés des normes de tous les vecteurs laplacien ;
double ChampPotentiels:: erreur(){
    double retour(0.0);
    for (int i(1); i < Nx; ++i){
        for (int j(1); j < Ny; ++j){
            for (int k(1); k < Nz; ++k){
                retour += collection3D[i][j][k].get_laplacien().norme2();
            }
        }
    }
    return retour;
}

// une méthode iteration() qui applique l'équation (6) du complément mathématique en tout point (u représente le potentiel vecteur)
void ChampPotentiels:: iteration() {
    //car à chaque fois, le vecteur laplacien change -> car le potentiel change;
    calcule_laplaciens();
    Vecteur2D u(0, 0);
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            for (int k = 1; k < Nz - 1; ++k) {
                u = collection3D[i][j][k].get_potentiel() + epsilon * collection3D[i][j][k].get_laplacien();
               //ici j'ai juste modifié que c'est un set_potentiel et non un set_laplacien
               collection3D[i][j][k].set_potentiel(u.get_x(), u.get_y());
            }
        }
    }
}
//une méthode resolution() qui répète l'itération précédente tant que l'erreur est plus grande qu'un seuil donné et le nombre d'itérations plus petit qu'un maximum donné.
void ChampPotentiels:: resolution(double seuil, unsigned int max_iterations, bool verbeuse){
//SHRI: l'affichage fait quelque chose de bizarre, il modifie les données au premier pas
    // et ne les modifie plus après, donc fait les 2000 iterations, il y a un problème avec
    //la méthode iterations.

    unsigned int max(0);
cout<<1<<" "<<erreur()<<endl;
    while((erreur()>seuil)&&(max<=max_iterations)){
    iteration();
    if (verbeuse) affichage(); lapla_affichage();
    max+=1;
    cout<<max+1<<" "<<erreur()<<endl;
}

}
//une méthode vitesse() qui prend trois paramètres i, j et k et retourne un tableau de trois double qui sont les coordonnées de la vitesse du vent en (xi, yj, zk)
vector<double> ChampPotentiels:: vitesse(unsigned int i, unsigned int j, unsigned int k){
    // attention aux conditions aux bords
    double v_infini;
    vector<double> vitesse(3);
    // A_3,i,j,k= v_infini/2 * Y[j] et A_2,i,j,k= -v_infini/2 * Z[k] => conditions aux bords mais au bord on force a v_infini donc:
    double ui(v_infini), uj(v_infini), uk(v_infini); // par defaut on met les valeurs du bord
    if((i==0)or(j==0)or(k==0)or(i==29)or(j==29)or(k==29)){
        vitesse[0]=0;
        vitesse[1]=0;
        vitesse[2]=0;
    }else{

    // si pas au bord alors:
    // par l'équation 7 du complément mathématique on a :
    //ui = 1/(2*lambda)*(A_3,i,j+1,k − A3_i,j−1,k − A2_i,j,k+1 + A2_i,j,k−1)
    //uj = 1/(2*lambda)*(−A3_i+1,j,k + A3_i−1,j,k)
    //uk = 1/(2*lambda)*(A2_i+1,j,k − A2_i−1,j,k)
    ui = (1/(2*lambda))*(collection3D[i][j+1][k].get_potentiel().get_y()
    -collection3D[i][j-1][k].get_potentiel().get_y()
    -collection3D[i][j][k+1].get_potentiel().get_x()
    -collection3D[i][j][k-1].get_potentiel().get_x());

    uj = (1/(2*lambda))*(-collection3D[i+1][j][k].get_potentiel().get_y()+collection3D[i-1][j][k].get_potentiel().get_y());
    uk = (1/(2*lambda))*(-collection3D[i+1][j][k].get_potentiel().get_x()+collection3D[i-1][j][k].get_potentiel().get_x());
    vitesse[0]=ui;
    vitesse[1]=uj;
    vitesse[2]=uk;
    }

    return vitesse;
}

double ChampPotentiels::norme3D_2(vector<double> test){
    return test[0]*test[0]+test[1]*test[1]+test[2]*test[2];
}

void ChampPotentiels::affiche_total() {
    //PAS LES BONNES RéPONSES POUR UNE QUELCONQUE RAISON ?!?!

    for (int i = 1; i < Nx-1; ++i) {
        for (int j = 1; j < Ny-1; ++j) {
            for (int k = 1; k < Nz-1; ++k) {
                cout<<i<<" "<<j<<" "<<k<<" "<<collection3D[i][j][k].get_potentiel().get_x()<<" "
                    <<collection3D[i][j][k].get_potentiel().get_y()<<" "<<vitesse(i,j,k)[0]<<" "
                    <<vitesse(i,j,k)[1]<<" "<<vitesse(i,j,k)[2]<<" "
                    <<norme3D_2(vitesse(i,j,k))<<endl;
            }
        }
    }
}
