#pragma once
#include "Vecteur2D.h"
#include "Montagne.h"
#include <vector>

class Potentiel{
private:
    Vecteur2D potentiel;
    Vecteur2D laplacien;
public:
    Vecteur2D get_potentiel() const;
    Vecteur2D get_laplacien() const;
    void set_potentiel(double,double);
    void set_laplacien(double,double);
};

class ChampPotentiels{
private:
    std::vector<std::vector<std::vector<Potentiel>>> collection3D;
    unsigned int Nx;
    unsigned int Ny;
    unsigned int Nz;
    double lambda;
    std::vector<double> X,Y,Z;
    static double epsilon;
public:
    //La taille du tableau a été initialisé et c'est constructeur par défaut et par valeur;

    ChampPotentiels(int Nx=0,int Ny=0, int Nz=0, double lambda=0.0)
            :Nx(Nx),Ny(Ny),Nz(Nz),lambda(lambda),collection3D(Nx, std::vector<std::vector<Potentiel>>(Ny, std::vector<Potentiel>(Nz))){
    set_xyz();
    }


    void initialise(double const& v_inf, Montagne const& everest);
    void calcule_laplaciens();
    void set_xyz();
    void affichage();
    void lapla_affichage();


    // méthode erreur() qui renvoie la somme des carrés des normes des vecteurs laplacien ;
    double erreur();
    //une méthode iteration() qui applique l'équation (6) du complément mathématique en tout point
    void iteration();
    //une méthode resolution() qui répète l'itération précédente tant que l'erreur est plus grande qu'un seuil donné et le nombre d'itérations plus petit qu'un maximum donné.
    void resolution(double seuil, unsigned int max_iterations, bool verbeuse);
    //une méthode vitesse() qui prend trois paramètres i, j et k et retourne un tableau de trois double qui sont les coordonnées de la vitesse du vent en (xi, yj, zk)
    std::vector<double> vitesse(unsigned int i, unsigned int j, unsigned int k);

    void affiche_total();

};