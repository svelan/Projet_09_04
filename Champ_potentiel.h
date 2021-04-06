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
public:
    //La taille du tableau a été initialisé et c'est constructeur par défaut et par valeur;

    ChampPotentiels(int Nx=0,int Ny=0, int Nz=0, double lambda=0.0)
            :Nx(Nx),Ny(Ny),Nz(Nz),lambda(lambda),collection3D(Nx, std::vector<std::vector<Potentiel>>(Ny, std::vector<Potentiel>(Nz))){
    set_xyz();
    }


    void initialise(double const& v_inf, Montagne const& everest);
    void calcule_laplaciens();
    void set_xyz();
    void set_value_x(double x, int i);
    void set_value_y(double y, int j);
    void set_value_z(double z, int k);
    void set_lambda(double x);
    void affichage();
};