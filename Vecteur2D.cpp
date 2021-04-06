#include <iostream>
#include <cmath>
#include "Vecteur2D.h"
using namespace std;

//prototypes de fonctions
void Vecteur2D::affiche(ostream& sortie) const{
    sortie<<x<<" "<<y<<" ";
}

void Vecteur2D::set_coord(double x1, double y1){
    x=x1;
    y=y1;
}


bool Vecteur2D::compare(Vecteur2D v1) const{
    return((v1.x==x)&&(v1.y==y));
}

Vecteur2D Vecteur2D::addition(Vecteur2D autre){
    Vecteur2D v1;
    v1.set_coord(x+autre.x,y+autre.y);
    return v1;
}

Vecteur2D Vecteur2D::soustraction(Vecteur2D autre){
    Vecteur2D v1;
    v1.set_coord(x-autre.x, y-autre.y);
    return v1;
}

double Vecteur2D::prod_scal(Vecteur2D autre) const{
    return x*autre.x+y*autre.y;
}


Vecteur2D Vecteur2D::oppose() const{
    Vecteur2D v1;
    v1.set_coord(-x,-y);
    return v1;
}

Vecteur2D Vecteur2D::mult(double scal) const{
    Vecteur2D v1;
    v1.set_coord((scal*x), scal*y);

    return v1;
}

Vecteur2D Vecteur2D::unitaire() {
    Vecteur2D v1;
    v1.set_coord(x/(sqrt(x*x+y*y)), y/(sqrt(x*x+y*y)));
    return v1;
}


//OPERATEUR D'AFFICHAGE

ostream& operator<<(ostream& sortie, const Vecteur2D& v){
    v.affiche(sortie);
    return sortie;
}

//OPERATEURS D'EGALITE

bool Vecteur2D:: operator==(const Vecteur2D& v){
    return ((x==v.x) && (y==v.y));
}

//OPERATEURS += ; -=  ; + ; - ; produit par un scalaire ; produit scalaire ; norme
//opposé
Vecteur2D& Vecteur2D::operator+=(const Vecteur2D& v1){
    x+=v1.x;
    y+=v1.y;
    return *this;
}
Vecteur2D& Vecteur2D::operator-=(const Vecteur2D& v1){
    x-=v1.x;
    y-=v1.y;
    return *this;
}
Vecteur2D& Vecteur2D::operator*=(double a){
    x=a*x;
    y=a*y;
    return *this;
}
const Vecteur2D operator*(double a, const Vecteur2D& v){
    Vecteur2D u;
    u+=v;
    u*=a;
    return u;
}
const Vecteur2D operator*(Vecteur2D v, double a){
    Vecteur2D u;
    u+=v;
    u*=a;
    return u;
}
const Vecteur2D operator+(Vecteur2D u, const Vecteur2D& v){
    Vecteur2D extra;
    extra+=v;
    extra+=u;
    return extra;
}

const Vecteur2D operator-(Vecteur2D u, const Vecteur2D& v){
    Vecteur2D extra(u);
    extra-=v;
    return extra;
}

double operator*(const Vecteur2D& un, const Vecteur2D& deux){
    return un.get_x()*deux.get_x()+un.get_y()*deux.get_y();
}

double Vecteur2D:: norme() const{
    return sqrt(*this * *this);
}

double Vecteur2D:: norme2() const{
    return norme()*norme();
}

const Vecteur2D operator-(Vecteur2D a) {
    Vecteur2D u;
    u.set_coord(-a.get_x(),-a.get_y());
    return u;
}
const Vecteur2D operator~(Vecteur2D a){
    return a.unitaire();
}
