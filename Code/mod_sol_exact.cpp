#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "mod_sol_exact.h"

using namespace std;


Exact::Exact()
{}


double Exact::f(double h) {
    return 2*(sqrt(_g*_h_gauche)-sqrt(_g*h)) + (_h_droite-h)*sqrt(_g*(h+_h_droite)/(2*h*_h_droite));
}


double Exact::false_position(double a, double b, double tol, int max_iter) {
    if (this->f(a) * this->f(b) >= 0) {
        throw std::invalid_argument("il faut que f(a) et f(b) soit de signe opposé");
    }


    for (int i = 0; i < max_iter; i++) {
        double c = (b * this->f(a) - a * this->f(b)) / (this->f(a) - this->f(b));

        if (std::abs(this->f(c)) < tol) {
            return c;
        }

        if (this->f(c) * this->f(a) < 0) {
            b = c;
        } else {
            a = c;
        }
    }

    throw std::runtime_error("Max d'iteration atteint");
}

void Exact::Initialisation() {
    ifstream fichier;
    fichier.open("parametre.txt");

    int nx;
    double tmax, x_barrage;

    fichier >> _xmin;
    fichier >> _xmax;
    fichier >> nx;
    _dx = (_xmax-_xmin)/(nx-1);

    fichier >> tmax;

    fichier >> x_barrage;
    fichier >> _h_gauche;
    fichier >> _h_droite;


    _g = 9.81;
    
    double tolerance = 1e-6;
    int max_iterations = 10000;
    double dx;
    double sigma;

    _n = 0;
    _he = this->false_position(0.0001, _xmax, tolerance, max_iterations);
    _ue = 2*(sqrt(_g*_h_gauche)-sqrt(_g*_he));
    _lambdag = -sqrt(_g*_h_gauche);
    _lambdae = _ue - sqrt(_g*_he);
    _sigma = (-_he*_ue)/(_h_droite-_he);
}



void Exact::Update(double t)
{
    std::string name_file ;
    std::ofstream mon_flux ;
    name_file = "/home/segal/Documents/MatMeca/TER-2A-bis/Solutions/exact_1D/exact_1D_" + std::to_string(_n) + ".dat";
    mon_flux.open(name_file);

    double x;    
    double h;
    double u;
    
    x = _xmin;
    while(x<=_xmax)
    {
        if (x<=_lambdag*t)
        {
            mon_flux << x <<  " " << _h_gauche << " " << 0 << std::endl;
        }
        else if (x>_lambdag*t && x<=_lambdae*t)
        {
            h = (1./(9*_g))*pow(2*sqrt(_g*_h_gauche)-(x/t),2);
            u = 2*(sqrt(_g*_h_gauche)-sqrt(_g*h));
            mon_flux << x <<  " " << h << " " << u << std::endl;
        }
        else if (x>_lambdae*t && x<=_sigma*t)
        {
            mon_flux << x <<  " " << _he << " " << _ue << std::endl;
        }
        else if (x>_sigma*t)
        {
            mon_flux << x <<  " " << _h_droite << " " << 0 << std::endl;
        }
        x += _dx;
    }
    _n += 1;
    mon_flux.close();
}