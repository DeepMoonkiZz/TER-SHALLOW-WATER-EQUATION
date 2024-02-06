#include "mod_schema.h"

#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>


using namespace std;

Schema::Schema()
{}


void Schema::Initialize_constant()
{
    ifstream fichier;
    fichier.open("parametre.txt");

    fichier >> _xmin;
    fichier >> _xmax;
    fichier >> _nx;
    _dx = (_xmax-_xmin)/(_nx-1);

    fichier >> _tmax;

    fichier >> _x_barrage;
    fichier >> _h_gauche;
    fichier >> _h_droite;

    _g = 9.81;
}


void Schema::Initialize_vector()
{
    _t = 0;

    _Flux.resize(_nx+1);

    _U.resize(_nx);

    for (int i=0; i < _nx; i++) {
        if (_xmin + i*_dx <= _x_barrage) {
            _U[i] = make_pair(_h_gauche, 0);
        }
        else {
            _U[i] = make_pair(_h_droite, 0);
        }
    }
}

void Schema::Update()
{
    // Fonction qui calcul la valeur bmax
    this->Update_B_max();
    // Calculer la condition CFL
    _dt = (0.9*_dx)/_bmax;
    // Calcul du nouveau t
    _t += _dt;
    // Calcul des flux
    this->Update_Flux();
    // Calcul de u(n+1)
    this->Update_U();
}


void Schema::Update_B_max()
{
    double u, h;
    _bmax = 0;
    for (int i=0; i < _nx; i++) {
        h = _U[i].first;
        u = _U[i].second/h;
        if (abs(u+sqrt(h*_g)) > _bmax) {
            _bmax = abs(u+sqrt(h*_g));
        }
        if (abs(u-sqrt(h*_g)) > _bmax) {
            _bmax = abs(u-sqrt(h*_g));
        }
    }
}


void Schema::Update_Flux()
{
    pair<double,double> U_droite, U_gauche;
    U_droite.first = _h_droite, U_droite.second = 0;
    U_gauche.first = _h_gauche, U_gauche.second = 0;

    _Flux[0] = this->Flux(_U[0], _U[0]);
    _Flux[_nx] = this->Flux(_U[_nx-1], _U[_nx-1]);

    for (int i=0; i < _nx-1; i++) {
        _Flux[i+1] = this->Flux(_U[i+1], _U[i]);
    }   
}


pair<double,double> Schema::Fonction_F(pair<double,double> Ui)
{
    pair<double,double> FUi;

    FUi.first = Ui.second;
    FUi.second = Ui.second*Ui.second/Ui.first + 0.5*_g*Ui.first*Ui.first;

    return FUi;
}



pair<double,double> Schema::Flux(pair<double,double> Up, pair<double,double> Um)
{
    pair<double,double> F, FUp, FUm;

    FUp = this->Fonction_F(Up);
    FUm = this->Fonction_F(Um);

    F.first = 0.5*(FUp.first + FUm.first) - 0.5*_bmax*(Up.first - Um.first);
    F.second = 0.5*(FUp.second + FUm.second) - 0.5*_bmax*(Up.second - Um.second);

    return F;
}


void Schema::Update_U()
{
    for (int i=0; i < _nx; i++) {
        _U[i].first = _U[i].first - (_dt/_dx)*(_Flux[i+1].first - _Flux[i].first);
        _U[i].second = _U[i].second - (_dt/_dx)*(_Flux[i+1].second - _Flux[i].second);
    }
}


void Schema::Save_Solution()
{
    char buffer[100];
    std::snprintf(buffer, sizeof(buffer), "/home/segal/Documents/MatMeca/TER-2A-bis/Solutions/valid_1D/valid_1D_%i.dat", _k);
    _k++;

    // Écriture du fichier
    std::ofstream file_out;
    file_out.open(buffer); // Utilisation de buffer pour ouvrir le fichier

    if (file_out.is_open()) {
        for (int i = 0; i < _nx; i++) {
            file_out << _xmin + i * _dx << " " << _U[i].first << " " << _U[i].second/_U[i].first << " " << this->GetT() << std::endl;
        }
        file_out.close();
        //std::cout << "Les données ont été écrites dans le fichier : " << buffer << std::endl;
    } else {
        std::cout << "Impossible d'ouvrir le fichier : " << buffer << std::endl;
    }
}