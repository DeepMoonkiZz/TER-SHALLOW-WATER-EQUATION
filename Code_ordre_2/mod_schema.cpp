#include "mod_schema.h"

#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>


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

    _Flux_prime.resize(_nx+1);

    _U.resize(_nx);

    _U_prime.resize(_nx);

    _U_second.resize(_nx);

    _UO2.resize(_nx);

    _Verification.resize(_nx, 1);

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
    _dt = (0.9*_dx)/(2.*_bmax);
    // Calcul du nouveau t
    _t += _dt;
    // boucle de verification
    while (!all_of(_Verification.begin(), _Verification.end(), [](int i) { return i == 0; })) {
        this->Update_UO2();
        this->Update_U();
    }
    _Verification.assign(_nx, 1);
}


void Schema::Update_U()
{
    for (int i=0; i < _nx; i++) {
        if (_UO2[i].first>=0) {
            _U[i].first = _UO2[i].first;
            _U[i].second = _UO2[i].second;  
            _Verification[i] = 0;
        }
    }
}


void Schema::Update_UO2()
{
    for (int i=0; i<_nx; i++) {
        if (_Verification[i]==1) {
            for (int j=max(0, i-2); j<min(_nx, j+3); j++) {
                this->Compute_Flux(j);
                this->Compute_U_prime(j);
                this->Compute_Flux_prime(j);
                this->Compute_U_second(j);
                this->Compute_UO2(j);
            }
        }
    }
}


void Schema::Compute_UO2(int j)
{
    _UO2[j].first = _U[j].first + _U_second[j].first;
    _UO2[j].second = _U[j].second + _U_second[j].second;
}


void Schema::Compute_U_prime(int j)
{
    for (int i=max(0, j-1); i<min(_nx, j+2); i++) {
        _U_prime[j].first = _U[j].first - (_dt/_dx)*(_Flux[j+1].first - _Flux[j].first);
        _U_prime[j].second = _U[j].second - (_dt/_dx)*(_Flux[j+1].second - _Flux[j].second);
    }
}


void Schema::Compute_U_second(int j)
{
    _U_second[j].first = _U_prime[j].first - (_dt/_dx)*(_Flux_prime[j+1].first - _Flux_prime[j].first);
    _U_second[j].second = _U_prime[j].second - (_dt/_dx)*(_Flux_prime[j+1].second - _Flux_prime[j].second);
}


void Schema::Compute_Flux(int j)
{
    pair<double,double> U_droite, U_gauche; 

    if (j==0 || j==1) {
        U_gauche.first = 2.*_U[0].first - _U[1].first, U_gauche.second = 2.*_U[0].second - _U[1].second;
        _Flux[0] = this->Flux(_U[0], U_gauche);
    }
    else if (j==_nx-1 || j==_nx-2) {
        U_droite.first = 2.*_U[_nx-1].first - _U[_nx-2].first, U_droite.second = 2.*_U[_nx-1].second - _U[_nx-2].second;
        _Flux[_nx] = this->Flux(U_droite, _U[_nx-1]);
    }
    for (int i=max(0, j-2); i<min(_nx-1, j+2); i++) {
        U_gauche.first = (_U[i].first + _U[i+1].first)/2., U_gauche.second = (_U[i].second + _U[i+1].second)/2.;
        U_droite.first = (3.*_U[i].first - _U[i+1].first)/2., U_droite.second = (3.*_U[i].second - _U[i+1].second)/2.;
        _Flux[i+1] = this->Flux(U_gauche, U_droite);
    }
}


void Schema::Compute_Flux_prime(int j)
{
    pair<double,double> U_droite, U_gauche; 

    if (j==0) {
        U_gauche.first = 2.*_U_prime[0].first - _U_prime[1].first, U_gauche.second = 2.*_U_prime[0].second - _U_prime[1].second;
        _Flux_prime[0] = this->Flux(_U_prime[0], U_gauche);
    }
    else if (j==_nx-1) {
        U_droite.first = 2.*_U_prime[_nx-1].first - _U_prime[_nx-2].first, U_droite.second = 2.*_U_prime[_nx-1].second - _U_prime[_nx-2].second;
        _Flux_prime[_nx] = this->Flux(U_droite, _U_prime[_nx-1]);
    }
    for (int i=max(0, j-1); i<min(_nx-1, j+1); i++) {
        U_gauche.first = (_U_prime[i].first + _U_prime[i+1].first)/2., U_gauche.second = (_U_prime[i].second + _U_prime[i+1].second)/2.;
        U_droite.first = (3.*_U_prime[i].first - _U_prime[i+1].first)/2., U_droite.second = (3.*_U_prime[i].second - _U_prime[i+1].second)/2.;
        _Flux_prime[i+1] = this->Flux(U_gauche, U_droite);
    }
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


pair<double,double> Schema::Fonction_F(pair<double,double> Ui)
{
    pair<double,double> FUi;

    FUi.first = Ui.second;
    FUi.second = Ui.second*Ui.second/Ui.first + 0.5*_g*Ui.first*Ui.first;

    return FUi;
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


void Schema::Save_Solution()
{
    char buffer[100];
    std::snprintf(buffer, sizeof(buffer), "/home/segal/Documents/MatMeca/TER-SHALLOW-WATER-EQUATION/Solutions/valid_1D/valid_1D_%i.dat", _k);
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