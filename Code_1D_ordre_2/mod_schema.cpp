#include "mod_schema.h"

#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <thread>
#include <chrono>


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

    _Verification_U.resize(_nx, 1);
    _Verification_Flux.resize(_nx+1, 1);
    _Verification_U_prime.resize(_nx, 1);
    _Verification_Flux_prime.resize(_nx+1, 1);

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

    // Calcul des flux
    _Flux = this->Update_Flux_O2(_U);
    // Calcul de U prime
    _U_prime = this->Update_U(_U, _Flux, _Verification_U_prime);
    // Calcul des flux prime
    _Flux_prime = this->Update_Flux_O2(_U_prime);
    // Calcul de U second
    _U_second = this->Update_U(_U_prime, _Flux_prime, _Verification_U);
    // Calcul de U ordre 2
    _UO2 = this->Update_UO2(_U, _U_second, _Verification_U);
    // Update les vecteurs vérification
    this->Update_verification();

    // Boucle sur l'ordre 1 si critère limiteur non atteint
    while (!all_of(_Verification_U.begin(), _Verification_U.end(), [](int i) { return i == 0; })) {
        // Calcul des flux
        _Flux = this->Update_Flux(_U, _Verification_Flux);
        // Calcul de U prime
        _U_prime = this->Update_U(_U, _Flux, _Verification_U_prime);
        // Calcul des flux prime
        _Flux_prime = this->Update_Flux(_U_prime, _Verification_Flux_prime);
        // Calcul de U second
        _U_second = this->Update_U(_U_prime, _Flux_prime, _Verification_U);
        // Calcul de U ordre 2
        _UO2 = this->Update_UO2(_U, _U_second, _Verification_U);
        // Update les vecteurs vérification
        this->Update_verification();
    }

    _U = _UO2;

        cout << "Verif U :";
    for (int i=0; i<_nx; i++) {
        cout << _U[i].first << " ";
    }
    cout << endl;

    _Verification_U.assign(_nx, 1);
    _Verification_Flux.assign(_nx+1, 1);
    _Verification_U_prime.assign(_nx, 1);
    _Verification_Flux_prime.assign(_nx+1, 1);
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


vector<pair<double,double>> Schema::Update_Flux_O2(vector<pair<double,double>> U)
{
    vector<pair<double,double>> F(_nx+1, {0, 0});
    pair<double,double> U_droite, U_gauche, Um, Up, Upp;

    Um.first = (2.*U[0].first) - U[1].first, Um.second = (2.*U[0].second) - U[1].second;
    Up.first = (2.*U[_nx-1].first) - U[_nx-2].first, Up.second = (2.*U[_nx-1].second) - U[_nx-2].second;
    Upp.first = (3.*U[_nx-1].first) - (2.*U[_nx-2].first), Upp.second = (3.*U[_nx-1].second) - (2.*U[_nx-2].second);



    for (int i=0; i<_nx+1; i++) {
        if (i==0) {
            U_gauche.first = (3.*U[0].first - U[1].first)/2., U_gauche.second = (3.*U[0].second - U[1].second)/2.;
            U_droite.first = (Um.first + U[0].first)/2., U_droite.second = (Um.second + U[0].second)/2.;
            F[0] = this->Flux(U[0], U[0]);
        }
        else if (i==_nx) {
            U_gauche.first = (3.*Up.first - Upp.first)/2., U_gauche.second = (3.*Up.second - Upp.second)/2.;
            U_droite.first = (U[_nx-1].first + Up.first)/2., U_droite.second = (U[_nx-1].second + Up.second)/2.;
            F[_nx] = this->Flux(U[_nx-1], U[_nx-1]);
        }
        else if (i==_nx-1) {
            U_gauche.first = (3.*U[_nx-1].first - Up.first)/2., U_gauche.second = (3.*U[_nx-1].second - Up.second)/2.;
            U_droite.first = (U[_nx-2].first + U[_nx-1].first)/2., U_droite.second = (U[_nx-2].second + U[_nx-1].second)/2.;
            F[_nx-1] = this->Flux(U_gauche, U_droite);
        }
        else {
            U_gauche.first = (3.*U[i].first - U[i+1].first)/2., U_gauche.second = (3.*U[i].second - U[i+1].second)/2.;
            U_droite.first = (U[i-1].first + U[i].first)/2., U_droite.second = (U[i-1].second + U[i].second)/2.;
            F[i] = this->Flux(U_gauche, U_droite);
        }
    }

    return F;
}


vector<pair<double,double>> Schema::Update_Flux(vector<pair<double,double>> U, vector<int> Verification_Flux)
{
    vector<pair<double,double>> F(_nx+1, {0, 0});

    for (int i=0; i<_nx+1; i++) {
        if (Verification_Flux[i]) {
            if (i==0) {
                F[0] = this->Flux(U[0], U[0]);
            }
            else if (i==_nx) {
                F[_nx] = this->Flux(U[_nx-1], U[_nx-1]);
            }
            else {
                F[i] = this->Flux(U[i], U[i-1]);
            }
        }
    }

    return F;
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


vector<pair<double,double>> Schema::Update_U(vector<pair<double,double>> U, vector<pair<double,double>> F, vector<int> Verification_U)
{
    vector<pair<double,double>> U_update(_nx, {0, 0});

    for (int i=0; i<_nx; i++) {
        if (Verification_U[i]) {
            U_update[i].first = U[i].first - (_dt/_dx)*(F[i+1].first - F[i].first);
            U_update[i].second = U[i].second - (_dt/_dx)*(F[i+1].second - F[i].second);
        }
    }

    return U_update;
}


vector<pair<double,double>> Schema::Update_UO2(vector<pair<double,double>> U, vector<pair<double,double>> U_second, vector<int> Verification_U)
{
    vector<pair<double,double>> UO2(_nx);

    for (int i=0; i<_nx; i++) {
        if (Verification_U[i]) {
            UO2[i].first = (U[i].first + U_second[i].first)/2.;
            UO2[i].second = (U[i].second + U_second[i].second)/2.;
        }
        else {
            UO2[i].first = _UO2[i].first;
            UO2[i].second = _UO2[i].second;
        }
    }

    return UO2;
}


void Schema::Update_verification()
{
    _Verification_U.assign(_nx, 0);
    _Verification_Flux.assign(_nx+1, 0);
    _Verification_U_prime.assign(_nx, 0);
    _Verification_Flux_prime.assign(_nx+1, 0);

    for (int i=0; i<_nx; i++) {
        if (_UO2[i].first<0) {
            _Verification_U[max(0, i-2)] = 1;
            _Verification_U[max(0, i-1)] = 1;
            _Verification_U[i] = 1;
            _Verification_U[min(_nx-1, i+1)] = 1;
            _Verification_U[min(_nx-1, i+2)] = 1;
            _Verification_Flux_prime[i] = 1;
            _Verification_Flux_prime[i+1] = 1;
            _Verification_U_prime[max(0, i-1)] = 1;
            _Verification_U_prime[i] = 1;
            _Verification_U_prime[min(_nx-1, i+1)] = 1;
            _Verification_Flux[max(0, i-1)] = 1;
            _Verification_Flux[i] = 1;
            _Verification_Flux[i+1] = 1;
            _Verification_Flux[min(_nx, i+2)] = 1;
        }
    }
}


void Schema::Save_Solution()
{
    char buffer[100];
    std::snprintf(buffer, sizeof(buffer), "Solutions/sol_1D_ordre_2/valid_1D/valid_1D_%i.dat", _k);
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