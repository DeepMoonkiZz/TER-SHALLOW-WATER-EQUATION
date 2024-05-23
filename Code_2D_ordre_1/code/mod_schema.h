#ifndef _MOD_SCHEMA_H

#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include "Mesh2D.h"

using namespace std;

class Schema
{
    protected:

    // ----------------------------------------
    //                CONSTANT
    // ----------------------------------------

        // // x minimum
        // double _xmin;
        // // x maximum
        // double _xmax;
        // pas d'espace delta
        double _delta;
        // Nombre de maille
        int _N;

        // valeur de tmax
        double _tmax;

        // position du barrage
        //double _x_barrage;
        // hauteur a gauche du barrage
        double _h_gauche;
        // hauteur a droite du barrage
        double _h_droite;

        // accélération gravitationnel
        double _g;
        // compteur d'itération
        int _k;

    // ----------------------------------------
    //                VARIABLE
    // ----------------------------------------

        // Valeur de t
        double _t;
        // Valeur de delta t 
        double _dt;
        // Valeur de bmax
        double _bmax;
        // Vecteur u
        vector<vector<double>> _U;
        // Vecteur flux
        //vector<vector<double>> _Flux;
        const Mesh2D* _msh;
        

    public:

    // ----------------------------------------
    //                FONCTION
    // ----------------------------------------

        // Constructeur de schema
        Schema(Mesh2D* mesh);
        //Destructeur du schéma
        //virtual ~Schema();
        // Initialisation de la classe schema
        void Initialize_constant();
        // Initialisation de la classe schema
        void Initialize_vector();
        // Update la temperature un dt plus tard
        void Update();
        // Update la valeur de b max
        void Update_B_max();
        // Update la valeur des flux ordre 1
        // Fonction du fond
        double B(double x, double y);
        // Vecteur source 
        vector<double> S(double x, double y);
        //virtual void Update_Flux();
        // Update la valeur de u
        void Update_U();
        // Applique la fonction F à la maille i
        vector<vector<double>> Fonction_F(vector<double> Ui);
        // Calcul le flux entre i et i+1
        vector<double> Flux(vector<double> Uk, vector<double>Ul, vector<double> ne);
        // fonction critere methode ...
        //virtual int critere(pair<double,double> Ui) = 0;
    // ----------------------------------------
    //          RECUPERER LES VARIABLES
    // ----------------------------------------

        // recuperer t
        double GetT() {return _t;};
        // recuperer t maximum
        double GetTmax() {return _tmax;};
        // recuperer U
        vector<vector<double>> GetU() {return _U;};
        // Sauvegarder la solution
        void Save_Solution_h();
        void Save_Solution_U();
        void Save_Solution_V();
};


#define _MOD_SCHEMA_H
#endif