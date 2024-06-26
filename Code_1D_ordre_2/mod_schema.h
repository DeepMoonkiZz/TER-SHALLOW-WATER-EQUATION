#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

using namespace std;

class Schema
{
    private:

    // ----------------------------------------
    //                CONSTANT
    // ----------------------------------------

        // x minimum
        double _xmin;
        // x maximum
        double _xmax;
        // pas d'espace en x
        double _dx;
        // Nombre de maille en x
        int _nx;

        // valeur de tmax
        double _tmax;

        // position du barrage
        double _x_barrage;
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
        vector<pair<double,double>> _U;
        // Vecteur U prime
        vector<pair<double,double>> _U_prime;
        // Vecteur U second
        vector<pair<double,double>> _U_second;
        // Vecteur flux
        vector<pair<double,double>> _Flux;
        // Vecteur flux prime
        vector<pair<double,double>> _Flux_prime;
        // Vecteur u ordre 2
        vector<pair<double,double>> _UO2;
        // Vecteur vérification de U
        vector<int> _Verification_U;
        // Vecteur vérification de U prime
        vector<int> _Verification_U_prime;
        // Vecteur vérification du flux
        vector<int> _Verification_Flux;
        // Vecteur vérification de flux prime
        vector<int> _Verification_Flux_prime;
        

    public:

    // ----------------------------------------
    //                FONCTION
    // ----------------------------------------

        // Constructeur de schema
        Schema();
        // Initialisation de la classe schema
        void Initialize_constant();
        // Initialisation de la classe schema
        void Initialize_vector();
        // Update la temperature un dt plus tard
        void Update();
        // Update la valeur de b max
        void Update_B_max();
        // Update la valeur des flux
        vector<pair<double,double>> Update_Flux(vector<pair<double,double>> U, vector<int> Verification_Flux);
        // Update la valeur des flux
        vector<pair<double,double>> Update_Flux_O2(vector<pair<double,double>> U);
        // Update la valeur de u
        vector<pair<double,double>> Update_U(vector<pair<double,double>> U, vector<pair<double,double>> F, vector<int> Verification_U);
        // Update la valeur de uo2
        vector<pair<double,double>> Update_UO2(vector<pair<double,double>> U, vector<pair<double,double>> U_second, vector<int> Verification_U);
        // Update la valeur des vecteur verificaiton
        void Update_verification();
        // Applique la fonction F à la maille i
        pair<double,double> Fonction_F(pair<double,double> Ui);
        // Calcul le flux entre i et i+1
        pair<double,double> Flux(pair<double,double> Up, pair<double,double> Um);

    // ----------------------------------------
    //          RECUPERER LES VARIABLES
    // ----------------------------------------

        // recuperer t
        double GetT() {return _t;};
        // recuperer t maximum
        double GetTmax() {return _tmax;};
        // recuperer t maximum
        int GetNX() {return _nx;};
        // recuperer U
        vector<pair<double,double>> GetU() {return _U;};
        // Sauvegarder la solution
        void Save_Solution();
};