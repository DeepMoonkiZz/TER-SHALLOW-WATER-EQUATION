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
        // Vecteur vérification
        vector<int> _Verification;
        

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
        void Compute_Flux(int j);
        // Update la valeur des flux
        void Compute_Flux_prime(int j);
        // Update la valeur des flux
        void Compute_U_prime(int j);
        // Update la valeur des flux
        void Compute_U_second(int j);
        // Update la valeur des flux
        void Compute_UO2(int j);
        // Update la valeur de uo2
        void Update_UO2();
        // Update la valeur de u
        void Update_U();
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