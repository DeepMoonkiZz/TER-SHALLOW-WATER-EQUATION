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
        // Vecteur flux
        vector<pair<double,double>> _Flux;
        

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
        void Update_Flux();
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