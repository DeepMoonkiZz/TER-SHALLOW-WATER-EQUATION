#include <iostream>

using namespace std;

class Exact
{
    private:

        // Acceleration grav 
        double _g;        
        // x minimum
        double _xmin;
        // x maximum
        double _xmax;
        // hauteur a gauche du barrage
        double _h_gauche;
        // hauteur a droite du barrage
        double _h_droite;
        // 
        double _dx;
        //
        int _n;
        // 
        double _lambdag, _lambdae;
        // 
        double _he, _ue;
        // 
        double _sigma;
        

    public:

    // ----------------------------------------
    //                FONCTION
    // ----------------------------------------

        // Constructeur de schema
        Exact();
        // Initialisation 
        void Initialisation();
        // Fonction dont on cherche le 0
        double f(double h);
        // Fausse position
        double false_position(double a, double b, double tol, int max_iter);
        // Update a chaque iteration en temps
        void Update(double t);    
};