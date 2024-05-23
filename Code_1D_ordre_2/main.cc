#include "mod_schema.h"
#include "mod_sol_exact.h"

#include <iostream>
#include <vector>
#include <fstream>


using namespace std;


pair<double, double> Compute_error(vector<pair<double, double>> U, vector<pair<double, double>> U_exact) 
{
    pair<double, double> error;
    double error_u(0), error_v(0);
    for(int i=0; i<U.size(); i++) {
        error_u += pow((U[i].first - U_exact[i].first), 2);
        error_v += pow((U[i].second/U[i].first - U_exact[i].second), 2);
    }
    error.first = error_u/U.size();
    error.second = error_v/U.size();

    return error;
}


int main(int argc,char** argv)
{
    double t, tmax;
    pair<double, double> error;
    Schema* schema(0);
    Exact* exact(0);

    schema = new Schema();
    exact = new Exact();

    // 
    
    char buffer[100];
    snprintf(buffer, sizeof(buffer), "/home/segal/Documents/MatMeca/TER-SHALLOW-WATER-EQUATION/Solutions/error_1D/error_%i.dat", schema->GetNX());



    // ----------------------------------------------------------- //

    // Solution numérique

    // ----------------------------------------------------------- //

    schema->Initialize_constant();
    schema->Initialize_vector(); 

    exact->Initialisation();

    tmax = schema->GetTmax();

    // Écriture du fichier
    ofstream file_out;
    file_out.open(buffer); // Utilisation de buffer pour ouvrir le fichier

    while(schema->GetT() < tmax)
    {
        schema->Update();
        schema->Save_Solution();
        exact->Update(schema->GetT());
        error = Compute_error(schema->GetU(), exact->GetU());
        file_out << schema->GetT() << " " << error.first << " " << error.second << std::endl;
    }
    file_out.close();

    
    return 0;
}


