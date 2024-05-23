#include "mod_schema.h"
#include "Mesh2D.h"

#include <iostream>
#include <vector>

using namespace std;


int main(int argc,char** argv)
{
    double t, tmax;
    Schema* schema(0);
    //maillage avec séparation
    const std::vector<int> ref {10,11,12,13,14};
    const std::vector<std::string> BC {"sortie","mur","sortie","mur","barrage"};
    Mesh2D* mesh = new Mesh2D(ref,BC);
    mesh->Read_mesh("../maillage/canal.mesh");

    //maillage sans séparation
    // const std::vector<int> ref {10,11,12,13};
    // const std::vector<std::string> BC {"sortie","mur","sortie","mur"};
    //Mesh2D* mesh = new Mesh2D(ref,BC);
    //mesh->Read_mesh("../maillage/canal2.mesh");

    schema = new Schema(mesh);

    // ----------------------------------------------------------- //

    // Solution numérique

    // ----------------------------------------------------------- //

    schema->Initialize_constant();
    schema->Initialize_vector(); 

    tmax = schema->GetTmax();

    schema->Save_Solution_h();
    schema->Save_Solution_U();
    schema->Save_Solution_V();

    while(schema->GetT() < tmax)
    {
        schema->Update();
        schema->Save_Solution_h();
        schema->Save_Solution_U();
        schema->Save_Solution_V();
    }
   
    return 0;
}