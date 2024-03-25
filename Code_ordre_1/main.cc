#include "mod_schema.h"
#include "mod_sol_exact.h"

#include <iostream>
#include <vector>

using namespace std;


int main(int argc,char** argv)
{
    double t, tmax;
    Schema* schema(0);
    Exact* exact(0);

    schema = new Schema();
    exact = new Exact();


    // ----------------------------------------------------------- //

    // Solution numérique

    // ----------------------------------------------------------- //

    schema->Initialize_constant();
    schema->Initialize_vector(); 

    exact->Initialisation();

    tmax = schema->GetTmax();

    while(schema->GetT() < tmax)
    {
        schema->Update();
        schema->Save_Solution();
        exact->Update(schema->GetT());
    }

    
    return 0;
}