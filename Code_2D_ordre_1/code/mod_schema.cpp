#ifndef _MOD_SCHEMA_CPP

#include "mod_schema.h"

#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>


using namespace std;

Schema::Schema(Mesh2D* mesh):_msh(mesh)
{}

// Schema::~Schema()
// {}

double Schema::B(double x, double y)
{
    return x+y;
}


vector<double> Schema::S(double x, double y)
{
    vector<double> dB(3);
    double h(0.001);

    dB[0] = 0;
    dB[1] = - _dt * 9.81 * (this->B(x+h, y) - this->B(x-h, y))/(2*h);
    dB[2] = - _dt * 9.81 * (this->B(x, y+h) - this->B(x, y-h))/(2*h);

    return dB;
}


void Schema::Initialize_constant()
{
    ifstream fichier;
    fichier.open("parametre.txt");

    fichier >> _tmax;
    fichier >> _h_gauche;
    fichier >> _h_droite;
    _k = 0;
    _N = this->_msh->Get_triangles().size();

    _g = 9.81;

    double deltak;
    _delta = 100;
    const Eigen::VectorXd &Aire_maille = this->_msh->Get_triangles_area();
    const Eigen::VectorXd &Perimetre_maille= this->_msh->Get_triangles_length();
    for (int i=0; i < _N; i++)
    {
        deltak = Aire_maille[i]/Perimetre_maille[i];
        if (deltak < _delta)
        {
            _delta = deltak;
        }
    }
}


void Schema::Initialize_vector()
{
    _t = 0;

    _U.resize(_N);

    //Maillage sans séparation
    //const Eigen::Matrix<double, Eigen::Dynamic, 2> &Centre_maille = this->_msh->Get_triangles_center();

    for (int i=0; i < _N; i++) 
    {
        if (this->_msh->Get_triangles()[i].Get_reference() == 15) 
        {
            _U[i].resize(3);
            _U[i][0] = _h_gauche;
            _U[i][1] = 0;
            _U[i][2] = 0;
        }
        else 
        {
            _U[i].resize(3);
            _U[i][0] = _h_droite;
            _U[i][1] = 0;
            _U[i][2] = 0;
        }

        // double Xi = Centre_maille(i,0);

        // if(Xi <= 1)
        // {
        //     _U[i].resize(3);
        //     _U[i][0] = _h_gauche;
        //     _U[i][1] = 0;
        //     _U[i][2] = 0;
        // }
        // else 
        // {
        //     _U[i].resize(3);
        //     _U[i][0] = _h_droite;
        //     _U[i][1] = 0;
        //     _U[i][2] = 0;
        // }


    }
}


void Schema::Update()
{
    // Fonction qui calcul la valeur bmax
    this->Update_B_max();
    // Calculer la condition CFL
    _dt = (0.9*_delta)/(2*_bmax);
    // Calcul du nouveau t
    _t += _dt;
    //std::cout<< "delta " <<delta<<" dt "<<_dt<<std::endl;

    // Calcul de u(n+1)
    this->Update_U();
}



void Schema::Update_B_max()
{
    double u,v,h;
    _bmax = 0;
    for (int i=0; i < _N; i++) {
        h = _U[i][0];
        u = _U[i][1]/h;
        v = _U[i][2]/h;
        if (fabs(u+sqrt(h*_g)) > _bmax) {
            _bmax = fabs(u+sqrt(h*_g));
        }
        if (fabs(u-sqrt(h*_g)) > _bmax) {
            _bmax = fabs(u-sqrt(h*_g));
        }
        if (fabs(v+sqrt(h*_g)) > _bmax) {
            _bmax = fabs(v+sqrt(h*_g));
        }
        if (fabs(v-sqrt(h*_g)) > _bmax) {
            _bmax = fabs(v-sqrt(h*_g));
        }
    }
}


vector<vector<double>> Schema::Fonction_F(vector<double> Ui)
{
    vector<vector<double>> FUi(3);

    FUi[0].resize(2);
    FUi[0][0] = Ui[1];
    FUi[0][1] = Ui[2];

    FUi[1].resize(2);
    FUi[1][0] = (Ui[1]*Ui[1]/Ui[0]) + (_g*Ui[0]*Ui[0]/2.);
    FUi[1][1] = Ui[1]*Ui[2]/Ui[0];

    FUi[2].resize(2);
    FUi[2][0] = Ui[1]*Ui[2]/Ui[0];
    FUi[2][1] = (Ui[2]*Ui[2]/Ui[0]) + (_g*Ui[0]*Ui[0]/2.);

    //std::cout<<_t<<" "<<Ui[0]<<" "<< FUi[0][0]<<" "<<FUi[0][1]<<" "<<FUi[1][0]<<" "<<FUi[1][1]<<" "<<FUi[2][0]<<" "<<FUi[2][1]<<std::endl;

    return FUi;
}



vector<double> Schema::Flux(vector<double> Uk, vector<double>Ul, vector<double> ne)
{
    vector<vector<double>> FUk, FUl;
    vector<double> F(3);

    FUk = this->Fonction_F(Uk);
    FUl = this->Fonction_F(Ul);

    for(int i=0; i<3; i++)
    {
        F[i] = (FUk[i][0]*ne[0])+(FUk[i][1]*ne[1])+(FUl[i][0]*ne[0])+(FUl[i][1]*ne[1])-(_bmax*(Uk[i]-Ul[i]));
        F[i] = F[i]/2.;
    }

    return F;
}


void Schema::Update_U()
{
    vector<vector<double>> Unp1(_N);
    vector<double> S(3);
    for (int i=0; i<_N; i++)
    {
        Unp1[i].resize(3);
        S = this->S(x, y);
        for(int j=0; j<3; j++)
        {
            Unp1[i][j] = _U[i][j] + S[j];
        }
    }
    const Eigen::VectorXd &Aire_maille = this->_msh->Get_triangles_area();
	const Eigen::VectorXd &Longueur_arrete = this->_msh->Get_edges_length();

    for (int i=0; i < this->_msh->Get_edges().size(); i++) 
    {
        int maille1 = this->_msh->Get_edges()[i].Get_T1();
		int maille2 = this->_msh->Get_edges()[i].Get_T2();

        vector<double> ne(2);
        ne[0] = -this->_msh->Get_edges_normal()(i,0);
        ne[1] = -this->_msh->Get_edges_normal()(i,1);

        vector<double> F(3);
        if (Longueur_arrete[i] != 0)
        {
        if (maille2 != -1)
        {
            F = this->Flux(_U[maille1],_U[maille2],ne);
            for(int j=0; j<3;j++)
            {
                Unp1[maille1][j] = Unp1[maille1][j] + (_dt*Longueur_arrete[i]/Aire_maille[maille1])*F[j];
                Unp1[maille2][j] = Unp1[maille2][j] - (_dt*Longueur_arrete[i]/Aire_maille[maille2])*F[j];

                //std::cout << Longueur_arrete[i]<< " " <<Aire_maille[maille2] << " " << F[j] << std::endl;
            }
        }
        else
        {
            // if (this->_msh->Get_edges()[i].Get_reference() == 10 || this->_msh->Get_edges()[i].Get_reference() == 12)
            // {
            //     vector<vector<double>> FUk;
            //     FUk = this->Fonction_F(_U[maille1]);
            //     for (int j=0; j<3;j++)
            //     {
            //         F[j] = -((FUk[j][0]*ne[0])+(FUk[j][1]*ne[1]));
            //         Unp1[maille1][j] = Unp1[maille1][j] - (_dt*Longueur_arrete[i]/Aire_maille[maille1])*F[j];
            //     }
                
            // }

            // if (this->_msh->Get_edges()[i].Get_reference() == 11 || this->_msh->Get_edges()[i].Get_reference() == 13)
            // {
            //     for (int j=0; j<3;j++)
            //     {
            //         F[j] = 0;
            //         Unp1[maille1][j] = Unp1[maille1][j] - (_dt*Longueur_arrete[i]/Aire_maille[maille1])*F[j];
            //     }
            //     std::cout << this->_msh->Get_edges()[i].Get_reference() << " " << i << std::endl;
                
            // }

            vector<vector<double>> FUk;
            FUk = this->Fonction_F(_U[maille1]);
            for (int j=0; j<3;j++)
            {
                F[j] = -((FUk[j][0]*ne[0])+(FUk[j][1]*ne[1]));
                Unp1[maille1][j] = Unp1[maille1][j] - (_dt*Longueur_arrete[i]/Aire_maille[maille1])*F[j];
            }
        }
        }
        
    }
    for (int i=0; i<_N; i++)
    {
        for(int j=0; j<3; j++)
        {
            _U[i][j] = Unp1[i][j];
        }
    }
}


void Schema::Save_Solution_h()
{
    // char buffer[100];
    // std::snprintf(buffer, sizeof(buffer), "../solution/valid_2D/valid_2D_%i.dat", _k);

    // const Eigen::Matrix<double, Eigen::Dynamic, 2> &Centre_maille = this->_msh->Get_triangles_center();

    // // Écriture du fichier
    // std::ofstream file_out;
    // file_out.open(buffer); // Utilisation de buffer pour ouvrir le fichier

    // if (file_out.is_open()) 
    // {
    //     for (int i = 0; i < _N; i++) 
    //     {
    //         double Xi = Centre_maille(i,0);
	// 		double Yi = Centre_maille(i,1);
    //         file_out << Xi << " " << Yi << " " << _U[i][0] << " " << _U[i][1]/_U[i][0] << " " << _U[i][2]/_U[i][0] << " " << this->GetT() << std::endl;
    //     }
    //     file_out.close();
    //     //std::cout << "Les données ont été écrites dans le fichier : " << buffer << std::endl;
    // } else {
    //     std::cout << "Impossible d'ouvrir le fichier : " << buffer << std::endl;
    // }





	string name_file =  "../solution/valid_h/sol_h_" + std::to_string(_k) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS hauteur float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,_U[i][0]) << endl;
	}
	solution << endl;

	solution.close();

}

void Schema::Save_Solution_U()
{
    string name_file =  "../solution/valid_U/sol_U_" + std::to_string(_k) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS vitesse_U float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	//double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		//solution << max(eps,_U[i][1]/_U[i][0]) << endl;
        solution << _U[i][1]/_U[i][0] << endl;
	}
	solution << endl;

	solution.close();
}


void Schema::Save_Solution_V()
{
    string name_file =  "../solution/valid_V/sol_V_" + std::to_string(_k) + ".vtk";
    _k++;
	unsigned int nb_vert = this->_msh->Get_vertices().size();

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS vitesse_U float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	//double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		//solution << max(eps,_U[i][1]/_U[i][0]) << endl;
        solution << _U[i][2]/_U[i][0] << endl;
	}
	solution << endl;

	solution.close();
}


#define _MOD_SCHEMA_CPP
#endif