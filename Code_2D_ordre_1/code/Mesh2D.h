#ifndef _MESH_2D_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Vertex
{
private:
   Eigen::Vector2d _v_coor;
   int _ref;
public:
   Vertex();
   Vertex(double x, double y, int ref);
   void Print() const;
   const Eigen::Vector2d Get_coor() const {return _v_coor;};
};

class Edge
{
private:
   Eigen::Vector2i _v_edge;
   int _ref;
   int _t1, _t2;
   std::string _BC;
public:
   Edge();
   Edge(int vertex1, int vertex2, int ref, std::string BC);
   void Print() const;
   void Add_triangle(int t)
   {
      if (_t1 == -1)
      _t1 = t;
      else
      _t2 = t;
   }
   const Eigen::Vector2i& Get_vertices() const { return _v_edge;}
   int Get_T1() const { return _t1; };
   int Get_T2() const { return _t2; };
   int Get_reference() const { return _ref;};
   std::string Get_BC() const { return _BC;};
};

class Triangle
{
private:
   Eigen::Vector3i _v_triangle;
   int _ref;
public:
   Triangle();
   Triangle(int vertex1, int vertex2, int vertex3, int ref);
   void Print() const;
   const Eigen::Vector3i& Get_vertices() const { return _v_triangle; }
   int Get_reference() const { return _ref;}; 
};

class Mesh2D
{
private:
   // liste de tous les sommets
   std::vector<Vertex> _vertices;
   // liste de tous les triangles
   std::vector<Triangle> _triangles;
   // centre de tous les triangles
   Eigen::Matrix<double, Eigen::Dynamic, 2> _tri_center;
   // aire de tous les triangles
   Eigen::VectorXd _tri_area;
   // dx de tous les triangles
   Eigen::VectorXd _tri_h;
   // liste de toutes les arêtes
   std::vector<Edge> _edges;
   // liste de toutes les normales unitaires !!!
   Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_normal;
   // liste de toutes les longueurs d'arêtes
   Eigen::VectorXd _edg_length;
   // centre des aretes
   Eigen::Matrix<double, Eigen::Dynamic, 2> _edg_center;
   // coordonnées des arrêtes
   std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2> > _edg_coord;
   // vecteur de référence des BC
   const std::vector<int> _BC_ref;
   // vecteur de type
   const std::vector<std::string> _BC_type;

public:
   Mesh2D(const std::vector<int> & BC_ref, const std::vector<std::string> & BC_type);
   void Read_mesh(std::string name_mesh);
   void Build_triangles_center_and_area();
   void Build_edges_normal_length_and_center();

   const std::vector<Vertex> & Get_vertices() const {return _vertices;};

   const std::vector<Triangle> & Get_triangles() const {return _triangles;};
   const Eigen::Matrix<double, Eigen::Dynamic, 2> & Get_triangles_center() const {return _tri_center;};
   const Eigen::VectorXd & Get_triangles_area() const  {return _tri_area;};
   const Eigen::VectorXd & Get_triangles_length() const  {return _tri_h;};

   const std::vector<Edge> & Get_edges() const {return _edges;};
   const Eigen::VectorXd & Get_edges_length() const {return _edg_length;};
   const Eigen::Matrix<double, Eigen::Dynamic, 2> & Get_edges_normal() const {return _edg_normal;};
   const Eigen::Matrix<double, Eigen::Dynamic, 2> & Get_edges_center() const {return _edg_center;};
   const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2> > & Get_edges_coord() const {return _edg_coord;};

protected:
   void Add_single_edge(const Edge& edge, int ne, std::vector<int>& head_minv,
      std::vector<int>& next_edge, int& nb_edges);
   };

   #define _MESH_2D_H
   #endif
