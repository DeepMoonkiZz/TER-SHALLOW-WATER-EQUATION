#ifndef _MESH_2D_CPP

#include "Mesh2D.h"
#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

Vertex::Vertex()
{
   this->_v_coor[0] = -10000; this->_v_coor[1] = -10000; this->_ref = -1;
}

Vertex::Vertex(double x, double y, int ref) : _ref(ref)
{
   this->_v_coor[0] = x; this->_v_coor[1] = y;
}

void Vertex::Print() const
{
   cout << "[x, y] = [" << this->_v_coor[0] << " " << this->_v_coor[1] << "];" << endl;
   cout << "ref = " << this->_ref << endl;
}

Edge::Edge()
{
   this->_v_edge[0] = -1; this->_v_edge[1] = -1; this->_ref = -1;
}

Edge::Edge(int vertex1, int vertex2, int ref, std::string BC) : _ref(ref), _BC(BC)
{
   // sort
   if (vertex1 > vertex2)
   {
      this->_v_edge[0] = vertex2;
      this->_v_edge[1] = vertex1;
   }
   else
   {
      this->_v_edge[0] = vertex1;
      this->_v_edge[1] = vertex2;
   }
   this->_t1 = -1;
   this->_t2 = -1;
}


void Edge::Print() const
{
   cout << "[pt1, pt2] = [" << this->_v_edge[0] << " " << this->_v_edge[1] << "];" << endl;
   cout << "[t1, t2] = [" << this->_t1 << " " << this->_t2 << "];" << endl;
   cout << "ref = " << this->_ref << endl;
}

Triangle::Triangle()
{
   this->_v_triangle[0] = -1; this->_v_triangle[1] = -1; this->_v_triangle[2] = -1; this->_ref = -1;
}

Triangle::Triangle(int vertex1, int vertex2, int vertex3, int ref) : _ref(ref)
{
   this->_v_triangle[0] = vertex1; this->_v_triangle[1] = vertex2; this->_v_triangle[2] = vertex3;
}

void Triangle::Print() const
{
   cout << "[pt1, pt2, pt3] = [" << this->_v_triangle[0] << " " << this->_v_triangle[1] << " " << this->_v_triangle[2] << "];" << endl;
   cout << "ref = " << this->_ref << endl;
}

Mesh2D::Mesh2D( const std::vector<int> & BC_ref, const std::vector<std::string> & BC_type ) :
_BC_ref(BC_ref), _BC_type(BC_type)
{
}

void Mesh2D::Build_triangles_center_and_area()
{
   this->_tri_center.resize(this->_triangles.size(),2);
   this->_tri_area.resize(this->_triangles.size());
   this->_tri_h.resize(this->_triangles.size());

   for (unsigned int i = 0; i < this->_triangles.size(); i++)
   {
      int n1 = this->_triangles[i].Get_vertices()(0);
      int n2 = this->_triangles[i].Get_vertices()(1);
      int n3 = this->_triangles[i].Get_vertices()(2);

      double x1 = this->_vertices[n1].Get_coor()(0), y1 = this->_vertices[n1].Get_coor()(1);
      double x2 = this->_vertices[n2].Get_coor()(0), y2 = this->_vertices[n2].Get_coor()(1);
      double x3 = this->_vertices[n3].Get_coor()(0), y3 = this->_vertices[n3].Get_coor()(1);

      this->_tri_area(i) = 0.5*fabs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
      this->_tri_h(i) = (sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
      +sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))
      +sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)))/3.;
      // centre du triangle
      this->_tri_center(i,0) = (x1 + x2 + x3)/3.0;
      this->_tri_center(i,1) = (y1 + y2 + y3)/3.0;
   }
}

void Mesh2D::Build_edges_normal_length_and_center()
{
   this->_edg_coord.resize(2);
   this->_edg_coord[0].resize(this->_edges.size(),2);
   this->_edg_coord[1].resize(this->_edges.size(),2);

   this->_edg_center.resize(this->_edges.size(),2);
   this->_edg_normal.resize(this->_edges.size(),2);
   this->_edg_length.resize(this->_edges.size());

   Eigen::Vector2d diff;

   for (unsigned int i = 0; i < this->_edges.size(); i++)
   {
      int t1 = this->_edges[i].Get_T1();
      int n1 = this->_edges[i].Get_vertices()(0);
      int n2 = this->_edges[i].Get_vertices()(1);

      if (t1 >= 0)
      {
         double x1 = this->_vertices[n1].Get_coor()(0), y1 = this->_vertices[n1].Get_coor()(1);
         double x2 = this->_vertices[n2].Get_coor()(0), y2 = this->_vertices[n2].Get_coor()(1);

         this->_edg_coord[0](i,0) = x1; this->_edg_coord[0](i,1) = y1;
         this->_edg_coord[1](i,0) = x2; this->_edg_coord[1](i,1) = y2;
         this->_edg_length(i) = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

         // centre de l'arete
         this->_edg_center(i,0) = 0.5*(x1+x2);
         this->_edg_center(i,1) = 0.5*(y1+y2);

         // vecteur entre le centre de l'arete et le centre de l'element e1
         diff = this->_edg_center.row(i) - this->_tri_center.row(t1);

         // normale suivant un sens arbitraire
         this->_edg_normal(i,0) = y1 - y2;
         this->_edg_normal(i,1) = x2 - x1;

         double scal = diff(0)*this->_edg_normal(i,0) + diff(1)*this->_edg_normal(i,1);
         if (scal < 0)
         {
            // on change le signe de la normale
            // pour avoir une normale sortante de e1 vers e2
            this->_edg_normal.row(i) = -this->_edg_normal.row(i);
         }
         this->_edg_normal.row(i) /= this->_edg_length(i);
      }
   }
}

// methode interne qui rajoute une arete
void Mesh2D::Add_single_edge(const Edge& edge, int ne, vector<int>& head_minv, vector<int>& next_edge, int& nb_edges)
{
   int n1 = edge.Get_vertices()(0);
   int n2 = edge.Get_vertices()(1);
   int ref = edge.Get_reference();
   std::string BC_type = edge.Get_BC();

   bool exist = false;
   // we look at the list of edges leaving from n1
   // if we find the same edge than n1->n2 we add the edge
   for (int e = head_minv[n1]; e != -1; e = next_edge[e])
   {
      if (this->_edges[e].Get_vertices()(1) == n2)
      {
         if (ne >= 0)
         {
            this->_edges[e].Add_triangle(ne);
         }
         exist = true;
      }
   }

   // if the edge has not been found, we create it
   if (!exist)
   {
      // we initialize the edge
      this->_edges[nb_edges] = Edge(n1, n2, ref, BC_type);
      if (ne >= 0)
      {
         this->_edges[nb_edges].Add_triangle(ne);
      }
      // we update the arrays next_edge and head_minv
      next_edge[nb_edges] = head_minv[n1];
      head_minv[n1] = nb_edges;
      nb_edges++;
   }
}

void Mesh2D::Read_mesh(string name_mesh)
{
   ifstream mesh_file(name_mesh.data());
   if (!mesh_file.is_open())
   {
      cout << "Unable to open file " << name_mesh << endl;
      exit(0);
   }
   else
   {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading mesh: " << name_mesh << endl;
   }

   string file_line;
   vector<Edge> edges_boundary;
   int dim = 3;

   while (!mesh_file.eof())
   {
      getline(mesh_file, file_line);
      if (file_line.find("Dimension") != std::string::npos)
      {
         mesh_file >> dim;
      }
      else if (file_line.find("Vertices") != std::string::npos)
      {
         int nb_vertices(0);
         mesh_file >> nb_vertices;
         cout << "Number of vertices  (" << nb_vertices << ")" << endl;
         this->_vertices.resize(nb_vertices);
         for (int i = 0 ; i < nb_vertices ; ++i)
         {
            double x,y,z; int ref;
            mesh_file >> x >> y >> z >> ref;
            this->_vertices[i] = Vertex(x, y, ref);
         }
      }
      else if (file_line.find("Edges") != std::string::npos)
      {
         int nb_edges(0);
         mesh_file >> nb_edges;
         cout << "Number of edges (" << nb_edges << ")" << endl;
         edges_boundary.resize(nb_edges);
         int n1, n2, ref;
         for (int i = 0 ; i < nb_edges ; ++i)
         {
            mesh_file >> n1 >> n2 >> ref;
            n1--; n2--;
            std::string BC_type("none");
            for (unsigned int i=0 ; i < this->_BC_ref.size() ; i++)
            {
               if (ref == this->_BC_ref[i])
               {
                  BC_type = this->_BC_type[i];
               }
            }
            if (BC_type == "none")
            {
               cout << "Problem with BC in your mesh (reference or type are wrong)" << endl;
               exit(0);
            }
            edges_boundary[i] = Edge(n1, n2, ref, BC_type);
         }
      }
      else if (file_line.find("Triangles") != std::string::npos)
      {
         int nb_triangles(0);
         mesh_file >> nb_triangles;
         cout << "Number of triangles (" << nb_triangles << ")" << endl;
         this->_triangles.resize(nb_triangles);
         for (int i = 0 ; i < nb_triangles ; ++i)
         {
            int vertex1, vertex2, vertex3, ref;
            mesh_file >> vertex1 >> vertex2 >> vertex3 >> ref;
            vertex1--; vertex2--; vertex3--;
            this->_triangles[i] = Triangle(vertex1, vertex2, vertex3, ref);
         }
      }
   }

   cout << "---------Edges and Associated Triangles----------" << endl;
   // Toutes les aretes exterieures du maillage sont presentes
   int nb_edges = (3*_triangles.size() + edges_boundary.size())/2;
   this->_edges.resize(nb_edges);

   int nb_vertices = this->_vertices.size();
   vector<int> head_minv(nb_vertices, -1);
   vector<int> next_edge(nb_edges, -1);

   // on rajoute d'abord les aretes du bord
   nb_edges = 0;
   for (unsigned int i = 0; i < edges_boundary.size(); i++)
   {
      this->Add_single_edge(edges_boundary[i], -1, head_minv, next_edge, nb_edges);
   }

   // ensuite les aretes interieures
   for (unsigned int i = 0; i < this->_triangles.size(); i++)
   {
      const Eigen::Vector3i& nv = this->_triangles[i].Get_vertices();
      for (unsigned int j = 0; j < 3; j++)
      {
         Edge edge(nv(j), nv((j+1)%3), 0, "none");
         this->Add_single_edge(edge, i, head_minv, next_edge, nb_edges);
      }
   }

   cout << "-----------Triangles center and area-------------" << endl;
   Build_triangles_center_and_area();

   cout << "------------Edges Normal --------------" << endl;
   Build_edges_normal_length_and_center();

   cout << "-------------------------------------------------" << endl;
}

#define _MESH_2D_CPP
#endif
