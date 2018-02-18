#pragma once

#include <vector>

#include <iostream>
#include <fstream>
#include "Vector.h"
#include "kdtree.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
using namespace std;


#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"

typedef Eigen::Triplet<double> Trip;
typedef Eigen::SparseMatrix<double> SpMat;

struct Triangle
{
	int idx; //tris[idx] is this triangle
	int v1i, v2i, v3i; //triangle formed by verts[v1i]-verts[v2i]-verts[v3i]
	float area;
	int e1, e2, e3; // edges
	Vector tNormal;
	Vector axisAngles;
	int normalId;
	float dihedral;
	std::vector<Triangle *> neighborhood;
	Triangle(int i, int a, int b, int c) : idx(i), v1i(a), v2i(b), v3i(c) {};
	
	
};

struct Vertex
{
	int idx=-1; //verts[idx]
	float* coords, //coords[0] ~ x coord, ..
		 * normal; //direction
	float diff;
	vector< int > triList;
	vector< int > edgeList;
	vector< int > vertList;
	float dihedral;
	Vertex(int i, float* c) : idx(i), coords(c) {};
	Vertex(float x, float y, float z){
		coords = new float[3];
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	}

	~Vertex(){
		delete coords;
	}

};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i;
	float length;
	vector<int> tris;

	Edge(int i, int a, int b, float l) : idx(i), v1i(a), v2i(b), length(l) {};
};


namespace graph{

	struct Edge;
	typedef struct Vertex{

		Vertex(int i) : id(i){}
		int id;
		vector<Edge*> ins;
		vector<Edge*> outs;
	}Vertex;

	typedef struct Edge {

		Vertex* from;
		Vertex* to;
		Vector* v;
	}Edge;

}

class Mesh
{
public:
	vector< Triangle* > tris;
	vector< Vertex* > verts;
	vector< Edge* > edges;
	std::vector<Trip> tripletList;
	SpMat L;
	SpMat W;
	Eigen::VectorXd vx;
	Eigen::VectorXd vy;
	Eigen::VectorXd vz;
	Eigen::VectorXd dx;
	Eigen::VectorXd dy;
	Eigen::VectorXd dz;

	Eigen::VectorXd LTLx;
	Eigen::VectorXd LTLy;
	Eigen::VectorXd LTLz;

	Eigen::VectorXd cx;
	Eigen::VectorXd cy;
	Eigen::VectorXd cz;


	Eigen::VectorXd bx;
	Eigen::VectorXd by;
	Eigen::VectorXd bz;

	~Mesh(){

		for (int i = 0; i < tris.size(); i++)
			delete tris[i];

		for (int i = 0; i < edges.size(); i++)
			delete edges[i];

		for (int i = 0; i < verts.size(); i++)
			delete verts[i];


	}

	void loadOff(char* fName);
	
	void loadxyz(char* fName);
	
	void createCube(float sl);
	void assignNormalsToTriangles();
	void findNeighborhoodTriangles();
	int findEdgeBetween(int v, int w);
	void assignEdgesToTriangles();
	void addExtraEdgesToTheMesh(float threshold);
	void deform2(Mesh* target, float weigth);
	void deform(Mesh* target, float weight);
	void deform(Mesh* target, float weigth,vector<int> corrVector);
	void write(char* filename, bool isMesh=true);
	void initialize();
	void deform(Mesh* target, float weigth,int iterations);
	void deform(Mesh* target, float weigth, int iterations, vector<int> corrVector);
	int addVertex(float* c);
	vtkPolyData* getVTKPolyData(bool isMesh=true);
	vtkPolyData* getVTKPolyDataSubsetByVertices(vector<int> vertices,int color=110);
	Mesh* getSubsetByVertices(vector<int> vertices);
	void calculateDihedrals();
	void shiftMesh(float x, float y){

		for (int i = 0; i < this->verts.size(); i++){
			this->verts[i]->coords[0] += x;
			this->verts[i]->coords[1] += y;
		}

	}

	Eigen::Vector3f center;
	float calculateScale();
	void scale(float s);
	int smoothNeighborsNormal(int tid, int numOfProcessed);
	int smoothAllNormals(int tid);
	void angleXY();
	void toOrigin();
	void calculateAreas();
	void addTriangle(int v1i, int v2i, int v3i);
	vector<graph::Vertex*> createVectorGraph();
	void removeEdge(int id);
	std::map<std::pair<int, int>, int> edgeMap;

	
private:
	
	
	void addEdge(int v1i, int v2i,int tid=-1);
	bool makeVertsNeighbors(int v, int w,int tid=-1);
	
};