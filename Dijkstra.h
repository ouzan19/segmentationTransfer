#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

#include "Mesh.h"
#include "FibonacciHeap.h"
#include <fstream>


namespace Dijsktra{
class DijsktraSP{


public:
	DijsktraSP();
	~DijsktraSP();
	void setMesh(Mesh* mesh);
	std::vector<Dijsktra::Node*> run(int end);
	void writeMatrix(vector<vector<float>> matrix, char* filename);
	vector<vector<float>> createMatrix();
private:
	Mesh* mesh;





};



}

#endif // DIJKSTRA_HPP