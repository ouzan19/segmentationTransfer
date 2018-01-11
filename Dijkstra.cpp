#include "Dijkstra.h"


using namespace Dijsktra;
DijsktraSP::DijsktraSP(){}

DijsktraSP::~DijsktraSP(){}



void DijsktraSP::setMesh(Mesh* mesh){


	this->mesh = mesh;

}


std::vector<Dijsktra::Node*> DijsktraSP::run(int end){


	
	std::vector<Dijsktra::Node*> vertices;
	std::vector<Dijsktra::Edge*> edges;

	int n = mesh->verts.size();
	for(int j=0; j<n-1 ; j++)
		{
			Dijsktra::Node* v = new Dijsktra::Node(j, 0);
			vertices.push_back(v);
		}

	vertices.push_back(new Dijsktra::Node(n-1, 0)); 
	vertices[n-1]->state = Dijsktra::LABELED;

	for(int i=0;i<mesh->edges.size();i++){
	
		int tail = mesh->edges[i]->v1i;
		int head = mesh->edges[i]->v2i;
		float length = mesh->edges[i]->length;

		    Dijsktra::Edge* edge = new Dijsktra::Edge(vertices[tail], vertices[head], length);
			edge->head->addIncomingEdge(edge);
			edge->tail->addOutgoingEdge(edge);
			edges.push_back(edge);

			Dijsktra::Edge* edge2 = new Dijsktra::Edge(vertices[head], vertices[tail], length);
			edge2->head->addIncomingEdge(edge2);
			edge2->tail->addOutgoingEdge(edge2);
			edges.push_back(edge2);
	
	
	}

	
	Dijsktra::FibonacciHeap* heap = new Dijsktra::FibonacciHeap();
	
	heap->insertVertex(vertices[end]);
	
	bool abort = false;
	long j = 0;
	// Scan
	do
	{
		// Delete minimum path
		Dijsktra::Node* v = heap->deleteMin();
		
		v->state = Dijsktra::SCANNED;
		
		for(int i = 0; i < v->incomingEdges.size(); i++)
		{
			Dijsktra::Edge* currentEdge = v->incomingEdges[i];
			Dijsktra::Node* headOfCurrentEdge = currentEdge->tail;

			if(headOfCurrentEdge->state != Dijsktra::SCANNED)
				{
				if(headOfCurrentEdge->state == Dijsktra::UNLABELED)
				{
					// Insert a vertex with infinite key
					headOfCurrentEdge->state = Dijsktra::LABELED;
					headOfCurrentEdge->pred = v;
					headOfCurrentEdge->key = v->key + currentEdge->length;
					heap->insertVertex(headOfCurrentEdge);
				}
				else if(headOfCurrentEdge->key > v->key + currentEdge->length )
				{
					// decrease the key of a vertex with finite key
					headOfCurrentEdge->pred = v;
					heap->decreaseKey(v->key + currentEdge->length, headOfCurrentEdge);
				}
			}
		}
	}
	while(!heap->isEmpty());

	for(int j=0;j<edges.size();j++)
			delete edges[j];

	
	return vertices;
	
}



void DijsktraSP::writeMatrix(vector<vector<float>> matrix,  char* filename){

	ofstream f;
	f.open(filename);
	int size = mesh->verts.size();

	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++)
			f << matrix[i][j] << " ";
		f << std::endl;
	}

	f.close();
}


vector<vector<float>> DijsktraSP::createMatrix(){

	vector<vector<float>> matrix;
	int size = mesh->verts.size();

	for(int i=0;i<size;i++){
	
		vector<float> temp;
		std::vector<Dijsktra::Node*> vertices = run(i);

		for(int j=0;j<size;j++)
			temp.push_back(vertices[j]->key);

		matrix.push_back(temp);
		//for(int j=0;j<vertices.size();j++)
			//delete vertices[j];

	}

	return matrix;
}











