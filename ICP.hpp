#ifndef ICP_HPP
#define ICP_HPP


#include "Mesh.h"
#include "Vector.h"
#include "kdtree.h"
#include "Eigen/Dense" 
//#include "Eigen/Sparse"
#include <iostream>
#include <fstream>




class ICP {
public:
	Eigen::Matrix4f globalTrans;
	Eigen::Matrix4f globalTrans1;
	Eigen::Matrix4f globalTrans2;
	Eigen::Matrix4f globalTrans3;
	Eigen::Matrix4f globalTrans4;
	bool svdBased;
	vector<Eigen::Matrix4f> Transformations;
	void align(Mesh *r, Mesh *i);
private:

	Mesh* reference;
	Mesh* input;

};




#endif //ICP_HPP