/*
 * Vector.h
 *
 *  Created on: Oct 12, 2013
 *      Author: tayfun
 */

#ifndef VECTOR_H_
#define VECTOR_H_

using namespace std;

#include <vector>
#include <cmath>

class Vector {
public:
	Vector(float x, float y, float z);
	Vector();
	virtual ~Vector();

	void crossProduct(Vector & v2,Vector & v3);
	float dotProduct(const Vector & v2) const;
	void Normalize();
	float Norm() const;
	Vector Subtract(Vector & v2);
	Vector Add(const Vector & v2) const;
	Vector Multiply(float n) const;
	float calculateAngleBetween(Vector & v2);
	void printVector() const;
	//attributes of Vector class
	float X;
	float Y;
	float Z;
	float operator [] (int index);

};

#endif /* VECTOR_H_ */
