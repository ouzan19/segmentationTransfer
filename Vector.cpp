/*
 * Vector.cpp
 *
 *  Created on: Oct 12, 2013
 *      Author: tayfun
 */

#include "Vector.h"
#include <iostream> // TO BE DELETED

using namespace std; //TO BE DELETED

Vector::Vector() {
	// TODO Auto-generated constructor stub
	X = 0;
	Y = 0;
	Z = 0;
}

Vector::Vector(float x,float y,float z) {
	// TODO Auto-generated constructor stub
	X = x;
	Y = y;
	Z = z;
}

Vector::~Vector() {
	// TODO Auto-generated destructor stub
}

void Vector::crossProduct(Vector & v2,Vector & v3)
{
	v3.X = this->Y*v2.Z - this->Z*v2.Y;
	v3.Y = this->Z*v2.X - this->X*v2.Z;
	v3.Z = this->X*v2.Y - this->Y*v2.X;
}


void Vector::Normalize()
{
	float norm;
	norm = sqrt(X*X + Y*Y + Z*Z);

	if(norm <= 0) return;

	X = X/norm;
	Y = Y/norm;
	Z = Z/norm;

}

float Vector::dotProduct(const Vector & v2) const
{
	return ( this->X * v2.X ) + ( this->Y * v2.Y ) + ( this->Z * v2.Z );
}

Vector Vector::Subtract(Vector & v2)
{
	Vector v3;
	v3.X = this->X - v2.X;
	v3.Y = this->Y - v2.Y;
	v3.Z = this->Z - v2.Z;
	return v3;
}

Vector Vector::Add(const Vector & v2) const
{
	Vector v3;
	v3.X = this->X + v2.X;
	v3.Y = this->Y + v2.Y;
	v3.Z = this->Z + v2.Z;
	return v3;

}

Vector Vector::Multiply(float n) const
{
	Vector v3;
	v3.X = this->X * n;
	v3.Y = this->Y * n;
	v3.Z = this->Z * n;
	return v3;
}

float Vector::operator [] (int index){

	if (index == 0)
		return X;
	if (index == 1)
		return Y;
	if (index == 2)
		return Z;

	return X;
}


void Vector::printVector() const //TO BE DELETED
{
	cout<<"X: "<<X<<" Y: "<<Y<<" Z: "<<Z<<endl;
}

float Vector::Norm() const
{
	float norm;
	norm = sqrt(X*X + Y*Y + Z*Z);
	return norm;
}

float Vector::calculateAngleBetween(Vector & v2){

	float cosx = (dotProduct(v2)) / (Norm() * v2.Norm());
	if (cosx >= 1)
		return 0;
	if (cosx <= -1)
		return 180;

	float angle = acos(cosx) * 180.0 / 3.14159265;
	return angle;
	


}