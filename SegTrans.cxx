/*=========================================================================

  Program:   Visualization Toolkit
  Module:    Cube.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This example shows how to manually create vtkPolyData.

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderer.h"
#include <iostream>
#include <fstream>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkIdTypeArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkRendererCollection.h>
#include <vtkProperty.h>
#include <vtkPlanes.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPointSource.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkAreaPicker.h>
#include <vtkExtractGeometry.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include <vtkVersion.h>
#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include "Mesh.h"
#include "ICP.hpp"
#include "Dijkstra.h"
#include <vtkStructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkArrowSource.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <time.h>
#include <vtkTextProperty.h>
#include <set>

class HighlightInteractorStyle;
#define VTKISRBP_ORIENT 0
#define VTKISRBP_SELECT 1
using namespace std;
using namespace Eigen;
vtkPolyData *cube,*currPolyData;
vtkSmartPointer<vtkPolyDataMapper> mapper;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<HighlightInteractorStyle>  style;
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
bool* isSelected;
int c = 0;
int faceNo = 1;
float s;
vector<int> segments;
int nVerts;
Mesh* mesh, *mesh2, *averageFaceTemp, *averageFace;
int numOfSegments = 36;
float colorPerSegment = 360 / numOfSegments;
Mesh *initialFace, *currFace;
bool pointCloudMode = false;
int numOfSegs;
string name;

pair<float,float> HammingDistance(Mesh* mesh1, vector<int> segments1, vector<int> segments2);
float randIndex(vector<int> segments1, vector<int> segments2);
pair<float, float> consistencyError(Mesh* mesh, vector<int> segments1, vector<int> segments2);
Mesh* getOnlySegmented(Mesh* mesh, vector<int> segments, vector<int> &oldIndices);
float cutDiscrepancy(Mesh* mesh, vector<int> segments1, vector<int> segments2);

float dist2Between(float* p1, float* p2){

	float sum = 0;
	for (int i = 0; i < 3; i++)
		sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);

	return sum;
}

float PCAError(float** meshSubSet, float** voxSubSet, int n1, int n2, float** p1Shifted, float** p2Shifted,float** R1,float** R2){

	int* closestPntInP1 = new int[n2];

#ifdef KDTREE_FOR_CLOSEST_PNTS
	//k-d tree holding the mesh vertex coordinates and their corresponding idxs to speed up the closest pnt search below
	kdtree* kd = kd_create(3);
	kdres* resultSet;
	int* vertIdx;
	float vertCoord[3];
#endif

	//apply R1 & R2 rotation coming from alignPrincipalAxes() to each spectral point (icp input)
	for (int bv = 0; bv < n1; bv++)
	{
		//note that, meshSubSet changes in each iteration, but not to worry 'cos p1Shifted constant and hence new meshSubSet are all valid
		//aligning rotation applied to mean-shifted k=3-d coords (works)
		meshSubSet[bv][0] = R1[0][0] * p1Shifted[bv][0] + R1[0][1] * p1Shifted[bv][1] + R1[0][2] * p1Shifted[bv][2];
		meshSubSet[bv][1] = R1[1][0] * p1Shifted[bv][0] + R1[1][1] * p1Shifted[bv][1] + R1[1][2] * p1Shifted[bv][2];
		meshSubSet[bv][2] = R1[2][0] * p1Shifted[bv][0] + R1[2][1] * p1Shifted[bv][1] + R1[2][2] * p1Shifted[bv][2];
#ifdef KDTREE_FOR_CLOSEST_PNTS
		int* tmp = new int[1];
		*tmp = bv;
		kd_insert3f(kd, meshSubSet[bv][0], meshSubSet[bv][1], meshSubSet[bv][2], tmp);
#endif
	}
	for (int bv = 0; bv < n2; bv++)
	{
		voxSubSet[bv][0] = R2[0][0] * p2Shifted[bv][0] + R2[0][1] * p2Shifted[bv][1] + R2[0][2] * p2Shifted[bv][2];
		voxSubSet[bv][1] = R2[1][0] * p2Shifted[bv][0] + R2[1][1] * p2Shifted[bv][1] + R2[1][2] * p2Shifted[bv][2];
		voxSubSet[bv][2] = R2[2][0] * p2Shifted[bv][0] + R2[2][1] * p2Shifted[bv][1] + R2[2][2] * p2Shifted[bv][2];
	}

	//////////////// error of current alignment ////////////////
	//before error computation, i need closest-pnt correspondence b/w meshSubSet & voxSubSet
	float dist;
#ifdef KDTREE_FOR_CLOSEST_PNTS
	for (int i2 = 0; i2 < n2; i2++) //for each voxel
	{
		resultSet = kd_nearest3f(kd, voxSubSet[i2][0], voxSubSet[i2][1], voxSubSet[i2][2]); //find the mesh vertex closest to voxSubSet[i2]
		vertIdx = (int*)kd_res_itemf(resultSet, vertCoord); //get/retrieve the data and position of the current result item; not interested in the vertCoord though
		closestPntInP1[i2] = *vertIdx;
	}
#else
	float minDista = 9999999999999999999.0f, maxDista = -99999999999999999.0f; //for dist normalization
	for (int i2 = 0; i2 < n2; i2++)
	{
		float minDist = 9999999999999999999.0f; //INF's inside the loop 'cos i fill closestPntInP1[] per element
		for (int i1 = 0; i1 < n1; i1++)
		{
			dist = //distanceBetweenKD(voxSubSet[i2], meshSubSet[i1], k); is the correct way but to save time assume that k = 3 always and also discard sqrt by using dist2 function
				dist2Between(voxSubSet[i2], meshSubSet[i1]);
			if (dist < minDist)
			{
				closestPntInP1[i2] = i1;
				minDist = dist;
			}
			if (dist < minDista) //global min
				minDista = dist;
			if (dist > maxDista)
				maxDista = dist;
		}
	} //end of i2
#endif
	//compute the new alignment error as the sum of squared distances b/w closest pnts (or squared avgGd2Roots differences)
	float alignErr = 0.0f;
	for (int bv2 = 0; bv2 < n2; bv2++)
	{
		//distance b/w (target) mesh2.bv and its correspondence in mesh1 for this iteration (closestBaseInP1)
		dist = pow(voxSubSet[bv2][0] - meshSubSet[closestPntInP1[bv2]][0], 2.0f) +
			pow(voxSubSet[bv2][1] - meshSubSet[closestPntInP1[bv2]][1], 2.0f) +
			pow(voxSubSet[bv2][2] - meshSubSet[closestPntInP1[bv2]][2], 2.0f);
		alignErr += sqrt(dist); //do sqrt() above 'cos dist needs to be normalized just-below and min/maxDista is sqrt'ed values
	}
	
#ifdef KDTREE_FOR_CLOSEST_PNTS
	kd_res_free(resultSet);
	kd_free(kd);
#endif

	return alignErr / n2;
}


void PCAposeNormalization(float** meshSubSet, float** voxSubSet, int n1, int n2, float** bestMeshRot, float ** bestVoxRot)
{

	int k = 3;
	//employs PCA based pose normalization which transforms k-dimensional n1-size meshSubSet s.t. its principle axes are aligned with the conventional x- y- z- axes; same done for n2-size voxSubSet hence they're aligned

	std::cout << "principle axes being aligned with the conventional x- y- z- axes..\n";

	//assumes meshSubSet & voxSubSet are at the ~same scale 'cos for each rotation candidate below i use alignErr based on closest point distances sum so they better be matched at the correct scale

	//////////////// rotation ///////////////////
	//meshSubSet already translated so that center at origin (0,0,0) so mean-shifted coords is meshSubSet itself but i want to put it in kxn1 mtrx ms so that cov is computed by 1 mtrx multiplication as ms x ms.transpose()
	//principal axes for meshSubSet
	MatrixXf ms(k, n1); //mean-shifted kxn1 matrix for the mesh1; X is for dynamic size f is for float
	for (int i = 0; i < k; i++) //for all k rows
		for (int j = 0; j < n1; j++) //fill n1 columns
			ms(i, j) = meshSubSet[j][i];

	MatrixXf mst = ms.transpose(); //transpose of ms
	
	Matrix3f cov1 = ms * mst; //covariance mtrx of mesh1 (info: since resulting cov is always real symmetric, its eigvals are nonnegative real numbers)

	//eigenvalues(): The eigenvalues are repeated according to their algebraic multiplicity, so there are as many eigenvalues as rows in the matrix. The eigenvalues are sorted in increasing order.
	//eigenvectors(): Returns a matrix whose columns are eigenvectors of the input mtrx corresponding to the eigenvalues, i.e. i'th col is the i'th value in the resulting eigenvalues()
	SelfAdjointEigenSolver<Matrix3f> eigensolver(cov1); //solve eigenvalue/vector for symmetric n by n matrices, a.k.a. selfadjoint matrices
	if (eigensolver.info() != Success) abort();
	Vector3f eig_vals = eigensolver.eigenvalues(); //given in ascending order
	Matrix3f principal_axes = eigensolver.eigenvectors();
	//cout << "The eigenvalues of cov are:\n" << eig_vals << "\nHere's a matrix whose columns are eigenvectors of cov corresponding to these eigenvalues:\n" << principal_axes << endl;

	//eig_vals in ascending order as we go from (0) to (n); last column of principal_axes = eigvector1 (corresponding to largest eigval), n-1'th column = eigvec2, and so on..
	float** principalAxes1 = new float*[k]; //eigvecs returned are already normalized; first row of principalAxes1 is the eigvector w/ the largest eigenvalue
	Matrix3f mat; //just for the inversion test of determinant == -1

	for (int j = 0; j < k; j++)
	{
		principalAxes1[j] = new float[k];
		for (int c = 0; c < k; c++)
			mat(j, c) = principalAxes1[j][c] = principal_axes(c, k - j - 1); //so the last column principal_axes(*, 2) is put in the first row principalAxes1[0][*], then principal_axes(*, 1) put in principalAxes1[1][*], and so on			
	}

	float det = mat.determinant(), tiny = 0.000001f; //thresholded det==-1 test
	if (det > -1.0f - tiny && det < -1.0f + tiny){ //improper rotation that involves reflection which causes inversion, i.e. cw faces flip to ccw faces; so negate 1 axis and get rid of it
		for (int c = 0; c < k; c++)
			principalAxes1[0][c] = -principalAxes1[0][c];



		std::cout << "improper determinant" << std::endl;
	}

	//aligning rotation applied to mean-shifted coords (ms = transpose of meshSubSet[i] - m1[i], or equivalently ms = kxn1)
	float** p1Shifted = new float*[n1];
	for (int i = 0; i < n1; i++)
	{
		p1Shifted[i] = new float[k];
		for (int c = 0; c < k; c++)
			p1Shifted[i][c] = ms(c, i);
	}

	//principal axes for voxSubSet
	MatrixXf ms2(k, n2);
	for (int i = 0; i < k; i++) //for all k rows
		for (int j = 0; j < n2; j++) //fill n2 columns
			ms2(i, j) = voxSubSet[j][i];
	MatrixXf mst2 = ms2.transpose(); //transpose of ms2
	Matrix3f cov2 = ms2 * mst2;
	SelfAdjointEigenSolver<Matrix3f> eigensolver2(cov2); //solve eigenvalue/vector for symmetric n by n matrices, a.k.a. selfadjoint matrices
	if (eigensolver2.info() != Success) abort();
	eig_vals = eigensolver2.eigenvalues();
	principal_axes = eigensolver2.eigenvectors();
	float** principalAxes2 = new float*[k];
	for (int j = 0; j < k; j++)
	{
		principalAxes2[j] = new float[k];
		for (int c = 0; c < k; c++)
			principalAxes2[j][c] = principal_axes(c, k - j - 1);
	}
	float** p2Shifted = new float*[n2];
	for (int i = 0; i < n2; i++)
	{
		p2Shifted[i] = new float[k];
		for (int c = 0; c < k; c++)
			p2Shifted[i][c] = ms2(c, i);
	}

	//align each principal axes to traditional xyz-axes via R1 & R2
	float** R1 = new float*[k], ** R2 = new float*[k]; //kxk rotation mtrx to be applied to meshSubSet & voxSubSet, respectively (won't be a direct mtrx multiplication though)
	for (int i = 0; i < k; i++)
	{
		R1[i] = new float[k];		R2[i] = new float[k];
	}


	
	int* closestPntInP1 = new int[n2]; //filled for voxSubSet entries; closestPntInP1[3] = 7 means voxSubSet[3] and meshSubSet[7] are closest		
	float minAlignErr = 999999999999999999.0f, //axes leading to minimum alignErr are the desired ones
		alignErr = 100.0f;
	/*	int* perms = new int[k], fact3 = 6; //3! = 6 and perms[] at r=0 call = 6 5 4 3 2 1, at r=717 call = 1 6 5 2 4 3, and so on, for k = 6 (6!.txt) case
	FILE* fPtr;
	if (! (fPtr = fopen("permutations\\0-start\\3!.txt", "r")))
	{
	cout << "permutation file not found\n";
	exit(0);
	}
	for (int r = 0; r < fact3; r++) //k! reordering (just cancel lines from here to principalAxes2[j][c] = .. (inclusive) to skip this reordering stuff)*/
	{
		/*		//refill principalAxes2 according to the current reordering (one of k! permutations)
		//one permutation of 3 elements (1..3) is 1 line of 3!.txt; get this from file
		for (int p = 0; p < k; p++)
		fscanf(fPtr, "%d ", &perms[p]);
		for (int j = 0; j < k; j++)
		for (int c = 0; c < k; c++)
		principalAxes2[j][c] = principal_axes(c, k-perms[j]-1); //*/


		float** principalAxes22 = new float*[k];
		for (int j = 0; j < k; j++)
		{
			principalAxes22[j] = new float[k];
			for (int c = 0; c < k; c++)
				principalAxes22[j][c] = principalAxes2[j][c];
		}


		float k1=1, k2=1, k3=1;
		//sign flip (it is enough to keep pax1 fixed, and try 8 pax2 possibilities; so, canceling this pax1 outer loop (similar to alignWith()))
		int jj = -1, maxJ = (int)pow(2.0f, k); //= 16 for k = 4, 32 for k = 5
		for (int j = 0; j < maxJ ; j++) //try each orientation of pax2
		{
			if (j == 0){

				k1 = 1;
				k2 = 1;
				k3 = 1;
			}

			else if (j == 1){

				k1 = 1;
				k2 = 1;
				k3 = -1;
			}

			else if (j == 2){

				k1 = 1;
				k2 = -1;
				k3 = 1;

			}

			else if (j == 3){

				k1 = 1;
				k2 = -1;
				k3 = -1;
			}

			else if (j == 4){

				k1 = -1;
				k2 = 1;
				k3 = 1;
			}

			else if (j == 5){

				k1 = -1;
				k2 = 1;
				k3 = -1;
			}

			else if (j == 6){

				k1 = -1;
				k2 = -1;
				k3 = 1;

			}

			else if (j == 7){

				k1 = -1;
				k2 = -1;
				k3 = -1;

			}


			for (int c = 0; c < 3; c++){
				principalAxes2[0][c] = k1*principalAxes22[0][c];
				principalAxes2[1][c] = k2*principalAxes22[1][c];
				principalAxes2[2][c] = k3*principalAxes22[2][c];
			}
			///////// +- eigvector handling for k=3 case ends ////////////

			//try rotation w/ current/new principalAxes
			for (int a = 0; a < k; a++)
				for (int b = 0; b < k; b++)
				{
					//transform unit vectors principalAxes[] onto traditional xyz-axes
					R1[a][b] = principalAxes1[a][b];
					//transform unit vectors principalAxes2[] onto traditional xyz-axes
					R2[a][b] = principalAxes2[a][b];
					mat(a, b) = R2[a][b];
				}
			//cout << "det of " << j << "'th voxel rotation: " << mat.determinant() << "\t\t" << mat.inverse().determinant() << endl; //if det=1 then inv.det=1; if =-1 then inv.det=-1
			float det = mat.determinant(), tiny = 0.000001f; //thresholded det==-1 test
			if (det > -1.0f - tiny && det < -1.0f + tiny) //improper rotation that involves reflection which causes inversion, i.e. cw faces flip to ccw faces; no need to negate 'cos negate version will come or already came in this loop; so just skip this reflected/bad one
			{
				//continue;

				for (int c = 0; c < k; c++)
					principalAxes2[0][c] = -principalAxes2[0][c];

				for (int a = 0; a < k; a++)
					for (int b = 0; b < k; b++)
						R2[a][b] = principalAxes2[a][b];


			}

			if (j == -1){


				R1[0][0] = 1;	R1[0][1] = 0;	R1[0][2] = 0;
				R1[1][0] = 0;	R1[1][1] = 1;	R1[1][2] = 0;
				R1[2][0] = 0;	R1[2][1] = 0;	R1[2][2] = 1;

				R2[0][0] = 1;	R2[0][1] = 0;	R2[0][2] = 0;
				R2[1][0] = 0;	R2[1][1] = 1;	R2[1][2] = 0;
				R2[2][0] = 0;	R2[2][1] = 0;	R2[2][2] = 1;
					
			}

			alignErr = PCAError(meshSubSet, voxSubSet, n1, n2, p1Shifted, p2Shifted, R1, R2);
			//alignErr = PCAError(voxSubSet, meshSubSet,n2, n1, p2Shifted, p1Shifted, R2, R1);

			//cout << j << "'th axes' closest-distances-sum = " << alignErr << "\n----------------------------------\n";
			if (alignErr < minAlignErr) //not -1 means a preknown patient loaded, so do a shortcut PCA using the preknown cool j'th axis
			{
				minAlignErr = alignErr;
				jj = j;
				//bestMeshRot = R1;	bestVoxRot = R2; this fails, i.e. bestMeshRot changes whenever R1 changes in a nonrelated iteration
				for (int a = 0; a < k; a++)
				for (int b = 0; b < k; b++)
				{
					bestMeshRot[a][b] = R1[a][b];
					bestVoxRot[a][b] = R2[a][b];
				}
			}


		} //end of j
		//cout << jj << "'th axes with a closest-distances-sum of " << minAlignErr << " selected\n";
	} //end of r

	//ofstream f1,f2;

	//f1.open("rot1.txt");
	//f2.open("rot2.txt");

	Matrix3f rot1, rot2;

	for (int a = 0; a < k; a++)
		for (int b = 0; b < k; b++)
		{
			rot1(a, b) = bestMeshRot[a][b];
			rot2(a, b) = bestVoxRot[a][b];
		}

	//f1 << rot1;
	//f2 << rot2;

	//f1.close();
	//f2.close();

	rot1 = rot1.inverse() * rot2;

	for (int a = 0; a < k; a++)
		for (int b = 0; b < k; b++)
		{

			bestVoxRot[a][b] = rot1(a, b);
		}

	//memo re-capture
	delete[] p1Shifted;	delete[] p2Shifted;
	delete[] principalAxes1;	delete[] principalAxes2;
	delete[] closestPntInP1;
	for (int i = 0; i < k; i++) delete[] R1[i]; delete[] R1;
	for (int i = 0; i < k; i++) delete[] R2[i]; delete[] R2;

	cout << "pca done!\n\n";
}

void PCAAlignment(Mesh* face1, Mesh* face2){

	Vector com1(0, 0, 0);
	Vector com2(0, 0, 0);
	for (int i = 0; i < face1->verts.size(); i++){


		com1.X += face1->verts[i]->coords[0];
		com1.Y += face1->verts[i]->coords[1];
		com1.Z += face1->verts[i]->coords[2];
	}

	for (int i = 0; i < face2->verts.size(); i++){


		com2.X += face2->verts[i]->coords[0];
		com2.Y += face2->verts[i]->coords[1];
		com2.Z += face2->verts[i]->coords[2];
	}

	com1 = com1.Multiply(1.0f / face1->verts.size());
	com2 = com2.Multiply(1.0f / face2->verts.size());

	//cout << "com1:	" << face1->verts[0]->coords[0] << "	" << face1->verts[0]->coords[1] << "	" << face1->verts[0]->coords[2] << endl;
	//cout << "com2:	" << face2->verts[0]->coords[0] << "	" << face2->verts[0]->coords[1] << "	" << face2->verts[0]->coords[2] << endl;

	float **coor1 = new float*[face1->verts.size()];
	for (int i = 0; i < face1->verts.size(); i++)
		coor1[i] = new float[3];


	float **coor2 = new float*[face2->verts.size()];
	for (int i = 0; i < face2->verts.size(); i++)
		coor2[i] = new float[3];


	for (int i = 0; i < face1->verts.size(); i++){

		for (int j = 0; j < 3; j++)
			coor1[i][j] = face1->verts[i]->coords[j] - com1[j];
	}

	for (int i = 0; i < face2->verts.size(); i++){

		for (int j = 0; j < 3; j++)
			coor2[i][j] = face2->verts[i]->coords[j] - com2[j];
	}

	float** rot1 = new float*[3], ** rot2 = new float*[3];
	for (int i = 0; i < 3; i++)
	{
		rot1[i] = new float[3];		rot2[i] = new float[3];
	}


	PCAposeNormalization(coor1, coor2, face1->verts.size(), face2->verts.size(), rot1, rot2);

	/*for (int i = 0; i < face1->verts.size(); i++){

	float x = face1->verts[i]->coords[0];
	float y = face1->verts[i]->coords[1];
	float z = face1->verts[i]->coords[2];

	face1->verts[i]->coords[0] = rot1[0][0] * x + rot1[0][1] * y + rot1[0][2] * z;
	face1->verts[i]->coords[1] = rot1[1][0] * x + rot1[1][1] * y + rot1[1][2] * z;
	face1->verts[i]->coords[2] = rot1[2][0] * x + rot1[2][1] * y + rot1[2][2] * z;


	}*/

	float x = com2[0];
	float y = com2[1];
	float z = com2[2];

	float xt = rot2[0][0] * x + rot2[0][1] * y + rot2[0][2] * z;
	float yt = rot2[1][0] * x + rot2[1][1] * y + rot2[1][2] * z;
	float zt = rot2[2][0] * x + rot2[2][1] * y + rot2[2][2] * z;

	xt = com1[0] - xt;
	yt = com1[1] - yt;
	zt = com1[2] - zt;


	//com2 = Vector(0, 0, 0);

	for (int i = 0; i < face2->verts.size(); i++){

		float x = face2->verts[i]->coords[0];
		float y = face2->verts[i]->coords[1];
		float z = face2->verts[i]->coords[2];

		face2->verts[i]->coords[0] = rot2[0][0] * x + rot2[0][1] * y + rot2[0][2] * z + xt;
		face2->verts[i]->coords[1] = rot2[1][0] * x + rot2[1][1] * y + rot2[1][2] * z + yt;
		face2->verts[i]->coords[2] = rot2[2][0] * x + rot2[2][1] * y + rot2[2][2] * z + zt;


	}





	/*float xt = com1[0] - com2[0];
	float yt = com1[1] - com2[1];
	float zt = com1[2] - com2[2];

	for (int i = 0; i < face2->verts.size(); i++){

	face2->verts[i]->coords[0] += xt;
	face2->verts[i]->coords[1] += yt;
	face2->verts[i]->coords[2] += zt;


	}*/

	for (int i = 0; i < 3; i++) delete[] rot1[i]; delete[] rot1;
	for (int i = 0; i < 3; i++) delete[] rot2[i]; delete[] rot2;

	for (int i = 0; i < face1->verts.size(); i++) delete[] coor1[i]; delete[] coor1;
	for (int i = 0; i < face2->verts.size(); i++) delete[] coor2[i]; delete[] coor2;

}


void initialFineRigidAlignment(Mesh* face1, Mesh* face2){

	ICP icp;
	icp.align(face1, face2);

}

void alignDatabase(){

	int lambda = 1;
	averageFace = new Mesh();

	if(!pointCloudMode)
		averageFace->loadOff("faces\\face1.off");
	else 
		averageFace->loadxyz("faces\\face1.xyz");

	s = averageFace->calculateScale();


	Mesh *initialFace = new Mesh();

	if (!pointCloudMode)
		initialFace->loadOff("faces\\face1.off");
	else
		initialFace->loadxyz("faces\\face1.xyz");
	

	int i = 2;
	string s1 = "faces\\face";
	string s2 = "faces\\corr";
	string s3 = "faces\\deformed";
	string s4 = "faces\\coarse";
	string s5 = "faces\\fine";
	
	float tt = 0;

	for (int face = faceNo; face <= faceNo; face++){
		
		averageFaceTemp = new Mesh();
		if (!pointCloudMode)
			averageFaceTemp->loadOff("faces\\face1.off");
		else
			averageFaceTemp->loadxyz("faces\\face1.xyz");

		averageFaceTemp->initialize();

		string path, path2, path3, path4, path5;
		if (!pointCloudMode){
			 path = s1 + to_string(face) + ".off";
			 path2 = s2 + to_string(face) + ".txt";
			 path3 = s3 + to_string(face) + ".off";
			 path4 = s4 + to_string(face) + ".off";
			 path5 = s5 + to_string(face) + ".off";
		}
		else{
			path = s1 + to_string(face) + ".xyz";
			path2 = s2 + to_string(face) + ".txt";
			path3 = s3 + to_string(face) + ".xyz";
			path4 = s4 + to_string(face) + ".xyz";
			path5 = s5 + to_string(face) + ".xyz";

		}
		mesh = new Mesh();
		
		if (!pointCloudMode)
			mesh->loadOff((char*)path.c_str());
		else
			mesh->loadxyz((char*)path.c_str());
		

		/*
		string path6 = "faces\\face" + to_string(face) + ".xyz";
		ofstream fp;
		fp.open(path6);
		fp << mesh->verts.size() << endl;
		for (int i = 0; i < mesh->verts.size(); i++){

			float* c = mesh->verts[i]->coords;
			fp << c[0] << "	" << c[1] << "	" << c[2] << endl;

		}
		fp.close();
		*/
		//continue;


		//float s1 = mesh->calculateScale();
		//mesh->scale(s/s1);

		const clock_t begin_time = clock();

		PCAAlignment(averageFaceTemp, mesh);
		//mesh->write((char*)path4.c_str(), !pointCloudMode);

		

		initialFineRigidAlignment(averageFaceTemp, mesh);
	
		//mesh->write((char*)path5.c_str(), !pointCloudMode);

		mesh->initialize();
		mesh->deform(averageFaceTemp, 100000, 10);
		//mesh->write((char*)path3.c_str(), !pointCloudMode);
		

		ofstream ff;
		ff.open(path2);
		kdtree* kd = kd_create(3);
		kdres* resultSet = NULL;
		int* voxelIdx;
		float voxelCoord[3];

		for (int i = 0; i < mesh->verts.size(); i++)
		{

			int* tmp = new int[1];
			*tmp = i;
			kd_insert3f(kd, mesh->verts[i]->coords[0], mesh->verts[i]->coords[1], mesh->verts[i]->coords[2], tmp);
			//delete tmp;
		}


		for (int i = 0; i < averageFaceTemp->verts.size(); i++)
		{

			resultSet = kd_nearest3f(kd, averageFaceTemp->verts[i]->coords[0], averageFaceTemp->verts[i]->coords[1], averageFaceTemp->verts[i]->coords[2]);
			voxelIdx = (int*)kd_res_itemf(resultSet, voxelCoord);

			ff << *voxelIdx << endl;

			//averageFace->verts[i]->coords[0] = (voxelCoord[0] + lambda*averageFace->verts[i]->coords[0]) / (lambda + 1);
			//averageFace->verts[i]->coords[1] = (voxelCoord[1] + lambda*averageFace->verts[i]->coords[1]) / (lambda + 1);
			//averageFace->verts[i]->coords[2] = (voxelCoord[2] + lambda*averageFace->verts[i]->coords[2]) / (lambda + 1);

		    //averageFaceTemp->verts[i]->coords[0] = averageFace->verts[i]->coords[0];
			//averageFaceTemp->verts[i]->coords[1] = averageFace->verts[i]->coords[1];
			//averageFaceTemp->verts[i]->coords[2] = averageFace->verts[i]->coords[2];
			kd_res_free(resultSet);

		}
		

		tt += float(clock() - begin_time) / CLOCKS_PER_SEC;
		lambda++;
		kd_clear(kd);
		kd_free(kd);
		ff.close();
		delete averageFaceTemp;
	}


	std::cout << tt/19 <<std::endl ;
	

	//averageFace->write("averageFaceInc.off");
}

bool isFirst = true;
// Define interaction style
class HighlightInteractorStyle : public vtkInteractorStyleRubberBandPick
{
public:
	
	static HighlightInteractorStyle* New();
	vtkTypeMacro(HighlightInteractorStyle, vtkInteractorStyleRubberBandPick);

	 HighlightInteractorStyle() : vtkInteractorStyleRubberBandPick()
	{
		
		this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		this->SelectedActor = vtkSmartPointer<vtkActor>::New();
		this->SelectedActor->SetMapper(SelectedMapper);
		extractPolyDataGeometry =
			vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
	}

	virtual void OnLeftButtonUp()
	{
		// Forward events
		vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

		if (this->CurrentMode == VTKISRBP_SELECT)
		{

			vtkSmartPointer<vtkSelectVisiblePoints> VisibleFilter = vtkSmartPointer<vtkSelectVisiblePoints>::New();
			VisibleFilter->SetInputData(this->PolyData);
			VisibleFilter->SetRenderer(renderer);
			VisibleFilter->Update();
			
			vtkPolyData *pp = VisibleFilter->GetOutput();

			std::cout << "Extracted " << pp->GetNumberOfPoints() << " cells." << std::endl;

			vtkPlanes* frustum = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetFrustum();

			extractPolyDataGeometry->SetInputData(this->PolyData);

			extractPolyDataGeometry->SetImplicitFunction(frustum);
			extractPolyDataGeometry->Update();

			int num = extractPolyDataGeometry->GetOutput()->GetNumberOfCells();
			std::cout << "Extracted " << num << " cells." << std::endl;

			vtkPolyData* selected = extractPolyDataGeometry->GetOutput();

			this->SelectedMapper->SetInputData(selected);

			this->SelectedMapper->ScalarVisibilityOff();

			num = selected->GetPoints()->GetNumberOfPoints();
			
			vtkDataArray * s = currPolyData->GetPointData()->GetScalars();
			int count = 0;
			vtkPoints* points = currPolyData->GetPoints();
			for (int i = 0; i < num; i++){

				double coords[3];
				selected->GetPoint(i, coords);

				for (int j = 0; j < pp->GetNumberOfPoints(); j++){

					double coords2[3];
					pp->GetPoint(j, coords2);

					if (abs(coords2[0] - coords[0]) < 0.001 && abs(coords2[1] - coords[1]) < 0.001 && abs(coords2[2] - coords[2]) < 0.001){


						for (int k = 0; k < points->GetNumberOfPoints(); k++){

							double coords3[3];
							points->GetPoint(k, coords3);

							if (abs(coords2[0] - coords3[0]) <0.001 && abs(coords2[1] - coords3[1]) < 0.001 && abs(coords2[2] - coords3[2]) < 0.001){

								if (!isSelected[k]){
									s->SetTuple1(k, c);
									segments[k] =  c / colorPerSegment;
									isSelected[k] =  true;
									count++;
								}
							}

						}
					}
				}
		
			}
			//std::cout <<"Updated:	" <<count << std::endl;
		
			currPolyData->GetPointData()->SetScalars(s);
			currPolyData->GetPointData()->GetScalars()->Modified();

			currPolyData->GetPointData()->SetScalars(s);
			currPolyData->GetPointData()->GetScalars()->Modified();
			
			
	
			//this->SelectedActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
			//this->SelectedActor->GetProperty()->SetPointSize(5);

			//this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(SelectedActor);
			this->GetInteractor()->GetRenderWindow()->Render();
			this->HighlightProp(NULL);
		}
	}

	virtual void OnRightButtonUp(){

		vtkInteractorStyleRubberBandPick::OnRightButtonUp();
		

	}

	virtual void OnKeyPress()
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();

		// Output the key that was pressed
		std::cout << "Pressed " << key << std::endl;

		// Handle an arrow key
		if (key == "space")
		{
			c += colorPerSegment;
			std::cout << "current color  "<<c << std::endl;
		}	

		if (key == "s"){
			/*
			int prevCount =0 ;
			while (true){
				int count = 0;
				for (int i = 0; i < currFace->verts.size(); i++)
					if (!isSelected[i])
						count++;

				std::cout << "Not matched:	" << count << std::endl;

				if (prevCount == count)
					break;

				prevCount = count;


				if (count == 0)
					break;

				for (int i = 0; i < currFace->verts.size(); i++){
					if (!isSelected[i]){


						int size = currFace->verts[i]->vertList.size();
						
						vector<int> candidates;
						for (int j = 0; j < size; j++){

							int n = currFace->verts[i]->vertList[j];

							int s = segments[n];
							if (s != -1){

								candidates.push_back(s);

							}
						}

						

						if (candidates.size()){
							std::map<int, int> m;
							for (auto c : candidates)
							{
								auto it = m.find(c);
								if (it != m.end())
									m[c]++;
								else
									m[c] = 1;
							}

							char mostFreq;
							int count2 = 0;
							for (auto mi : m)
								if (mi.second >= count2)
								{
									mostFreq = mi.first;
									count2 = mi.second;
								}

							isSelected[i] = true;
							segments[i] = mostFreq;
						}

						else {
							segments[i] = -1;
							isSelected[i] = false;

						}



					}
				}

				

			
			}

			
			*/
			ofstream f;
			f.open("segmentation_"+name+"_"+to_string(faceNo)+".txt");

			f << c / colorPerSegment << std::endl;

			for (int i = 0; i < segments.size(); i++)
				f << segments[i] << std::endl;

			f.close();

			

			/*
			vtkDataArray * s = cube->GetPointData()->GetScalars();
			for (int k = 0; k < currFace->verts.size(); k++){
				if (segments[k] != -1)
				s->SetTuple1(k, segments[k]*colorPerSegment );
			}

			cube->GetPointData()->SetScalars(s);
			cube->GetPointData()->GetScalars()->Modified();
			this->GetInteractor()->GetRenderWindow()->Render();
			*/


		}

		if (key == "m"){
			/*
			vector<int> oldVertices;
			Mesh* m = getOnlySegmented(currFace, segments, oldVertices);
			renderer->RemoveActor(renderer->GetActors()->GetLastActor());


			vector<int> bound = cutDiscrepancy(currFace, segments, segments);

			vtkPolyData *face;
			//cube->Delete();
			face = m->getVTKPolyData(!pointCloudMode);
			vtkDataArray * s = face->GetPointData()->GetScalars();


			for (int i = 0; i < bound.size(); i++){

				

					int v = find(oldVertices.begin(), oldVertices.end(), bound[i]) - oldVertices.begin();

					s->SetTuple1(v, 50 );

				
			}

			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(face);
			mapper->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);


			renderer->AddActor(actor);
			renderer->Modified();

			this->GetInteractor()->GetRenderWindow()->Render();
			*/
		}

		if (key == "l"){

			
			vtkPoints *points = currPolyData->GetPoints();
			segments.clear();
			ifstream f;
			f.open("segmentation_" + name + "_" + to_string(1) + ".txt");
			int segSize;
			f >> segSize;
			
			for (int i = 0; i < points->GetNumberOfPoints(); i++){

				int temp;
				
				f >> temp;
				
				segments.push_back(temp);

			}


			vtkDataArray * s = currPolyData->GetPointData()->GetScalars();
		
			for (int k = 0; k < points->GetNumberOfPoints(); k++){

					if (segments[k] != -1){
						s->SetTuple1(k, colorPerSegment*segments[k]);
						
					}
					else 
						s->SetTuple1(k, 360);
			}


			currPolyData->GetPointData()->SetScalars(s);
			currPolyData->GetPointData()->GetScalars()->Modified();

			

			//this->SelectedActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
			//this->SelectedActor->GetProperty()->SetPointSize(5);

			//this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(SelectedActor);
			this->GetInteractor()->GetRenderWindow()->Render();

		}

		if (key == "o"){

			int c = 0;
			for (int i = 0; i < segments.size();i++)
			if (segments[i] != -1)
				c++;

			std::cout << "segments	" << c << "		" << segments.size() << std::endl;
		}


		if (key == "n"){

			std::cout << faceNo + 1 << "  loading..." << std::endl;
			if (faceNo > 1)
				renderer->RemoveActor(renderer->GetActors()->GetLastActor());

			

			string path;
			if (!pointCloudMode)
				path = "faces\\face" + std::to_string(faceNo+1) + ".off";
			else
				path = "faces\\face" + std::to_string(faceNo+1) + ".xyz";

			Mesh *currentFace = new Mesh();

			if (!pointCloudMode)
				currentFace->loadOff((char*)path.c_str());
			else
				currentFace->loadxyz((char*)path.c_str());

			int x = 400 * ((2 - 1) % 5);
			int y = -400 * ((2 - 1) / 5);

			currentFace->shiftMesh(x, y);

			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);

			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(face);
			mapper->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);
			
			
			renderer->AddActor(actor);
			renderer->Modified();
			
			faceNo++;

			currFace = currentFace;
			currPolyData = face;

			style->SetPolyData(face);
			//style->Modified();
			renderWindowInteractor->SetInteractorStyle(style);
			//renderWindowInteractor->Start();
			//renderWindowInteractor->Modified();


			
			segments.clear();
			for (int i = 0; i < currentFace->verts.size(); i++){

				segments.push_back(-1);
				isSelected[i] = false;
				

			}
			
			this->GetInteractor()->GetRenderWindow()->Render();

		}

		if (key == "d"){
			
			
			//if (!isFirst){
				//renderer->RemoveActor(renderer->GetActors()->GetLastActor());
				//renderer->RemoveActor(renderer->GetActors()->GetLastActor());
			//}
			//isFirst = false;
			string path;
			if (!pointCloudMode)
				path= "faces\\face" + std::to_string(faceNo) + ".off";
			else
				path = "faces\\face" + std::to_string(faceNo) + ".xyz";

			Mesh *currentFace = new Mesh();

			if (!pointCloudMode)
				currentFace->loadOff((char*)path.c_str());
			else 
				currentFace->loadxyz((char*)path.c_str());
			


			currentFace->findNeighborhoodTriangles();
			currentFace->assignNormalsToTriangles();
			//initialFace->calculateDihedrals();
			//numOfSegs = initialFace->smoothAllNormals(10000);
			//currentFace->angleXY();

			int x = 400 * ((faceNo) % 5);
			int y = -900 * ((faceNo) / 5);

			currentFace->shiftMesh(x,y-500);

		
			
			vector<int> currSegments;

			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);
		
			int nVerts = face->GetPoints()->GetNumberOfPoints();
			vtkDataArray * s = face->GetPointData()->GetScalars();

			ifstream ff("faces\\corr"+to_string(faceNo)+".txt");
			ofstream fff("segmentation_" + name + "_" + to_string(faceNo) + ".txt");
			int temp = segments.size();
			fff << "13" << std::endl;
			for (int i = 0; i < currentFace->verts.size(); i++){

				int cor;
				
				ff >> cor;
				//ff >> diff;
				//cout << cor << endl;

				fff << segments[cor] << std::endl;
				currSegments.push_back(segments[cor]);
				
				
				if (segments[cor] != -1 ){
					s->SetTuple1(i, colorPerSegment*segments[cor]);
				}
				else
					s->SetTuple1(i, 360);
			}
			fff.close();
			face->GetPointData()->SetScalars(s);
			face->GetPointData()->GetScalars()->Modified();
		
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(face);
			mapper->SetScalarRange(0, 360);
			
			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);
		
			renderer->AddActor(actor);
			renderer->Modified();
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			Mesh *currentFace2 = new Mesh();

			if (!pointCloudMode)
				currentFace2->loadOff((char*)path.c_str());
			else
				currentFace2->loadxyz((char*)path.c_str());

			currentFace2->findNeighborhoodTriangles();
			currentFace2->assignNormalsToTriangles();

			currentFace2->shiftMesh(x, y);

			vtkPolyData *face2;
			//cube->Delete();
			face2 = currentFace2->getVTKPolyData(!pointCloudMode);

			
			vtkDataArray * s2 = face2->GetPointData()->GetScalars();


			for (int i = 0; i < currentFace->verts.size(); i++){


				if (segments[i] != -1){
					s2->SetTuple1(i, colorPerSegment*segments[i]);
				}
				else
					s2->SetTuple1(i, 360);

			}

			face2->GetPointData()->SetScalars(s2);
			face2->GetPointData()->GetScalars()->Modified();

			vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper2->SetInputData(face2);
			mapper2->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
			actor2->SetMapper(mapper2);
			actor2->GetProperty()->SetPointSize(5);

			renderer->AddActor(actor2);
			renderer->Modified();

			//////////////////////////////////////////////////////////////////////////////////////////////////
			this->GetInteractor()->GetRenderWindow()->Render();
			faceNo++;
			
			/*pair<float,float> hd = HammingDistance(currentFace, segments, currSegments);
			float randIndexError = randIndex(segments, currSegments);
			pair<float, float> ce = consistencyError(currentFace, segments, currSegments);
			float cutDiscError = cutDiscrepancy(currentFace, segments, currSegments) + cutDiscrepancy(currentFace, currSegments, segments);

			
			
			cout << "HD-RM	 :" << hd.first << endl;
			cout << "HD-RF	 :" << hd.second << endl;
			cout << "HD		 :" << (hd.first + hd.second) / 2 << endl;
			cout << "RI		 :" << randIndexError << endl;
			cout << "GCE	 :" << ce.first << endl;
			cout << "LCE	 :" << ce.second << endl;
			cout << "CutDisc :" << cutDiscError << endl;


			
			std::cout << "--------------------------" << std::endl;


			ofstream of("Error.txt", fstream::app);
			of << hd.first << " " <<
				hd.second << " " <<
				(hd.first + hd.second) / 2 << " " <<
				randIndexError << " " <<
				ce.first << " " <<
				ce.second << " " <<
				cutDiscError << endl;

			of.close();
			*/
		}

		if (key == "Shift_L"){

			alignDatabase();

		}


		if (key == "u"){

		
			ifstream f;
			f.open("query/mount1.xyz");
			float x, y, z;
			vector<Vertex*> vertices;

			for (int i = 0; i < 37; i++){

				float* c = new float[3];
				f >> c[0] >> c[1] >> c[2];

				Vertex* v = new Vertex(i, c);
				vertices.push_back(v);

			}

		}

	
		
		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

	void SetPolyData(vtkSmartPointer<vtkPolyData> polyData) { this->PolyData = polyData; }
private:
	vtkSmartPointer<vtkExtractPolyDataGeometry> extractPolyDataGeometry;
	vtkSmartPointer<vtkPolyData> PolyData;
	vtkSmartPointer<vtkActor> SelectedActor;
	vtkSmartPointer<vtkDataSetMapper> SelectedMapper;

};
vtkStandardNewMacro(HighlightInteractorStyle);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill(vector<Vertex*> verts){














}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vtkSmartPointer<vtkActor> getVectorActor(double* v, double* startPoint,double* color){


	//Create an arrow.
	vtkSmartPointer<vtkArrowSource> arrowSource =
		vtkSmartPointer<vtkArrowSource>::New();

	vtkSmartPointer<vtkMatrix4x4> matrix =
		vtkSmartPointer<vtkMatrix4x4>::New();

	// Create the direction cosine matrix
	matrix->Identity();
	for (unsigned int i = 0; i < 3; i++)
	{
		matrix->SetElement(i, 0, v[i]);
		//matrix->SetElement(i, 1, normalizedY[i]);
		//matrix->SetElement(i, 2, normalizedZ[i]);
	}

	// Apply the transforms
	vtkSmartPointer<vtkTransform> transform =
		vtkSmartPointer<vtkTransform>::New();
	transform->Translate(startPoint);
	transform->Concatenate(matrix);
	transform->Scale(100,100,100);

	// Transform the polydata
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPD =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	
	transformPD->SetTransform(transform);
	
	transformPD->SetInputConnection(arrowSource->GetOutputPort());

	//Create a mapper and actor for the arrow
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();

	mapper->SetInputConnection(arrowSource->GetOutputPort());
	actor->SetUserMatrix(transform->GetMatrix());
	actor->SetMapper(mapper);

	actor->GetProperty()->SetColor(color[0],color[1],color[2]);
	return actor;

}



void showPCA(vtkSmartPointer<vtkRenderer> renderer, Mesh* mesh, double* startPoint){

	Vector com1(0, 0, 0);
	
	for (int i = 0; i < mesh->verts.size(); i++){


		com1.X += mesh->verts[i]->coords[0];
		com1.Y += mesh->verts[i]->coords[1];
		com1.Z += mesh->verts[i]->coords[2];
	}

	
	com1 = com1.Multiply(1.0f / mesh->verts.size());

	int k = 3;
	int n2 = mesh->verts.size();
	MatrixXf ms2(k, n2);
	for (int i = 0; i < k; i++) //for all k rows
	for (int j = 0; j < n2; j++) //fill n2 columns
		ms2(i, j) = mesh->verts[j]->coords[i] - com1[i];
	MatrixXf mst2 = ms2.transpose(); //transpose of ms2
	Matrix3f cov2 = ms2 * mst2;
	SelfAdjointEigenSolver<Matrix3f> eigensolver2(cov2); //solve eigenvalue/vector for symmetric n by n matrices, a.k.a. selfadjoint matrices
	if (eigensolver2.info() != Success) abort();
	Vector3f eig_vals = eigensolver2.eigenvalues(); //given in ascending order
	Matrix3f principal_axes = eigensolver2.eigenvectors();

	double** principalAxes2 = new double*[k];
	for (int j = 0; j < k; j++)
	{
		principalAxes2[j] = new double[k];
		for (int c = 0; c < k; c++)
			principalAxes2[j][c] = principal_axes(c, k - j - 1);
	}

	double colors[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	for (int i = 0; i < 3;i++)
	renderer->AddActor(getVectorActor(principalAxes2[i], startPoint,colors[i]));
}



pair<float,float> HammingDistance(Mesh* mesh1, vector<int> segments1, vector<int> segments2){

	mesh1->calculateAreas();
	int segSize = *std::max_element(segments1.begin(), segments1.end());
	vector<set<int>> segmentTris1, segmentTris2;
	
	for (int i = 0; i < segSize + 1; i++){

		set<int> temp1,temp2;

		segmentTris1.push_back(temp1);
		segmentTris2.push_back(temp2);
	}


	float totalArea1 = 0;
	float totalArea2 = 0;
	for (int i = 0; i < mesh1->tris.size(); i++){
		

		int v1 = mesh1->tris[i]->v1i;
		int v2 = mesh1->tris[i]->v2i;
		int v3 = mesh1->tris[i]->v3i;

		int v1s1 = segments1[v1];
		int v2s1 = segments1[v2];
		int v3s1 = segments1[v3];

		int v1s2 = segments2[v1];
		int v2s2 = segments2[v2];
		int v3s2 = segments2[v3];

		if (v1s1 != -1 && v1s1 == v2s1 && v1s1 == v3s1){

			segmentTris1[v1s1].insert(mesh1->tris[i]->idx);

			totalArea1 += mesh1->tris[i]->area;

		}

		if (v1s2 != -1 && v1s2 == v2s2 && v1s2 == v3s2){

			segmentTris2[v1s2].insert(mesh1->tris[i]->idx);

			totalArea2 += mesh1->tris[i]->area;

		}


	}



	float HD1 = 0;

	for (int i = 0; i < segSize + 1; i++){

		if (segmentTris1[i].size() == 0)
			continue;

		float areaDiff = 0;
		std::vector<int> v;
		std::set_difference(segmentTris1[i].begin(), segmentTris1[i].end(), segmentTris2[i].begin(), segmentTris2[i].end(), std::inserter(v, v.begin()));
		for (auto it = v.begin(); it != v.end(); ++it){
			areaDiff += mesh1->tris[*it]->area;
		}

		//std::cout << areaDiff << std::endl;

		HD1 += areaDiff;
	}

	cout << "HD1:	" << HD1 / totalArea1 << endl;

	float HD2 = 0;

	for (int i = 0; i < segSize + 1; i++){

		if (segmentTris2[i].size() == 0)
			continue;

		float areaDiff = 0;
		std::vector<int> v;
		std::set_difference(segmentTris2[i].begin(), segmentTris2[i].end(), segmentTris1[i].begin(), segmentTris1[i].end(), std::inserter(v, v.begin()));
		for (auto it = v.begin(); it != v.end(); ++it){
			areaDiff += mesh1->tris[*it]->area;
		}

		//std::cout << areaDiff << std::endl;

		HD2 += areaDiff;
	}

	cout << "HD2:	" << HD2 / totalArea2 << endl;

	return make_pair<float, float>(HD1 / totalArea1, HD2 / totalArea2);

}


float randIndex(vector<int> segments1, vector<int> segments2){

	float rand = 0;
	int totalPairs = 0;

	for (int i = 0;i< segments1.size(); i++){

		for (int j = i + 1; j < segments1.size(); j++){

			if (segments1[i] == -1 && segments2[i] == -1)
				continue;

			int c = (segments1[i] == segments1[j]) ? 1 : 0;
			int p = (segments2[i] == segments2[j]) ? 1 : 0;

			rand += c*p + (1 - c)*(1 - p);
			totalPairs++;


		}


	}

	return 1 - (rand / totalPairs);

}


float localRefError(Mesh* mesh, vector<set<int>> segmentTris1, vector<set<int>> segmentTris2,int faceId) {

	int seg1 = -1;
	int seg2 = -1;

	for (int i = 0; i < segmentTris1.size(); i++){
		if (std::find(segmentTris1[i].begin(), segmentTris1[i].end(), faceId) != segmentTris1[i].end()){

			seg1 = i;
			break;
		}
	}

	if (seg1 == -1)
		return -1;

	for (int i = 0; i < segmentTris2.size(); i++){
		if (std::find(segmentTris2[i].begin(), segmentTris2[i].end(), faceId) != segmentTris2[i].end()){

			seg2 = i;
			break;
		}
	}

	if (seg2 == -1)
		return - 1;

	set<int> R_S1 = segmentTris1[seg1];
	set<int> R_S2 = segmentTris2[seg2];

	float areaDiff = 0;
	std::vector<int> v;
	std::set_difference(R_S1.begin(), R_S1.end(), R_S2.begin(), R_S2.end(), std::inserter(v, v.begin()));
	for (auto it = v.begin(); it != v.end(); ++it){
		areaDiff += mesh->tris[*it]->area;
	}

	float areaTotal = 0;
	for (auto it = R_S1.begin(); it != R_S1.end(); ++it){
		areaTotal += mesh->tris[*it]->area;
	}


	return areaDiff / areaTotal;


}


pair<float,float> consistencyError(Mesh* mesh, vector<int> segments1, vector<int> segments2){


	mesh->calculateAreas();
	int segSize = *std::max_element(segments1.begin(), segments1.end());
	vector<set<int>> segmentTris1, segmentTris2;

	for (int i = 0; i < segSize + 1; i++){

		set<int> temp1, temp2;

		segmentTris1.push_back(temp1);
		segmentTris2.push_back(temp2);
	}

	for (int i = 0; i < mesh->tris.size(); i++){


		int v1 = mesh->tris[i]->v1i;
		int v2 = mesh->tris[i]->v2i;
		int v3 = mesh->tris[i]->v3i;

		int v1s1 = segments1[v1];
		int v2s1 = segments1[v2];
		int v3s1 = segments1[v3];

		int v1s2 = segments2[v1];
		int v2s2 = segments2[v2];
		int v3s2 = segments2[v3];

		if (v1s1 != -1 && v1s1 == v2s1 && v1s1 == v3s1){

			segmentTris1[v1s1].insert(mesh->tris[i]->idx);

		}

		if (v1s2 != -1 && v1s2 == v2s2 && v1s2 == v3s2){

			segmentTris2[v1s2].insert(mesh->tris[i]->idx);

		}


	}

	vector<float> refError1, refError2;

	for (int i = 0; i < mesh->tris.size(); i++){

		float E1 = localRefError(mesh, segmentTris1, segmentTris2, i);
		float E2 = localRefError(mesh, segmentTris2, segmentTris1, i);

		if (E1 < 0 || E2 < 0)
			continue;

		refError1.push_back(E1);
		refError2.push_back(E2);

	}

	float GCE;
	float GCE1 = 0;
	float GCE2 = 0;
	float LCE = 0;

	for (int i = 0; i < refError1.size(); i++){

		GCE1 += refError1[i];
		GCE2 += refError2[i];
		LCE += min(refError1[i], refError2[i]);

	}


	GCE = min(GCE1, GCE2) / refError1.size();
	LCE /= refError1.size();

	pair<float, float> ce;
	ce.first = GCE;
	ce.second = LCE;
	return ce;

}


Mesh* getOnlySegmented(Mesh* mesh, vector<int> segments, vector<int> &oldIndices){


	Mesh* newMesh = new Mesh();

	
	for (int i = 0; i < mesh->verts.size(); i++){

		if (segments[i] != -1){

			oldIndices.push_back(i);
			newMesh->addVertex(mesh->verts[i]->coords);
		}
	}

	for (int i = 0; i < mesh->tris.size(); i++){
	
		int v1 = mesh->tris[i]->v1i;
		int v2 = mesh->tris[i]->v2i;
		int v3 = mesh->tris[i]->v3i;

		if (segments[v1] != -1 && segments[v2] != -1 && segments[v3] != -1){

			int v1n = find(oldIndices.begin(), oldIndices.end(), v1) - oldIndices.begin();
			int v2n = find(oldIndices.begin(), oldIndices.end(), v2) - oldIndices.begin();
			int v3n = find(oldIndices.begin(), oldIndices.end(), v3) - oldIndices.begin();

			newMesh->addTriangle(v1n, v2n, v3n);


		}
	
	}

	return newMesh;
}

float cutDiscrepancy(Mesh* mesh, vector<int> segments1, vector<int> segments2){

	
	vector<int> segBound1, segBound2;


	for (int i = 0; i < segments1.size(); i++){

		if (segments1[i] == -1)
			continue;

		bool isInside = true;
		for (int j = 0; j < mesh->verts[i]->vertList.size(); j++){

			int n = mesh->verts[i]->vertList[j];

			isInside = isInside && (segments1[i] == segments1[n]);


		}

		if (!isInside)
			segBound1.push_back(i);


	}

	for (int i = 0; i < segments2.size(); i++){

		if (segments2[i] == -1)
			continue;

		bool isInside = true;
		for (int j = 0; j < mesh->verts[i]->vertList.size(); j++){

			int n = mesh->verts[i]->vertList[j];

			isInside = isInside && (segments2[i] == segments2[n]);


		}

		if (!isInside)
			segBound2.push_back(i);


	}
	

	Dijsktra::DijsktraSP sp;
	sp.setMesh(mesh);
	
	int i = 0;
	float avDist = 0;
	int count = 0;

	for (int i = 0; i < segBound1.size(); i++){

		float minDist = 9999;
		int v2index;
		int v1 = segBound1[i];// find(oldIndices.begin(), oldIndices.end(), segBound1[i]) - oldIndices.begin();
		
		//if (v1 == oldIndices.size())
			//continue;

		vector<Dijsktra::Node*> r = sp.run(v1);
		

		for (int j = 0; j < segBound2.size(); j++){

			int v2 = segBound2[j];// find(oldIndices.begin(), oldIndices.end(), segBound2[j]) - oldIndices.begin();
			
			//if (v2 == oldIndices.size())
				//continue;


			float cost = dist2Between(mesh->verts[v1]->coords, mesh->verts[v2]->coords) + r[v2]->key;
			if (cost < minDist){
				minDist = cost;
				v2index = v2;
			}

		}
		//cout << v1 << "  -->  " << v2index << "   =   " << minDist << endl;
	//	cout << onlySegments->verts[v1]->coords[0] << "  " << onlySegments->verts[v1]->coords[1] << " " << onlySegments->verts[v1]->coords[2] << endl;
	//	cout << onlySegments->verts[v2index]->coords[0] << "  " << onlySegments->verts[v2index]->coords[1] << " " << onlySegments->verts[v2index]->coords[2] << endl;
		//cout << minDist << endl;

		for (int k = 0; k < r.size(); k++){

			delete r[k];

		}


		avDist += minDist;
		count++;

	}
	
	avDist /= count;

	

	float avgRad = 0;
	int avgCount = 0;
	for (int i = 0; i < mesh->verts.size(); i++){


		if (segments1[i] == -1)
			continue;

		float xd = mesh->verts[i]->coords[0] - mesh->center(0);
		float yd = mesh->verts[i]->coords[1] - mesh->center(1);
		float zd = mesh->verts[i]->coords[2] - mesh->center(2);

		avgRad += sqrt(xd*xd + yd*yd + zd*zd);
		avgCount++;

	}
	
	avgRad /= avgCount;
	avDist /= avgRad;
	return avDist;


	

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{



	//std::cout << "give user name" << std::endl;
	//std::cin >> name;
	std::cin >> faceNo;
	
	name = "initial";
	initialFace = new Mesh();

	if (!pointCloudMode)
		initialFace->loadOff("faces\\face1.off");
	else 
		initialFace->loadxyz("faces\\face1.xyz");

	 s = initialFace->calculateScale();
	 cube = initialFace->getVTKPolyData(!pointCloudMode);
	int nVerts = cube->GetPoints()->GetNumberOfPoints();

	
	isSelected = new bool[nVerts];
	segments.clear();
	for (int i = 0; i < nVerts; i++){
		isSelected[i] = false;
		segments.push_back(-1);

	}

	initialFace->findNeighborhoodTriangles();
	initialFace->assignNormalsToTriangles();
	//initialFace->calculateDihedrals();
	//numOfSegs = initialFace->smoothAllNormals(10000);
	//initialFace->angleXY();

	currFace = initialFace;
	
	vtkDataArray* s = cube->GetPointData()->GetScalars();

	for (int i = 0; i < initialFace->verts.size(); i++)
		s->SetTuple1(i,360);

	cube->GetPointData()->SetScalars(s);
	cube->GetPointData()->GetScalars()->Modified();
	
	currPolyData = cube;

	 mapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
   mapper->SetInputData(cube);
  mapper->SetScalarRange(0, 360);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(5);
  
  // Visualize
  renderer =
	  vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
	  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkAreaPicker> areaPicker =
	  vtkSmartPointer<vtkAreaPicker>::New();
   renderWindowInteractor =
	  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetPicker(areaPicker);
  renderWindowInteractor->SetRenderWindow(renderWindow);


  /////////////////////////////////////////////////////////////////////////////////////////
 /* Mesh* mesh = new Mesh();
  mesh->loadOff("faces\\fine2.off");
  float d1 = mesh->calculateScale();
  mesh->shiftMesh(200, 0);
  vtkPolyData *cube2 = mesh->getVTKPolyData();
  nVerts = cube2->GetPoints()->GetNumberOfPoints();
  mesh->shiftMesh(-200, 0);
  vector<float> distances = calcClosests(initialFace, mesh);

  float min = 99999999;
  float max = 0;

  for (int i = 0; i < distances.size(); i++){

	  if (distances[i] > max)
		  max = distances[i];

	  if (distances[i] < min)
		  min = distances[i];
  }

  for (int i = 0; i < distances.size(); i++)
	  distances[i] -= min;


  vtkDataArray* s2 = cube2->GetPointData()->GetScalars();

  for (int i = 0; i < mesh->verts.size(); i++)
	  s2->SetTuple1(i,  distances[i]);

  cube2->GetPointData()->SetScalars(s2);
  cube2->GetPointData()->GetScalars()->Modified();



  vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper2->SetInputData(cube2);
  mapper2->SetScalarRange(0, 75.9);
  


  vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
  actor2->SetMapper(mapper2);
  actor2->GetProperty()->SetPointSize(5);
  renderer->AddActor(actor2);
  ////////////////////////////////////////////////////////////////////////////////////////////
  
  Mesh* mesh2 = new Mesh();
  mesh2->loadOff("faces\\deformed2.off");
  mesh2->shiftMesh(400, 0);
  vtkPolyData *cube3 = mesh2->getVTKPolyData();
  nVerts = cube3->GetPoints()->GetNumberOfPoints();
  mesh2->shiftMesh(-400, 0);
  distances = calcClosests(initialFace, mesh2);

   min = 99999999;
   max = 0;

  for (int i = 0; i < distances.size(); i++){

	  if (distances[i] > max)
		  max = distances[i];

	  if (distances[i] < min)
		  min = distances[i];
  }

  for (int i = 0; i < distances.size(); i++)
	  distances[i] -= min;


   s2 = cube3->GetPointData()->GetScalars();

  for (int i = 0; i < mesh2->verts.size(); i++)
	  s2->SetTuple1(i, distances[i]);

  cube3->GetPointData()->SetScalars(s2);
  cube3->GetPointData()->GetScalars()->Modified();



  vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper3->SetInputData(cube3);
  mapper3->SetScalarRange(0, 75.9);



  vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
  actor3->SetMapper(mapper3);
  actor3->GetProperty()->SetPointSize(5);
  renderer->AddActor(actor3);

  ///////////////////////////////////////////////////////////////////////////////////////////
  double startPoint[3] = { 0, 200, 0 };
  double startPoint2[3] = { 200, 200, 0 };
 // showPCA(renderer, initialFace, startPoint);
  //showPCA(renderer, mesh, startPoint2);
  
  
  vtkSmartPointer<vtkScalarBarActor> scalarBar =
	  vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(mapper2->GetLookupTable());
  
  double pos[2] = { 100, 100};
  scalarBar->SetNumberOfLabels(4);
  scalarBar->SetWidth(scalarBar->GetWidth()/3);
  scalarBar->SetHeight(scalarBar->GetHeight()*0.75);
 
  //scalarBar->SetLabelTextProperty();

  double cc[3] = { 0, 0, 1 };
  vtkTextProperty* tp =  scalarBar->GetLabelTextProperty();
  tp->SetFontSize(100);
  

  tp->SetColor(cc);
  //tp->Modified();
  scalarBar->SetLabelTextProperty(tp);
  //scalarBar->Modified();


  // Create a lookup table to share between the mapper and the scalarbar
  vtkSmartPointer<vtkLookupTable> hueLut =
	  vtkSmartPointer<vtkLookupTable>::New();
  hueLut->SetTableRange(0, 75.9);
  //hueLut->SetHueRange(0, 1);
  //hueLut->SetSaturationRange(1, 1);
  //hueLut->SetValueRange(0, 1);
  hueLut->Build();

  mapper2->SetLookupTable(hueLut);
  scalarBar->SetLookupTable(hueLut);

 renderer->AddActor2D(scalarBar);*/
 /////////////////////////////////////////////////////////////////////////////////////////////

renderer->AddActor(actor);
  renderer->SetBackground(1,1,1);
  
  vtkSmartPointer<vtkAxesActor> axes =
	  vtkSmartPointer<vtkAxesActor>::New();

  vtkSmartPointer<vtkOrientationMarkerWidget> widget =
	  vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  widget->SetOrientationMarker(axes);
  widget->SetInteractor(renderWindowInteractor);
  widget->SetViewport(0.0, 0.0, 0.4, 0.4);
  widget->SetEnabled(0);
 
 
 // renderer->ResetCamera();
  renderWindow->Render();

  style =  vtkSmartPointer<HighlightInteractorStyle>::New();
  style->SetPolyData(cube);
  renderWindowInteractor->SetInteractorStyle(style);
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
