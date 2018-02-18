
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
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
#include <queue>

#define NUM_CONTOUR_POINTS 40

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
int numOfSegments = 8;
float colorPerSegment = 360 / numOfSegments;
Mesh *initialFace, *currFace;
bool pointCloudMode = false;
int numOfSegs;
string name;
vtkSmartPointer<vtkScalarBarActor> scalarBar;

//vector<int>  mouthPoints = { 553, 2239, 552, 4642, 559, 2252, 558, 4868, 4867, 4915, 4916, 5188, 3423, 3053, 936, 3067, 943, 3069, 944, 1850, 1854, 375, 1840, 377, 2767, 802, 3957, 1164, 4378, 1849, 3143, 3171, 4375, 3117, 3114, 3115, 4909, 4910, 6121, 4958, 4934, 2248, 6124, 1685, 5712, 5603, 5602 ,553};
//vector<int>   nosePoints  = { 4988, 5480, 5481, 5485, 5474, 5473, 5475, 5564, 5562, 5561, 5395, 5394, 2541, 2539, 2540, 4624, 4626, 5400, 5397, 5396, 5399, 4590, 4588, 4592, 5159, 3382, 3381, 2787, 2786, 2789, 3819, 3636, 3635, 3638, 3639, 2824, 2823, 2135, 2132, 500, 2812, 828, 2846, 845, 3805, 1116, 3631, 760, 2697, 763, 3034, 364, 3036, 361, 1807, 4988 };


vector<int>  mouthPoints = { 4749,4874,5000,5512,6024,6538,7054,7699,8215,8602,9118,9892,10408,10799,11062,11325,11585,11201,10308,9793,9148,8374,7600,7084,6437,5919,5400 ,4749};
vector<int>   eye1 = {2088,2478,2865,3124,3385,3772,4159,4545,4932,5320,5706,6091,6215,5954,5567,5180,4790,4403,4016,3758,3371,2985,2599,2214,2088 };
vector<int>   nosePoints = { 6981,6986,6993,7000,7011,7019,7026,6252,5607,5614,5749,5751,6525,7170,7815,8718,9621,10653,10907,11158,11026,9736,9731,9725,9719,9711,9705,9699,9692,9689,8915,8272,7759,6981};
vector<int>   eye2 = {10349,10998,11772,12417,13191,13706,14220,14733,14467,13691,13174,12658,12142,11497,10983,10726,10212,10349};
vector<int>  chick1 = {5200,5216,5235,4203,3171,3185,3457,3857,4254,46271,46509,46737,46049,45187,44020,23129,22471,22196,22430,22418,21757,22913,24200,25746,27423,2746,5200};
vector<int>  chick2 = {11523,11531,11539,11547,11559,12333,13107,13115,13124,12877,12504,50153,50005,49792,50667,51568,52569,33192,33439,33425,32756,32869,32859,31954,31049,30015,29240,28594,27690,15910,15136,14361,13199,11523};
vector<int> chin = { 5675,46794,46798,46807,47673,48258,48778,49366,49357,49353,10703,9671,8639,7091,5672,5675};
vector<vector<int>> faceParts = { mouthPoints, nosePoints, eye1, eye2, chick1, chick2, chin };


vector<int>  mouthPoints2 = { 1627,553,552,559,558,1509,1390,833,3422,3067,3069,944,1850,2844,375,377,3958,1846,3171,3114,1507,4933,4960,5714,5603,1627 };
vector<int>   eye12 = { 527,4820,1601,4713,1431,4718,4717,4723,4726,4730,4733,4734,4735,4741,2386,2383,4809,4808,5502,4806,4804,2209,2207,2211,527 };
vector<int>   nosePoints2 = { 4794,4835,4502,5393,5562,5395,2541,4624,5572,2581,2548,4581,809,2816,2585,2147,2769,3816,832,2823,500,828,845,3716,3714,3715,3721,3201,1806,2992,3775,3026,3024,4824,5532,4792,4794 };
vector<int>   eye22 = { 1975,2942,2937,2935,2933,2931,2926,2924,2919,2917,2915,2912,2913,3020,2152,1815,1811,1813,3005,3006,3743,3011,3008,3010,3019,1974,1973,1975 };
vector<int>  chick12 = {4994,5650,5627,5349,5348,4706,1633,5470,5467,4688,5388,5406,4508,4510,1340,5423,5303,4433,1289,5655,605,1677,5049,613,2409,4994  };
vector<int>  chick22 = { 3208,3894,3871,3593,3615,3679,3711,3745,800,497,3644,767,2702,2713,1122,791,2744,3567,2621,3899,1966,4197,1018,3306,1996,3208 };
vector<int> chin2 = { 5717,5463,4493,2684,3700,3704,3694,3674,446,2020,1951,2692,784,4499,2366,2432,5429,5719,5717 };
vector<vector<int>> faceParts2 = { mouthPoints2, nosePoints2, eye12, eye22, chick12, chick22, chin2 };
vector<vector<int>> completedFaceParts;
vector<int> facePoints;

void search(vector<graph::Vertex*> verts, vector<Vector*> path, vector<vector<graph::Edge*>> &allCandidates, vector<float> &costs);

pair<float,float> HammingDistance(Mesh* mesh1, vector<int> segments1, vector<int> segments2);
float randIndex(vector<int> segments1, vector<int> segments2);
pair<float, float> consistencyError(Mesh* mesh, vector<int> segments1, vector<int> segments2);
Mesh* getOnlySegmented(Mesh* mesh, vector<int> segments, vector<int> &oldIndices);
float cutDiscrepancy(Mesh* mesh, vector<int> segments1, vector<int> segments2);

vector<vector<int>> completeContoursByGeodesic(Mesh* cF, vector<vector<int>> inFaceParts);
vector<vector<int>> findInnerContourVertices(Mesh* cF, vector<vector<int>> inFaceParts);
vector<int**> extractContourFeatures(Mesh *face, vector<vector<int>> faceParts);
vector<float*> extractContourFeatures2(Mesh *face, vector<vector<int>> faceParts);
float matchFaceParts(vector<int**> f1, vector<int**> f2);

float matchContour(int**, int**);
float matchFaceParts(float* f1, float* f2);

float triangleArea(vector<Vertex*> vertices, std::vector<int> hole, int v0, int v1, int v2);
float holeCostFunc(vector<Vertex*> vertices, std::vector<int> hole, int i, int k, vector<int> &lamda);
void trace(vector<Vertex*> vertices, std::vector<int> hole, int i, int k, vector<int> &lamda, vector<int> &newTriangles);

vector<Vertex*> completeContour(vector<Vertex*> contourPoints);
float nonRigidAlignment(vector<Vertex*> &pntSet1, vector<Vertex*> &pntSet2, int k);
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
	/*averageFace = new Mesh();

	if(!pointCloudMode)
		averageFace->loadOff("faces\\face1.off");
	else 
		averageFace->loadxyz("faces\\face1.xyz");
	*/
	//s = averageFace->calculateScale();

	/*
	Mesh *initialFace = new Mesh();

	if (!pointCloudMode)
		initialFace->loadOff("faces\\face1.off");
	else
		initialFace->loadxyz("faces\\face1.xyz");
	*/

	ofstream of("faces\\corr1.txt");

	for (int i = 0; i < initialFace->verts.size(); i++)
		of << i << std::endl;

	of.close();


	int i = 2;
	string s1 = "faces\\face";
	string s2 = "faces\\corr";
	string s3 = "faces\\deformed";
	string s4 = "faces\\coarse";
	string s5 = "faces\\fine";
	
	float tt = 0;

	for (int face = 2; face <= 100; face++){
		std::cout << "			" << face << std::endl;
		averageFaceTemp = new Mesh();
		if (!pointCloudMode)
			averageFaceTemp->loadOff("faces\\face1.off");
		else
			averageFaceTemp->loadxyz("faces\\face1.xyz");

		//averageFaceTemp->scale(0.001);
		//averageFaceTemp->write("faces\\scaled1.off");
		//averageFaceTemp->initialize();
		

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
		
		//mesh->scale(0.001);
		//mesh->write((char*)("faces\\scaled" + to_string(face) + ".off").c_str());
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

		//PCAAlignment(averageFaceTemp, mesh);
		//mesh->write((char*)path4.c_str(), !pointCloudMode);

		

		//initialFineRigidAlignment(averageFaceTemp, mesh);
	
		//mesh->write((char*)path5.c_str(), !pointCloudMode);

		mesh->initialize();
		mesh->deform(averageFaceTemp, 1000, 50);
		mesh->write((char*)path3.c_str(), !pointCloudMode);
		

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
		delete mesh;
	}


	//std::cout << tt/19 <<std::endl ;
	

	//averageFace->write("averageFaceInc.off");
}
vector<float*> query;
Mesh* input;
bool isFirst = true;


void createInputFile(Mesh* m, vector<vector<int>> indexVector){

	ofstream off("inputPoints.xyz");
	for (int fp = 0; fp < indexVector.size(); fp++){
		for (int p = 0; p < indexVector[fp].size(); p++){

			int vid = indexVector[fp][p];
			float x, y, z;
			x = m->verts[vid]->coords[0];
			y = m->verts[vid]->coords[1];
			z = m->verts[vid]->coords[2];
			off << x << " " << y << " " << z << endl;

		}
	}

	off.close();

}















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

		if (key == "j"){
			
			//vtkDataArray * s = currPolyData->GetPointData()->GetScalars();

			ofstream of("inputPoints.xyz");
			for (int fp = 0; fp < faceParts.size(); fp++){
				for (int p = 0; p < faceParts[fp].size(); p++){

					int vid = faceParts[fp][p];
					float x, y, z;
					x = initialFace->verts[vid]->coords[0];
					y = initialFace->verts[vid]->coords[1];
					z = initialFace->verts[vid]->coords[2];
					of << x <<" "<< y <<" "<< z << endl;

					//s->SetTuple1(vid, fp * 50 + 10 );

				}
			}

			of.close();


			//currPolyData->GetPointData()->SetScalars(s);
			//currPolyData->GetPointData()->GetScalars()->Modified();
			//this->GetInteractor()->GetRenderWindow()->Render();
			//this->HighlightProp(NULL);

		}

		if (key == "h"){


			vector<vector<int>> innerPoints = { {}, {}, {}, {}, {}, {}, {} };
			vtkDataArray * s = currPolyData->GetPointData()->GetScalars();
			std::ifstream f("segmented2.xyz");

			for (int i = 0; i < initialFace->verts.size(); i++){

				float x, y, z;
				int   r, g, b;
				f >> x >> y >> z;
				f >> r >> g >> b;

				int colorId = -1;

				if (r == 0 && g == 0 && b == 0) colorId = 0;  // mouth
				if (r == 255 && g == 255 && b == 0) colorId = 1;  // chin
				if (r == 255 && g == 0 && b == 255) colorId = 2;  // nose
				if (r == 0 && g == 255 && b == 255) colorId = 3;  // left eye
				if (r == 0 && g == 255 && b == 0) colorId = 4;  // left chick
				if (r == 0 && g == 0 && b == 255) colorId = 5;  // right chick
				if (r == 255 && g == 0 && b == 0) colorId = 6;  // right eye

				segments[i] = colorId;
				isSelected[i] = colorId != -1;

			}
			f.close();
			for (int i = 0; i < initialFace->verts.size(); i++){

				int partId;
				int cid = segments[i];
				if (cid == -1)
					continue;


				if (cid == 0) partId = 0;
				if (cid == 1) partId = 6;
				if (cid == 2) partId = 1;
				if (cid == 3) partId = 2;
				if (cid == 4) partId = 4;
				if (cid == 5) partId = 5;
				if (cid == 6) partId = 3;

				innerPoints[partId].push_back(i);
				s->SetTuple1(i, 50 * partId + 10);
				
			}


			currPolyData->GetPointData()->SetScalars(s);
			currPolyData->GetPointData()->GetScalars()->Modified();
			this->GetInteractor()->GetRenderWindow()->Render();
			this->HighlightProp(NULL);
			
		}

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
			/*
			ofstream f;
			f.open("segmentation_"+name+"_"+to_string(faceNo)+".txt");

			f << c / colorPerSegment << std::endl;

			for (int i = 0; i < segments.size(); i++)
				f << segments[i] << std::endl;

			f.close();

			*/

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
		
			int i = 8;
			string meshPath = "faces\\face" + to_string(i) + ".off";
			input = new Mesh();
			input->loadOff((char*)meshPath.c_str());
			input->assignNormalsToTriangles();
			input->findNeighborhoodTriangles();
			input->calculateDihedrals();

			vtkPolyData* face = input->getVTKPolyData();

			vtkDataArray * s = face->GetPointData()->GetScalars();
			int nVerts = face->GetPoints()->GetNumberOfPoints();

			for (int i = 0; i < nVerts; i++) {


				if (input->verts[i]->dihedral > 20)
					s->InsertTuple1(i, 200);
				else s->InsertTuple1(i, 0);
			}


			face->GetPointData()->SetScalars(s);
			face->GetPointData()->GetScalars()->Modified();

			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(face);
			mapper->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);

			renderer->RemoveActor(renderer->GetActors()->GetLastActor());

			renderer->AddActor(actor);
			renderer->Modified();
			this->GetInteractor()->GetRenderWindow()->Render();

		}

		if (key == "l"){

			std::cout << "give user name" << std::endl;
			std::cin >> name;
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
			f.close();


			
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


			/*
			initialFace->shiftMesh(200, 0);
			vector<int> segment;
			for (int k = 0; k < segments.size(); k++)
				if (segments[k] == 0)
					segment.push_back(k);

			vtkPolyData *segmentPolyData = initialFace->getVTKPolyDataSubsetByVertices(segment);

			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(segmentPolyData);
			mapper->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);

			renderer->AddActor(actor);
			renderer->Modified();
			style->SetPolyData(segmentPolyData);
			*/
			
			renderWindowInteractor->SetInteractorStyle(style);
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
			faceNo++;
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

			int x = 100000 * ((faceNo) % 5);
			int y = -200000 * ((faceNo) / 5);

			currentFace->shiftMesh(x,y-100000);

		
			
			vector<int> currSegments;

			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);
		
			int nVerts = face->GetPoints()->GetNumberOfPoints();
			vtkDataArray * s = face->GetPointData()->GetScalars();

			ifstream ff("faces\\corr"+to_string(faceNo)+".txt");
			//ofstream fff("segmentation_" + name + "_" + to_string(faceNo) + ".txt");
			int temp = segments.size();
			//fff << "13" << std::endl;
			for (int i = 0; i < currentFace->verts.size(); i++){

				int cor;
				
				ff >> cor;
				//ff >> diff;
				//cout << cor << endl;

				//fff << segments[cor] << std::endl;
				currSegments.push_back(segments[cor]);
				
				
				if (segments[cor] != -1 ){
					s->SetTuple1(i, colorPerSegment*segments[cor]);
				}
				else
					s->SetTuple1(i, 360);
			}
			//fff.close();
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
			/*
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
			*/
			//////////////////////////////////////////////////////////////////////////////////////////////////
			this->GetInteractor()->GetRenderWindow()->Render();
			//faceNo++;
			
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

			string filename = "query\\random.xyz";
			std::ifstream inFile(filename);
			int size =  std::count(std::istreambuf_iterator<char>(inFile),
				std::istreambuf_iterator<char>(), '\n');

			ifstream f;
			f.open(filename);
			query.clear();
			for (int i = 0; i < size; i++){

				float* c = new float[3];
				f >> c[0] >> c[1] >> c[2];
				query.push_back(c);
	
			}

			std::vector<Vector*> path;
			for (int i = 1; i < query.size(); i++){

				path.push_back(new Vector(query[i - 1][0] - query[i][0],
					query[i - 1][1] - query[i][1],
					query[i - 1][2] - query[i][2])
					);

			}

			int matchId = -1;
			float globalMinCost = 999999999;
			vector<int> minVertices;

			for (int i = 12; i < 13;i++){

					if (i == 4 || i==11)
						continue;
					string meshPath =  "faces\\face" + to_string(i) + ".off";
					string outPath = "query\\out" + to_string(i) + ".xyz";
					 input = new Mesh();
					input->loadOff((char*) meshPath.c_str());

					vector<graph::Vertex*> vectorGraph = input->createVectorGraph();

					vector<vector<graph::Edge*>> allCandidates;
					vector<float> costs;
					search(vectorGraph, path, allCandidates, costs);

					for (int can = 0; can < allCandidates.size(); can++){

						if (allCandidates[can].size() > 10){
							cout << allCandidates[can].size() << " " << costs[can] << endl;

							
							string outPath = "candidate" + to_string(can) + ".xyz";
							std::ofstream fff(outPath);
							for (int i = 0; i < allCandidates[can].size(); i++){

								int ind = allCandidates[can][i]->from->id;
								fff << input->verts[ind]->coords[0] << " " << input->verts[ind]->coords[1] << " " << input->verts[ind]->coords[2] << std::endl;
								cout << "	" << allCandidates[can][i]->from->id;

							}

							int ind = allCandidates[can][allCandidates[can].size() - 1]->to->id;
							fff << input->verts[ind]->coords[0] << " " << input->verts[ind]->coords[1] << " " << input->verts[ind]->coords[2] << std::endl;
							cout << "	" << allCandidates[can][allCandidates[can].size() - 1]->to->id;
							fff.close();
						

						}

						
					}

					float minCost = 999999;
					int minIndex = -1;
					std::cout << i<<"  ";
					if (allCandidates.size()){

						
						std::cout << "match found with ";
						for (int j = 0; j < costs.size(); j++){

							int score = allCandidates[j].size() - size;
							score = abs(score);

							if (score < minCost){

								minCost = score;
								minIndex = j;
							}
							//cout << costs[i] << endl;
						}

						if (minCost < globalMinCost){

							globalMinCost = minCost;
							matchId = i;

							minVertices.clear();
							for (int i = 0; i < allCandidates[minIndex].size(); i++){

								minVertices.push_back(allCandidates[minIndex][i]->from->id);
						
							}

							minVertices.push_back(allCandidates[minIndex][allCandidates[minIndex].size() - 1]->to->id);
							

						}

						//std::cout << minIndex << "	is chosen" << std::endl;
						std::cout  << minCost << std::endl;
						std::ofstream fff(outPath);
						for (int i = 0; i < allCandidates[minIndex].size(); i++){

							int ind = allCandidates[minIndex][i]->from->id;
							fff << input->verts[ind]->coords[0] << " " << input->verts[ind]->coords[1] << " " << input->verts[ind]->coords[2] << std::endl;
							cout << "	" << allCandidates[minIndex][i]->from->id;

						}

						int ind = allCandidates[minIndex][allCandidates[minIndex].size() - 1]->to->id;
						fff << input->verts[ind]->coords[0] << " " << input->verts[ind]->coords[1] << " " << input->verts[ind]->coords[2] << std::endl;
						cout << "	" << allCandidates[minIndex][allCandidates[minIndex].size() - 1]->to->id;
						fff.close();

					}
					else
					std::cout << "no match" << std::endl;
		}
			
			if (matchId != -1){
				Mesh* bestFace = new Mesh();
				bestFace->loadOff((char*)("faces\\face" + to_string(matchId) + ".off").c_str());

				ifstream ff;
				ff.open("segmentation_" + name + "_" + to_string(matchId) + ".txt");

				int segSize;
				ff >> segSize;

				vector<int> part;
				for (int i = 0; i < bestFace->verts.size(); i++){

					int temp;
					ff >> temp;
					if (temp == 0)
						part.push_back(i);

				}
				ff.close();

				renderer->RemoveActor(renderer->GetActors()->GetLastActor());

				vtkPolyData *segmentPolyData = bestFace->getVTKPolyDataSubsetByVertices(part);

				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(segmentPolyData);
				mapper->SetScalarRange(0, 360);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetPointSize(5);

				renderer->AddActor(actor);
				renderer->Modified();
				style->SetPolyData(segmentPolyData);


				renderWindowInteractor->SetInteractorStyle(style);
				this->GetInteractor()->GetRenderWindow()->Render();
			}
			cout << "done" << endl;
		}
		
		if (key == "p"){
		
			string path;
			if (!pointCloudMode)
				path =  "faces\\face" + std::to_string(faceNo) + ".off";
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

			int x = 250000 * ((faceNo-1) % 5);
			int y = -250000 * ((faceNo - 1) / 5);

			currentFace->shiftMesh(x, y );



			vector<int> currSegments;

			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);

			
	        


			int nVerts = face->GetPoints()->GetNumberOfPoints();
			vtkDataArray * s = face->GetPointData()->GetScalars();

			for (int i = 0; i < currentFace->verts.size(); i++)
				s->SetTuple1(i, 360);

			ifstream ff("faces\\corr" + to_string(1) + ".txt");

			int temp = segments.size();
			vector<int> corrs;
			for (int i = 0; i < currentFace->verts.size(); i++){

				int cor;
				ff >> cor;
				corrs.push_back(cor);

			}

			
			vector<vector<int>> correspondedFaceParts;


			for (int i = 0; i < faceParts.size(); i++){

				vector<int> temp;
				for (int p = 0; p < faceParts[i].size(); p++){

					int refIndex = faceParts[i][p];
					int inputIndex = corrs[refIndex];
					//cout << refIndex << " " << inputIndex << std::endl;
					//s->SetTuple1(refIndex, 0);
					//s->SetTuple1(inputIndex, 360);
					temp.push_back(inputIndex);

				}

				correspondedFaceParts.push_back(temp);

			}
	
			completedFaceParts = completeContoursByGeodesic(currentFace, faceParts);
			
			vector<vector<int>> innerPoints = findInnerContourVertices(currentFace,completedFaceParts);
			
			

			for (int fp = 0; fp < innerPoints.size(); fp++){
				for (int p = 0; p < innerPoints[fp].size(); p++){

					s->SetTuple1(innerPoints[fp][p], 10+fp*50);
				}
			}

			for (int fp = 0; fp < faceParts.size(); fp++){
				for (int p = 0; p < completedFaceParts[fp].size(); p++){

					s->SetTuple1(completedFaceParts[fp][p], 100);

					
				}
			}

			/*
			int rs[7] = { 0, 255, 255,   0,   0,   0, 255 };
			int gs[7] = { 0, 255,   0, 255, 255,   0,   0 };
			int bs[7] = { 0,   0, 255, 255,   0, 255,   0 };
			std::ofstream off("segmented.xyz");
			for (int pp = 0; pp < nVerts; pp++){

				bool isContour = false;

				for (int fp = 0; fp < faceParts.size(); fp++){
					for (int p = 0; p < completedFaceParts[fp].size(); p++){

						if (pp == completedFaceParts[fp][p]){
							isContour = true;

							double *coords;
							coords = face->GetPoints()->GetPoint(pp);

							off << coords[0] << " " << coords[1] << " " << coords[2] << " " << rs[fp] << " " << gs[fp] << " " << bs[fp] << endl;

							break;
						}
						
					}
				}

				if (!isContour){

					double *coords;
					coords = face->GetPoints()->GetPoint(pp);

					off << coords[0] << " " << coords[1] << " " << coords[2] << " " << 255 << " " << 255 << " " << 255 << endl;
				}

			}

			off.close();
			*/
			s->Modified();
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
			this->GetInteractor()->GetRenderWindow()->Render();
			
			
			
			

			
			faceNo++;

		
		}

		if (key == "comma"){


			string path;
			if (!pointCloudMode)
				path = "faces\\face" + std::to_string(faceNo) + ".off";
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

			int x = 400 * ((faceNo - 1) % 5);
			int y = -400 * ((faceNo - 1) / 5);

			currentFace->shiftMesh(x, y);

			float avgEdgeLen = 0.0f;
			for (int i = 0; i < currentFace->edges.size(); i++)
				avgEdgeLen += currentFace->edges[i]->length;

			avgEdgeLen /= currentFace->edges.size();

			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);

			int nVerts = face->GetPoints()->GetNumberOfPoints();
			vtkDataArray * s = face->GetPointData()->GetScalars();

			ifstream ff("faces\\corr" + to_string(faceNo) + ".txt");

			int temp = segments.size();
			vector<int> corrs;
			for (int i = 0; i < currentFace->verts.size(); i++){

				int cor;
				ff >> cor;
				corrs.push_back(cor);

			}

			

			vector<float> dists;
			vector<float*> distColors;
			vector<Vertex*> verts2 = currentFace->verts;
			float maxDist=0;
			for (int i = 0; i < nVerts; i++){

				face->GetPoint(i);
				float dist = sqrt(pow(verts2[i]->coords[0] - verts2[corrs[i]]->coords[0], 2)
					+ pow(verts2[i]->coords[1] - verts2[corrs[i]]->coords[1], 2)
					+ pow(verts2[i]->coords[2] - verts2[corrs[i]]->coords[2], 2));

				dist /= avgEdgeLen;

				if (dist > maxDist)
					maxDist = dist;

				dists.push_back(dist);

			}

			for (int i = 0; i < initialFace->verts.size(); i++)
				s->SetTuple1(i,dists[i]/maxDist);

			face->GetPointData()->SetScalars(s);
			face->GetPointData()->GetScalars()->Modified();


			

			// Create a lookup table to share between the mapper and the scalarbar
			vtkSmartPointer<vtkLookupTable> hueLut =
				vtkSmartPointer<vtkLookupTable>::New();
			hueLut->SetTableRange(0, maxDist);
			hueLut->SetHueRange(0, 1);
			hueLut->SetSaturationRange(1, 1);
			hueLut->SetValueRange(0, 1);
			hueLut->Build();


			vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper2->SetLookupTable(hueLut);
			mapper2->SetInputData(face);


			vtkSmartPointer<vtkLookupTable> hueLut2 =
				vtkSmartPointer<vtkLookupTable>::New();
			hueLut2->SetTableRange(0, maxDist);
			hueLut2->SetHueRange(0, 1);
			hueLut2->SetSaturationRange(1, 1);
			hueLut2->SetValueRange(0, 1);
			hueLut2->Build();


			if (scalarBar){
				
				renderer->RemoveActor2D(scalarBar);
			}

			scalarBar =	vtkSmartPointer<vtkScalarBarActor>::New();
			scalarBar->SetLookupTable(hueLut2);

			double pos[2] = { 100, 100 };
			scalarBar->SetNumberOfLabels(4);
			scalarBar->SetWidth(scalarBar->GetWidth() / 3);
			scalarBar->SetHeight(scalarBar->GetHeight()*0.75);

			//scalarBar->SetLabelTextProperty();

			double cc[3] = { 0, 0, 1 };
			vtkTextProperty* tp = scalarBar->GetLabelTextProperty();
			tp->SetFontSize(100);
			tp->SetColor(cc);
			//tp->Modified();
			scalarBar->SetLabelTextProperty(tp);
			//scalarBar->Modified();
			//scalarBar->SetLookupTable(hueLut);
			scalarBar->Modified();
			renderer->AddActor2D(scalarBar);

			
			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper2);
			actor->GetProperty()->SetPointSize(5);
			
			
			renderer->AddActor(actor);
			renderer->Modified();
			this->GetInteractor()->GetRenderWindow()->Render();
			faceNo++;
		}

		if(key == "b"){
			vector<vector<int**>> faceHists;
			for(int i=1;i<=1000;i++){

				string path;
				if (!pointCloudMode)
					path = "faces\\face" + std::to_string(i) + ".off";
				else
					path = "faces\\face" + std::to_string(i) + ".xyz";

				Mesh *currentFace = new Mesh();

				if (!pointCloudMode)
					currentFace->loadOff((char*)path.c_str());
				else
					currentFace->loadxyz((char*)path.c_str());
				
				/*
				int x = 200000  * ((i - 1) % 5);
				int y = -200000 * ((i - 1) / 5);

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
				*/

				ifstream ff("faces\\corr" + to_string(1) + ".txt");
				vector<int> corrs;
				for (int j = 0; j < 53490; j++){

					int cor;
					ff >> cor;
					corrs.push_back(cor);

				}
				ff.close();

				vector<vector<int>> correspondedFaceParts;

				for (int j = 0; j < faceParts.size(); j++){

					vector<int> temp;
					for (int p = 0; p < faceParts[j].size(); p++){

						int refIndex = faceParts[j][p];
						int inputIndex = corrs[refIndex];
						//cout << refIndex << " " << inputIndex << std::endl;
						//s->SetTuple1(refIndex, 0);
						//s->SetTuple1(inputIndex, 360);
						temp.push_back(refIndex);

					}

					correspondedFaceParts.push_back(temp);

				}

				vector<float*> partHists;
				partHists = extractContourFeatures2(currentFace, correspondedFaceParts);

				ofstream hf("faces\\contourFeatures2_" + to_string(i) + ".txt");

				hf << partHists.size() << " " << NUM_CONTOUR_POINTS << " " << 80 << endl;

				for (int h = 0; h < partHists.size(); h++){

					for (int p = 0; p < NUM_CONTOUR_POINTS; p++){

						//for (int q = 0; q < 80; q++)
							hf << partHists[h][p] << " ";

						

					}
					hf << endl;
				}
				hf.close();

				//faceHists.push_back(partHists);
				/*for (int h = 0; h < partHists.size(); h++){

					for (int p = 0; p < NUM_CONTOUR_POINTS; p++)
					{
						delete partHists[h][p];

					}

					delete partHists[h];
				}*/
				delete currentFace;
				std::cout << i << std::endl;
			}

			cout << "---------------------------------------\n";
			/*for (int j = 0; j < 100; j++){

				std::cout << matchFaceParts(faceHists[0], faceHists[j]) << std::endl;
			}*/

			std::cout << "done\n";
		}

		if (key == "c"){

			string path;
			if (!pointCloudMode)
				path = "faces\\face" + std::to_string(faceNo) + ".off";
			else
				path = "faces\\face" + std::to_string(faceNo) + ".xyz";

			Mesh *currentFace = new Mesh();

			if (!pointCloudMode)
				currentFace->loadOff((char*)path.c_str());
			else
				currentFace->loadxyz((char*)path.c_str());

			

			int x = 250000 * ((faceNo - 1) % 5);
			int y = -250000 * ((faceNo - 1) / 5);

			currentFace->shiftMesh(x, y);

			vector<int> currSegments;
			vtkPolyData *face;
			//cube->Delete();
			face = currentFace->getVTKPolyData(!pointCloudMode);

			int nVerts = face->GetPoints()->GetNumberOfPoints();
			vtkDataArray * s = face->GetPointData()->GetScalars();

			for (int i = 0; i < currentFace->verts.size(); i++)
				s->SetTuple1(i, 360);

			s->Modified();
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
			this->GetInteractor()->GetRenderWindow()->Render();


			vector<int**> faceHist, referenceHist;
			int numParts, numHists, lenHists;

			ifstream f("faces\\contourFeatures" + to_string(faceNo) + ".txt");

			f >> numParts >> numHists >> lenHists;

			for (int i = 0; i < numParts; i++){

				int** hist = new int*[numHists];
				for (int j = 0; j < numHists; j++){

					hist[j] = new int[lenHists];
					for (int k = 0; k < lenHists; k++){

						int temp;
						f >> temp;
						hist[j][k] = temp;

					}
				}

				faceHist.push_back(hist);
			}

			f.close();


			ifstream f2("faces\\contourFeatures1.txt");

			f2 >> numParts >> numHists >> lenHists;

			for (int i = 0; i < numParts; i++){

				int** hist = new int*[numHists];
				for (int j = 0; j < numHists; j++){

					hist[j] = new int[lenHists];
					for (int k = 0; k < lenHists; k++){

						int temp;
						f2 >> temp;
						hist[j][k] = temp;

					}
				}

				referenceHist.push_back(hist);
			}

			f2.close();


			for (int j = 0; j < numParts; j++){

				std::cout << matchContour(faceHist[j], referenceHist[j]) << std::endl;
			}

			for (int i = 0; i < faceHist.size(); i++){
				for (int j = 0; j < numHists; j++){

					delete faceHist[i][j];
				}
				delete faceHist[i];
			}

			for (int i = 0; i < referenceHist.size(); i++){
				for (int j = 0; j < numHists; j++){

					delete referenceHist[i][j];
				}
				delete referenceHist[i];
			}


			faceNo++;
		}

		if (key == "f"){
			
			//Mesh* aface = new Mesh();
			//aface->loadOff("faces\\face1.off");
			//aface->findNeighborhoodTriangles();
			//completedFaceParts = completeContoursByGeodesic(aface, faceParts2);
			//vector<vector<int>> innerPoints = findInnerContourVertices(aface, completedFaceParts);

			vector<vector<int>> innerPoints={{}, {}, {}, {}, {}, {}, {} };

			std::ifstream f("segmented2.xyz");

			for (int i = 0; i < initialFace->verts.size(); i++){

				float x, y, z;
				int   r, g, b;
				f >> x >> y >> z;
				f >> r >> g >> b;

				int colorId = -1;

				if (r == 0 && g == 0 && b == 0) colorId = 0;  // mouth
				if (r == 255 && g == 255 && b == 0) colorId = 1;  // chin
				if (r == 255 && g == 0 && b == 255) colorId = 2;  // nose
				if (r == 0 && g == 255 && b == 255) colorId = 3;  // left eye
				if (r == 0 && g == 255 && b == 0) colorId = 4;  // left chick
				if (r == 0 && g == 0 && b == 255) colorId = 5;  // right chick
				if (r == 255 && g == 0 && b == 0) colorId = 6;  // right eye

				segments[i] = colorId;
				isSelected[i] = colorId != -1;

			}
			f.close();
			for (int i = 0; i < initialFace->verts.size(); i++){

				int partId;
				int cid = segments[i];
				if (cid == -1)
					continue;


				if (cid == 0) partId = 0;
				if (cid == 1) partId = 6;
				if (cid == 2) partId = 1;
				if (cid == 3) partId = 2;
				if (cid == 4) partId = 4;
				if (cid == 5) partId = 5;
				if (cid == 6) partId = 3;

				innerPoints[partId].push_back(i);

				
			}


			//renderer->RemoveActor(renderer->GetActors()->GetLastActor());
			Mesh *m = new Mesh();
			vector<vector<int>> inputFaceParts;
			ifstream of("inputPoints.xyz");
			for (int fp = 0; fp < faceParts.size(); fp++){
				
				vector<int> facePart;
				for (int p = 0; p < faceParts[fp].size(); p++){

					float x, y, z;
					of >> x >> y >> z;
					float* coords = new float[3];
					coords[0] = x;
					coords[1] = y;
					coords[2] = z;
					int vid = m->addVertex(coords);
					facePart.push_back(vid);
					

				}
				inputFaceParts.push_back(facePart);
			}

			of.close();

		
			vector<int**> inputHist;
			inputHist = extractContourFeatures(m, inputFaceParts);
			/*
			ofstream hf("faces\\InputContourFeatures.txt");

			hf << partHists.size() << " " << NUM_CONTOUR_POINTS << " " << 80 << endl;

			for (int h = 0; h < partHists.size(); h++){

				for (int p = 0; p < NUM_CONTOUR_POINTS; p++){

					for (int q = 0; q < 80; q++)
						hf << partHists[h][p][q] << " ";

					hf << endl;

				}
			}
			hf.close();
			*/
			int numParts = faceParts.size();
			int* bestMatches = new int[numParts];
			float* minDists = new float[numParts];
			for (int i = 0; i < numParts; i++)
				minDists[i] = 999999999;

			for (int fn = 2; fn < 1000; fn++){

				vector<int**> faceHist;
				int numParts, numHists, lenHists;

				ifstream f("faces\\contourFeatures" + to_string(fn) + ".txt");

				f >> numParts >> numHists >> lenHists;

				for (int i = 0; i < numParts; i++){

					int** hist = new int*[numHists];
					for (int j = 0; j < numHists; j++){

						hist[j] = new int[lenHists];
						for (int k = 0; k < lenHists; k++){

							int temp;
							f >> temp;
							hist[j][k] = temp;

						}
					}

					faceHist.push_back(hist);
				}

				f.close();

				for (int part = 0; part < numParts; part++){

					float dist = matchContour(faceHist[part], inputHist[part]);

					if (dist < minDists[part]){
						minDists[part] = dist;
						bestMatches[part] = fn;

					}

				}

				for (int i = 0; i < numParts; i++){
					for (int j = 0; j < numHists; j++)
						delete faceHist[i][j];
					delete faceHist[i];
				}

			}

			Mesh* outFace = new Mesh();
			outFace->loadOff("faces\\baseface.off");

			for (int i = 0; i < faceParts.size(); i++){

				cout << bestMatches[i] << " " << minDists[i] << endl;
				
				int fid = bestMatches[i];

				Mesh* tempMesh = new Mesh();
				tempMesh->loadOff((char*)("faces\\face" + to_string(fid) + ".off").c_str());

				vector<int> vList = innerPoints[i];

				for (int v = 0; v < vList.size(); v++){

					int vid = vList[v];

					outFace->verts[vid]->coords[0] = tempMesh->verts[vid]->coords[0];
					outFace->verts[vid]->coords[1] = tempMesh->verts[vid]->coords[1];
					outFace->verts[vid]->coords[2] = tempMesh->verts[vid]->coords[2];
				}

				
				
				/*tempMesh->shiftMesh(250000,0);
				
				vtkPolyData* tempPoly = tempMesh->getVTKPolyDataSubsetByVertices(innerPoints[i], 360 );
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(tempPoly);
				mapper->SetScalarRange(0, 360);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetPointSize(5);

				renderer->AddActor(actor);
				renderer->Modified();*/
				
				

			}

			outFace->shiftMesh(250000, 0);
			
			vector<int> innerPoints2;

			for (int ii = 0; ii < innerPoints.size(); ii++){

				for (int jj = 0; jj < innerPoints[ii].size(); jj++){

					innerPoints2.push_back(innerPoints[ii][jj]);


				}

			}

			/*for (int i = 0; i < faceParts.size(); i++){
			
				vtkPolyData* tempPoly = outFace->getVTKPolyDataSubsetByVertices(innerPoints[i], 360);
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(tempPoly);
				mapper->SetScalarRange(0, 360);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetPointSize(5);

				renderer->AddActor(actor);
				renderer->Modified();

			}*/

			
			vtkPolyData* tempPoly = outFace->getVTKPolyDataSubsetByVertices(innerPoints2, 360);
			
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(tempPoly);
			mapper->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetPointSize(5);

			renderer->AddActor(actor);
			renderer->Modified();
			
			/*
			vtkPolyData* points = m->getVTKPolyData(false);
			vtkDataArray * s = points->GetPointData()->GetScalars();
			
			int vid = 0;
			for (int fp = 0; fp < inputFaceParts.size();fp++)
			for (int p = 0; p < inputFaceParts[fp].size(); p++)
				s->SetTuple1(vid++, 10 + 50 * fp);

			vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper2->SetInputData(points);
			mapper2->SetScalarRange(0, 360);

			vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
			actor2->SetMapper(mapper2);
			actor2->GetProperty()->SetPointSize(5);

			renderer->AddActor(actor2);
			renderer->Modified();
			*/
			this->GetInteractor()->GetRenderWindow()->Render();
			

			
		}

		if (key == "v"){

			
			
			Mesh* testFace;

			for (int tf = 1; tf < 2; tf++){
				const clock_t begin_time = clock();
				std::cout << tf << std::endl;
				testFace = new Mesh();

				int testFaceId = tf; // tf;// rand() % 100;
				testFace->loadOff((char*)("faces\\testface" + to_string(testFaceId) + ".off").c_str());
			
				createInputFile(testFace, faceParts);

				
				vector<vector<int>> innerPoints = { {}, {}, {}, {}, {}, {}, {} };
				int colors[53490];
				std::ifstream f("faces\\segmented.xyz");

				for (int i = 0; i < 53490; i++){

					float x, y, z;
					int   r, g, b;
					f >> x >> y >> z;
					f >> r >> g >> b;

					int colorId = -1;

					if (r == 0 && g == 0 && b == 0) colorId = 0;  // mouth
					if (r == 255 && g == 255 && b == 0) colorId = 1;  // chin
					if (r == 255 && g == 0 && b == 255) colorId = 2;  // nose
					if (r == 0 && g == 255 && b == 255) colorId = 3;  // left eye
					if (r == 0 && g == 255 && b == 0) colorId = 4;  // left chick
					if (r == 0 && g == 0 && b == 255) colorId = 5;  // right chick
					if (r == 255 && g == 0 && b == 0) colorId = 6;  // right eye

					colors[i] = colorId;
					

				}
				f.close();

				for (int i = 0; i < 53490; i++){

					int partId;
					int cid = colors[i];
					if (cid == -1)
						continue;


					if (cid == 0) partId = 0;
					if (cid == 1) partId = 6;
					if (cid == 2) partId = 1;
					if (cid == 3) partId = 2;
					if (cid == 4) partId = 4;
					if (cid == 5) partId = 5;
					if (cid == 6) partId = 3;

					innerPoints[partId].push_back(i);


				}

				
				//renderer->RemoveActor(renderer->GetActors()->GetLastActor());
				Mesh *m = new Mesh();
				vector<vector<int>> inputFaceParts;
				ifstream of("inputPoints.xyz");
				for (int fp = 0; fp < faceParts.size(); fp++){

					vector<int> facePart;
					for (int p = 0; p < faceParts[fp].size(); p++){

						float x, y, z;
						of >> x >> y >> z;
						float* coords = new float[3];
						coords[0] = x;
						coords[1] = y;
						coords[2] = z;
						int vid = m->addVertex(coords);
						facePart.push_back(vid);


					}
					inputFaceParts.push_back(facePart);
				}

				of.close();


				vector<float*> inputHist;
				inputHist = extractContourFeatures2(m, inputFaceParts);


				ofstream hf("faces\\contourFeatures_input" + to_string(testFaceId) + ".txt");

				hf << inputHist.size() << " " << NUM_CONTOUR_POINTS << " " << 80 << endl;

				for (int h = 0; h < inputHist.size(); h++){

					for (int p = 0; p < NUM_CONTOUR_POINTS; p++){

						//for (int q = 0; q < 80; q++)
							hf << inputHist[h][p] << " ";

						

					}
					hf << endl;
				}
				hf.close();

				int numParts = faceParts.size();
				int* bestMatches = new int[numParts];
				float* minDists = new float[numParts];
				for (int i = 0; i < numParts; i++)
					minDists[i] = 999999999;
				int numHists;
				

				for (int fn = 1; fn < 1000; fn++){

					//if (!(fn == 694 || fn == 316 || fn == 28 || fn == 993 || fn == 688 || fn == 838 || fn == 389))
						//continue;

					//if (fn == 821 || fn == 79 || fn==141 || fn==35 || fn==240)
						//continue;

					vector<float*> faceHist;
					int numParts, numHists, lenHists;

					ifstream f("faces\\contourFeatures2_" + to_string(fn) + ".txt");

					f >> numParts >> numHists >> lenHists;

					for (int i = 0; i < numParts; i++){

						float* hist = new float[numHists];
						for (int j = 0; j < numHists; j++){

							/*hist[j] = new int[lenHists];
							for (int k = 0; k < lenHists; k++){

								int temp;
								f >> temp;
								hist[j][k] = temp;

							}*/

							float temp;
							f >> temp;
							hist[j] = temp;
						}

						faceHist.push_back(hist);
					}

					f.close();

					for (int part = 0; part < numParts; part++){

						float dist = matchFaceParts(faceHist[part], inputHist[part]);

						if (dist < minDists[part]){
							minDists[part] = dist;
							bestMatches[part] = fn;

						}

					}

					/*for (int i = 0; i < numParts; i++){
						for (int j = 0; j < numHists; j++)
							delete faceHist[i][j];
						delete faceHist[i];
					}*/

				}

				
				Mesh* outFace = new Mesh();
				Mesh* outFace2 = new Mesh();
				outFace->loadOff("faces\\baseface.off");
				outFace2->loadOff("faces\\baseface.off");

				vector<Mesh*> bestMatchFaces;

				for (int i = 0; i < faceParts.size(); i++){

					cout << bestMatches[i] << " " << minDists[i] << endl;

					int fid = bestMatches[i];

					Mesh* tempMesh = new Mesh();
					tempMesh->loadOff((char*)("faces\\face" + to_string(fid) + ".off").c_str());

					bestMatchFaces.push_back(tempMesh);

					

				}


				
				for (int i = 0; i < faceParts.size(); i++){


					float xd = 0;
					float yd = 0;
					float zd = 0;

					
					vector<int> vList = innerPoints[i];

					for (int v = 0; v < vList.size(); v++){

						int vid = vList[v];

						//bestMatchFaces[i]->verts[vid]->coords[0] += 0.999*xd;
						//bestMatchFaces[i]->verts[vid]->coords[1] += 0.999*yd;
						//bestMatchFaces[i]->verts[vid]->coords[2] += 0.999*zd;

						outFace->verts[vid]->coords[0] = bestMatchFaces[i]->verts[vid]->coords[0];
						outFace->verts[vid]->coords[1] = bestMatchFaces[i]->verts[vid]->coords[1];
						outFace->verts[vid]->coords[2] = bestMatchFaces[i]->verts[vid]->coords[2];
					}

				}

				
				

				vector<int> innerPoints2;
				for (int ii = 0; ii < innerPoints.size(); ii++)
				for (int jj = 0; jj < innerPoints[ii].size(); jj++)
					innerPoints2.push_back(innerPoints[ii][jj]);

				vector<int> vertList;
				for (int ii = 0; ii < faceParts.size(); ii++)
				for (int jj = 0; jj < faceParts[ii].size(); jj++)
					vertList.push_back(faceParts[ii][jj]);

				outFace2->initialize();
				outFace2->deform(outFace, 10000, 5, innerPoints2);
				outFace2->shiftMesh(250000,0);
				outFace2->write((char*)("faces\\testface" + to_string(testFaceId) + "_out.off").c_str());


				//outFace->shiftMesh(250000,0);
				//outFace->write((char*)("faces\\testface" + to_string(testFaceId) + "_out2.off").c_str());

				

				//Mesh* partsOnly = outFace->getSubsetByVertices(innerPoints2);
				//partsOnly->write((char*)("faces\\testface" + to_string(testFaceId) + "_out2.off").c_str());
				std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << "  secs\n\n";

				/*
				outFace->shiftMesh(250000, 250000*(tf-1));
				vector<int> innerPoints2;

				for (int ii = 0; ii < innerPoints.size(); ii++)
					for (int jj = 0; jj < innerPoints[ii].size(); jj++)
						innerPoints2.push_back(innerPoints[ii][jj]);
					

				vtkPolyData* outPoly = outFace->getVTKPolyDataSubsetByVertices(innerPoints2, 360);

				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(outPoly);
				mapper->SetScalarRange(0, 360);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetPointSize(5);

				renderer->AddActor(actor);

				testFace->shiftMesh(0, 250000 * (tf - 1));
				vtkPolyData* inPoly = testFace->getVTKPolyDataSubsetByVertices(innerPoints2, 360);
				
				vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper2->SetInputData(inPoly);
				mapper2->SetScalarRange(0, 360);

				vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
				actor2->SetMapper(mapper2);
				actor2->GetProperty()->SetPointSize(5);

				renderer->AddActor(actor2);

				renderer->Modified();
				*/
				/*for (int i = 0; i < numParts; i++){
					for (int j = 0; j < numHists; j++)
						delete inputHist[i][j];
					delete inputHist[i];
				}*/
				for (int i = 0; i < bestMatchFaces.size(); i++)
					delete bestMatchFaces[i];
				delete testFace;
				delete outFace;
				delete outFace2;
				delete m;
				delete minDists;
				//delete partsOnly;
				delete bestMatches;
				this->GetInteractor()->GetRenderWindow()->Render();


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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////     TPS-RPM   ////////////////////////////////////////////////////////////
#define INF							90000000000000.0f
#define EXP							2.718281828f
#include "MtrxOps.h"
/*float distanceBetween(float* v1, float* v2)
{
	//In fact, the keyword inline is not necessary. If the function is defined with its body directly and the function
	//has a smaller block of code, it will be automatically treated as inline by the compiler
	return (float)sqrt((v2[0] - v1[0])*(v2[0] - v1[0]) + (v2[1] - v1[1])*(v2[1] - v1[1]) + (v2[2] - v1[2])*(v2[2] - v1[2]));
}*/
inline float distanceBetween2(float* v1, float* v2)
{
	
	return sqrt((v2[0] - v1[0])*(v2[0] - v1[0]) + (v2[1] - v1[1])*(v2[1] - v1[1]) + (v2[2] - v1[2])*(v2[2] - v1[2]));
}

float computeMSEe(vector<Vertex*> pntSet1, vector<Vertex*> pntSet2)
{
	//computes mean squared error (mse) b/w 2 pnts sets of size n1 & n2, that is sum of squared distances of closest pairs
	int n1 = pntSet1.size();
	int n2 = pntSet2.size();
	//before computing mse, i need closest-pnt correspondence b/w mesh1 baseVerts & mesh2 baseVerts
	int* closestPntInSet1 = new int[n2]; //closest pnt in pntSet1 for given pntSet2 vertex (therefore n2 elements)
	float dist;
	for (int bv2 = 0; bv2 < n2; bv2++)
	{
		closestPntInSet1[bv2] = bv2;
		float minDist = INF;
		for (int bv1 = 0; bv1 < n1; bv1++)
		{
			dist = distanceBetween2(pntSet2[bv2]->coords, pntSet1[bv1]->coords); //since just comparison, not use sqrt to save time
			if (dist < minDist)
			{
				closestPntInSet1[bv2] = bv1; //closest pnt in pntSet1 for pntSet2 vertex bv2 is bv1
				minDist = dist;
			}
		}
	}
	float mse = 0.0f; //now i can compute this
	for (int bv2 = 0; bv2 < n2; bv2++)
	{
		//squared distance b/w (target) mesh2.bv and its correspondence in mesh1 for this iteration (closestBaseInM1)
		dist = pow(pntSet2[bv2]->coords[0] - pntSet1[closestPntInSet1[bv2]]->coords[0], 2.0f) +
			pow(pntSet2[bv2]->coords[1] - pntSet1[closestPntInSet1[bv2]]->coords[1], 2.0f) +
			pow(pntSet2[bv2]->coords[2] - pntSet1[closestPntInSet1[bv2]]->coords[2], 2.0f);
		//		if (! rootBasedSymSoln)
		mse += sqrt(dist); //no sqrt() above, hence dist = squaredDist indeed		
		//		else
		//			mse += pow(gd2r2[bv2] - gd2r1[ closestPntInSet1[bv2] ], 2.0f); //replace w/ += dist above for symmetry-problem solution
	}
	mse /= n2;
	delete[] closestPntInSet1;
	return mse;
}

float findT0(vector<Vertex*> pntSet1, vector<Vertex*> pntSet2)
{

	int n1 = pntSet1.size();
	int n2 = pntSet2.size();

	float maxDist2 = -INF; //this T0 best so far
	for (int bv = 0; bv < n1; bv++)
	{
		for (int bv2 = 0; bv2 < n2; bv2++)
		{
			float* v1 = pntSet1[bv]->coords;
			float* v2 = pntSet2[bv2]->coords; //pntSet2[bv]->coords;
			//wrongly v2 = pntSet2[bv]->coords was in use while i was running paper tests; now, i corrected it!!!!!!

			float dist2 = (v2[0] - v1[0])*(v2[0] - v1[0]) + (v2[1] - v1[1])*(v2[1] - v1[1]) + (v2[2] - v1[2])*(v2[2] - v1[2]);
			dist2 = sqrt(dist2);
			if (dist2 > maxDist2)
				maxDist2 = dist2;
		}
	}
	return maxDist2; //pairFrom2SetsWithSameIndex*/

}


void mtrxMult(float** result, int m1Row, int m1Col, int m2Row, int m2Col, float** m1, float** m2)
{
	//result = m1 * m2

	//init result to 0
	for (int r = 0; r < m1Row; r++)
	for (int c = 0; c < m2Col; c++)
		result[r][c] = 0.0f;

	/*pseudocode:
	for each row i in matrix1 do:
	for each column j in matrix 2 do:
	for each column k in matrix1 do:
	increase  matrix3[i,j] by matrix1[i,k]*matrix2[k,j] */
	for (int m1r = 0; m1r < m1Row; m1r++)
	for (int m2c = 0; m2c < m2Col; m2c++)
	for (int m1c = 0; m1c < m1Col; m1c++)
		result[m1r][m2c] += m1[m1r][m1c] * m2[m1c][m2c]; //call-by-ref to result
}

float minMSE = INF;
float nonRigidAlignment(vector<Vertex*> &pntSet1, vector<Vertex*> &pntSet2)
{
	int n1 = pntSet1.size();
	int n2 = pntSet2.size();

	//employs thin-plate spline, robust pnt mathing (TPS-RPM) non-rigid transformation to transform pntSet1 to target mesh2

	float initialMSE = computeMSEe(pntSet1, pntSet2); //dellaterrrrrrrrrrr (similarly, int k above is required for fprint only)

	cout << "non-rigid TPS-RPM to align 2 patches....\n";

	int N = n2, K = n1; //N & K not necessarily equal for this function

	//center of masses (or means) that will act as outlier cluster centers for last row and last col of m mtrx
	float m1[3];
	float m2[3];
	float sum[3] = { 0.0f, 0.0f, 0.0f };
	for (int bv = 0; bv < N; bv++)
	{
		sum[0] += pntSet2[bv]->coords[0];
		sum[1] += pntSet2[bv]->coords[1];
		sum[2] += pntSet2[bv]->coords[2];
	}
		
	m2[0] = sum[0] / (float)N;
	m2[1] = sum[1] / (float)N;
	m2[2] = sum[2] / (float)N;

	float sum2[3] = { 0.0f, 0.0f, 0.0f };
	for (int bv = 0; bv < K; bv++)
	{
		sum2[0] += pntSet1[bv]->coords[0];
		sum2[1] += pntSet1[bv]->coords[1];
		sum2[2] += pntSet1[bv]->coords[2];
	}
	
	m1[0] = sum2[0] / (float)N;
	m1[1] = sum2[1] / (float)N;
	m1[2] = sum2[2] / (float)N;

	//no need for initial its and/or principalAxAlignment in alignWith() 'cos TPS-RPM handles them
	//actually it's better to start w/ aligned center-of-masses; doing that in caller

	//////////// TPS-RPM algo ////////////
	//init control parameters
	int perT_maxit = 500;
	float anneal_rate = 0.93, //0.97f
		lamda1_init = 100000.000f, lamda2_init = 0.01f;
	//other params
	float T0 = findT0(pntSet1, pntSet2); //computeCorrespondenceError() = 3.3 x 10^7 (largestX), = 3.2 x 10^7 (pairFrom2SetsWithSameIndex)
		//initialMSE * 10.0f, //correspErr = 2.8 x 10^7
		//initialMSE * 7.0f, //correspErr = 2.8 x 10^7
		//initialMSE * 1.0f, //correspErr = 3.5 x 10^7
		//2500.0f, //correspErr = 2.3 x 10^7
		//[correspErr of input\animation\5-->212 = 1.8 x 10^11 but this data is huge, e.g. bbox.(w,h,d) = (5225.22, 7940.11, 18179.8), whereas jumping.bbox = (601.231, 1611.88, 622.9)
		//best is pairFrom2SetsWithSameIndex; others' correspErr lower but no worry 'cos they produce >badMSE nonRigidAlignments and hence cancelled most of the time
		//note also that pairFrom2SetsWithSameIndex used in findT0() corresponds to ground-truth correspondences (so, you may use L2-closest correspondence for pairing)
	float Tfinal = 0.1f;
	float T = T0; //start temparature at Tinit and gradually decrease it till Tfinal in determenistic annealing loop below
	//Tfinal = is chosen to be equal to the average of the squared distance between the nearest neighbors within the set of points which are being deformed; interpretation
	//is that at Tfinal, the Gaussian clusters for all the points will then barely overlap with one another (currently just use nIters or mse to quit)

	//K+1 x N+1 fuzzy correspondence matrix (+1 for outliers row & col)
	float** m = new float*[K + 1];
	for (int a = 0; a < K + 1; a++)
		m[a] = new float[N + 1];
	//auxiliary arrays for iterated row/col normalization below
	float* denomR = new float[K],
		*denomC = new float[N], threshold = 0.005f;

	//Y matrix that will be updated by m and stable mesh2.coords below and then guide mesh1.coords
	float** Y = new float*[K]; //Kx4, e.g. row7 = [1, 40, 45, 48] (1st col always 1)
	for (int a = 0; a < K; a++)
	{
		Y[a] = new float[4];
		Y[a][0] = 1.0f; //will never change again; so, do it here for once
	}
	float** V = new float*[K], //Kx4, e.g. row7 = [1, 40, 45, 48] (1st col always 1) V = initial mesh1.coords and used only for qr & F[][]
		** Vc = new float*[K]; //copy of V, which will not be modified by mtrxInv() and therefore usable inside the annealing loop
	for (int a = 0; a < K; a++)
	{
		V[a] = new float[4];
		Vc[a] = new float[3]; //no need to keep [0] = 1 for Vc
		V[a][0] = 1.0f;
		for (int c = 1; c < 4; c++)
		{
			V[a][c] = pntSet1[a]->coords[c - 1];
			Vc[a][c - 1] = V[a][c]; //so, Vc[7] is the coords of m1.baseVerts7
		}
	}

	//KxK F (phi, or kernel) mtrx, constant through annealing
	float** F = new float*[K];
	for (int a = 0; a < K; a++)
	{
		F[a] = new float[K]; //F[a][c] = dist between V[c] & V[a]
		for (int c = 0; c < K; c++)
			F[a][c] = -sqrt((V[c][1] - V[a][1])*(V[c][1] - V[a][1]) + (V[c][2] - V[a][2])*(V[c][2] - V[a][2]) + (V[c][3] - V[a][3])*(V[c][3] - V[a][3])); //[*][0] = 1, so no effect
	}

	//QR-decompose V to fetch Q1, Q2, R (same for all iterations) to be used in computation of transformation f = (d, w) below
	MtrxOps* mtrx = new MtrxOps();
	float* tau = mtrx->rmatrixqr(V, K, 4); //V is also updated via call-by-ref
	//get KxK Q
	float** Q = mtrx->rmatrixqrunpackq(V, K, 4, tau, K); //last K: i want first K, i.e. all, columns of Q
	//orthonormal = orthogonal (perpendicular) + normal (magnitude 1); so, turn orthogonal vectors of Q into orthonormal basis
	for (int ac = 0; ac < K; ac++) //i guess qr.rmatrixqrunpackq returns Q already normalized but be safe
	{
		float len = 0.0f;
		for (int ar = 0; ar < K; ar++) //find normalizer
			len += Q[ar][ac] * Q[ar][ac];
		len = sqrt(len);
		for (int ar = 0; ar < K; ar++) //normalize
			Q[ar][ac] /= len;
	}
	//get Kx4 RR
	float** RR = mtrx->rmatrixqrunpackr(V, K, 4);

	//get Q1 & Q2 from KxK Q (also keep their transposes)
	float** Q1 = new float*[K]; //Kx4, where 4 = D+1, where D = 3, i.e. 3D pnts
	for (int a = 0; a < K; a++)
	{
		Q1[a] = new float[4];
		for (int c = 0; c < 4; c++)
			Q1[a][c] = Q[a][c];
	}
	float** Q2 = new float*[K]; //Kx(K-4), where 4 = D+1, where D = 3, i.e. 3D pnts
	for (int a = 0; a < K; a++)
	{
		Q2[a] = new float[K - 4];
		for (int c = 0; c < K - 4; c++)
			Q2[a][c] = Q[a][c + 4];
	}
	float** Q1t = new float*[4]; //4xK
	for (int a = 0; a < 4; a++)
	{
		Q1t[a] = new float[K];
		for (int c = 0; c < K; c++)
			Q1t[a][c] = Q1[c][a];
	}
	float** Q2t = new float*[K - 4]; //(K-4)xK
	for (int a = 0; a < K - 4; a++)
	{
		Q2t[a] = new float[K];
		for (int c = 0; c < K; c++)
			Q2t[a][c] = Q2[c][a];
	}

	//get 4x4 R from Kx4 RR (by omitting all-zero part below row 4)
	float** R = new float*[4]; //4x4, where 4 = D+1, where D = 3, i.e. 3D pnts
	for (int a = 0; a < 4; a++)
	{
		R[a] = new float[4];
		for (int c = 0; c < 4; c++)
			R[a][c] = RR[a][c];
	}
	float** Rt = new float*[4], //4x4 (transpose of R)
		** RtR = new float*[4]; //4x4 (Rt * R)
	for (int a = 0; a < 4; a++)
	{
		RtR[a] = new float[4];
		Rt[a] = new float[4];
		for (int c = 0; c < 4; c++)
			Rt[a][c] = R[c][a];
	}
	mtrxMult(RtR, 4, 4, 4, 4, Rt, R);

	//auxiliary matrices that holds temporary mtrxMult results
	//constants ones for w computation
	float** tmpK4K = new float*[K - 4]; //K-4xK mtrx
	for (int i = 0; i < K - 4; i++)
		tmpK4K[i] = new float[K];
	float** tmpK4K4 = new float*[K - 4], ** tmpK4K4a = new float*[K - 4]; //K-4xK-4 mtrx
	for (int i = 0; i < K - 4; i++)
	{
		tmpK4K4[i] = new float[K - 4];
		tmpK4K4a[i] = new float[K - 4]; //lambda1 added version of tmpK4K, hence the name tmpK4K4a
	}
	float** tmpKK4 = new float*[K]; //KxK-4 mtrx
	for (int i = 0; i < K; i++)
		tmpKK4[i] = new float[K - 4];
	float** gamma = new float*[K - 4]; //K-4x4 mtrx
	for (int i = 0; i < K - 4; i++)
		gamma[i] = new float[4];

	mtrxMult(tmpK4K, K - 4, K, K, K, Q2t, F);
	mtrxMult(tmpK4K4, K - 4, K, K, K - 4, tmpK4K, Q2);

	//constants ones for d computation (to save time in the annealing loop)
	float** RtQ1t = new float*[4]; //4xK
	for (int a = 0; a < 4; a++)
		RtQ1t[a] = new float[K];
	mtrxMult(RtQ1t, 4, 4, 4, K, Rt, Q1t); //Rt * Q1t
	float** FQ2 = new float*[K]; //KxK-4
	for (int a = 0; a < K; a++)
		FQ2[a] = new float[K - 4];
	mtrxMult(FQ2, K, K, K, K - 4, F, Q2); //F * Q2
	float** tmpK4 = new float*[K]; //Kx4
	for (int a = 0; a < K; a++)
		tmpK4[a] = new float[4]; //will hold FQ2*gamma in the annealing loop below
	float** tmp44a = new float*[4], ** tmp44b = new float*[4]; //4x4
	for (int a = 0; a < 4; a++)
	{
		tmp44a[a] = new float[4];
		tmp44b[a] = new float[4]; //d = tmp44a * tmp44b
	}

	//desired transformation f = (w, d)
	float** w = new float*[K]; //Kx4 mtrx
	for (int i = 0; i < K; i++)
		w[i] = new float[4];
	float** d = new float*[4]; //4x4 mtrx
	for (int i = 0; i < 4; i++)
		d[i] = new float[4];

	float lambda1, lambda2, lambda1_init = 1.0f, lambda2_init = lambda1_init * 0.01f;
	//////////// outlier row/col handling ///////////////
	float moutlier = (1.0f / sqrt(T0)) * pow(EXP, -1), *sy = new float[N];
	for (int i = 0; i < N; i++) //all N columns (N+1'th column not updated (< N) 'cos never used anyway (lower-right corner of mtrx))
		m[K][i] = moutlier;
	for (int a = 0; a < K; a++) //all K rows (K+1'th row not updated 'cos never used anyway)
		m[a][N] = moutlier; //these outlier values will never change again (MatLab code)
	//////////// outlier row/col handling ends ///////////////
	int nIters = 0;
	float mse = INF, prevMSE; //loop terminator, minimum squared err between currently aligned pnt sets
	float** prevCoords1 = new float*[K]; //previous coords of deforming pntSet1 in case i need to roll-back to them after loop below
	for (int k = 0; k < K; k++)
		prevCoords1[k] = new float[3];
	while (true) //T >= Tfinal) //deterministic annealing iterations as long as T is hot/high enough
	{
		for (int iterT = 0; iterT < perT_maxit; iterT++)
		{
			//alternating update in 2 steps

			//step1: update the correspondence (given transforming mesh1.coords and stable mesh2.coords, update m)
			//the inner KxN part is updated
			for (int a = 0; a < K; a++)
			for (int i = 0; i < N; i++)
			{
				float top = -(pow(pntSet2[i]->coords[0] - pntSet1[a]->coords[0], 2.0f) +
					pow(pntSet2[i]->coords[1] - pntSet1[a]->coords[1], 2.0f) +
					pow(pntSet2[i]->coords[2] - pntSet1[a]->coords[2], 2.0f)) / (2 * T),
					left = (1.0f / sqrt(T));
				m[a][i] = left * pow(EXP, top);
			}
			//replacing outlier handling and normalization previously-above w/ the stuff from paper author's MatLab code
			//% normalize accross the outliers as well:
			for (int i = 0; i < N; i++)
				sy[i] = 0.0f; //init sy
			for (int i = 0; i < N; i++)
			for (int r = 0; r <= K; r++) //include outlier row as well
				sy[i] += m[r][i]; //fill sy, e.g. sy[7] = sum of elements (including outliers) in column7 of m
			//the inner KxN part is normalized by sy
			for (int a = 0; a < K; a++)
			for (int i = 0; i < N; i++)
				m[a][i] /= sy[i];

			//step2: update the transformation (f = (d, w))
			//update all K entries of Y
			for (int a = 0; a < K; a++)
			{
				float est[3] = { 0.0f, 0.0f, 0.0f };
				for (int i = 0; i < N; i++)
				{
					est[0] += m[a][i] * pntSet2[i]->coords[0];
					est[1] += m[a][i] * pntSet2[i]->coords[1];
					est[2] += m[a][i] * pntSet2[i]->coords[2];
				}
					
				for (int c = 1; c < 4; c++)
					Y[a][c] = est[c - 1]; //Y[*][0] = 1 set above and stays as is
				//so, Y[a] can be regarded as a mixture of mesh2 coords that corresponds to mesh1.baseVerts[a]
			}

			//new annealing parameters that'll affect matrices below
			lambda1 = lamda1_init * K * T;	lambda2 = lamda2_init * K * T; //these (big tmpK4K4a diagonals)
			//lambda1 = lamda1_init * T;		lambda2 = lamda2_init * T; //or these? (smaller diagonals for tmpK4K4a)

			//compute w
			for (int a = 0; a < K - 4; a++)
			for (int c = 0; c < K - 4; c++)
				tmpK4K4a[a][c] = tmpK4K4[a][c] + (a == c ? lambda1 : 0.0f); //lambda1 added version of tmpK4K4
			/*FILE* fPtr2;
			if (iterT == perT_maxit-1){
			fPtr2 = fopen("invinput.dat", "w");
			for (int a = 0; a < K-4; a++) for (int i = 0; i < K-4; i++) fprintf(fPtr2, "%f%s", tmpK4K4a[a][i], (i == K-5 ? ";\n" : ", ")); fclose(fPtr2);
			}*/
			mtrx->mtrxInv(tmpK4K4a, K - 4); //tmpK4K4a = inv(tmpK4K4a) via call-by-ref
			/*if (iterT == perT_maxit-1){
			fPtr2 = fopen("invoutput.dat", "w");
			for (int a = 0; a < K-4; a++) for (int i = 0; i < K-4; i++) fprintf(fPtr2, "%f%s", tmpK4K4a[a][i], (i == K-5 ? ";\n" : ", "));fclose(fPtr2);
			}*/

			mtrxMult(tmpK4K, K - 4, K - 4, K - 4, K, tmpK4K4a, Q2t);
			mtrxMult(gamma, K - 4, K, K, 4, tmpK4K, Y); //K-4x4 gamma
			mtrxMult(w, K, K - 4, K - 4, 4, Q2, gamma); //w = Q2 * gamma is ready now

			//compute d
			for (int a = 0; a < 4; a++)
			for (int c = 0; c < 4; c++)
				tmp44a[a][c] = RtR[a][c] + (a == c ? lambda2 : 0.0f); //lambda2 added version of RtR
			mtrx->mtrxInv(tmp44a, 4); //tmp44a = inv(tmp44a) via call-by-ref
			mtrxMult(tmpK4, K, K - 4, K - 4, 4, FQ2, gamma); //FQ2 * gamma
			for (int a = 0; a < K; a++)
			for (int c = 0; c < 4; c++)
				tmpK4[a][c] = Y[a][c] - tmpK4[a][c]; //tmpK4 = Y - tmpK4 effect
			mtrxMult(tmp44b, 4, K, K, 4, RtQ1t, tmpK4);
			for (int a = 0; a < 4; a++)
			for (int c = 0; c < 4; c++)
				tmp44b[a][c] -= RtR[a][c]; //tmp44b = tmp44b - RtR effect
			mtrxMult(d, 4, 4, 4, 4, tmp44a, tmp44b);
			for (int a = 0; a < 4; a++)
			for (int c = 0; c < 4; c++)
			if (a == c)
				d[a][c] += 1.0f; //1 added version of d (off-diagonal entries are already set by 'mtrxMult(d, ..' above)

			//remember previous values before starting to change them to current values (in case i roll-back after outer while-loop)
			prevMSE = mse;

			//apply transformation to mesh1.coords (which will affect m -> Y -> gamma -> d -> w (effect-chain) in next iteration)
			mtrxMult(tmpK4, K, K, K, 4, F, w); //F*w is added to affine-transformed coords below, i.e. after d-transformation
			for (int bv = 0; bv < K; bv++)
			{
				prevCoords1[bv][0] = pntSet1[bv]->coords[0]; //remember previous values before starting to change them to current values (in case i roll-back after outer while-loop)
				prevCoords1[bv][1] = pntSet1[bv]->coords[1];
				prevCoords1[bv][2] = pntSet1[bv]->coords[2];

				//apply affine transformation d
				float n1 = d[0][0] + Vc[bv][0] * d[1][0] + Vc[bv][1] * d[2][0] + Vc[bv][2] * d[3][0]; //Vc[*][0] = x, [1] = y, [2] = z
				pntSet1[bv]->coords[0] = d[0][1] + Vc[bv][0] * d[1][1] + Vc[bv][1] * d[2][1] + Vc[bv][2] * d[3][1];
				pntSet1[bv]->coords[1] = d[0][2] + Vc[bv][0] * d[1][2] + Vc[bv][1] * d[2][2] + Vc[bv][2] * d[3][2];
				pntSet1[bv]->coords[2] = d[0][3] + Vc[bv][0] * d[1][3] + Vc[bv][1] * d[2][3] + Vc[bv][2] * d[3][3];
				//apply warping coefficient matrix w (actually tmpK4 = F*w is applied)
				float n2 = tmpK4[bv][0];
				pntSet1[bv]->coords[0] += tmpK4[bv][1];
				pntSet1[bv]->coords[1] += tmpK4[bv][2];
				pntSet1[bv]->coords[2] += tmpK4[bv][3];
				/*				float wn = n1+n2;
				if (wn != 1) //never entering this if-block, which is good
				{
				cout << wn << " != 1; do homogenous normalization?\n"; //by dividing coords by wn?
				pntSet1[bv]->coords /= wn; //never coming here :)
				}*/
			} //end of bv

			/*//new center of masses to be used as outlier cluster centers below (m2 not change)
			sum = SbVec3f(0.0f, 0.0f, 0.0f);
			for (int bv = 0; bv < K; bv++)
			sum += pntSet1[bv]->coords;
			m1 = sum / (float) K; //???????? i guess should not be here*/
		} //end of perT_maxit

		//gradually decrease temperature; when it is cold enough, i.e. T ~ 0, fuzzy m will be desired binary correspondence mtrx
		T = T * anneal_rate;
		

		//compute new mse value as an iteration condition
		mse = computeMSEe(pntSet1, pntSet2);

		cout << "T = " << T << ", lambda1 = " << lambda1 << ", lambda2 = " << lambda2 <<", mse= "<< mse << "\n"; //"\n";

		if (nIters++ > 20000 || (nIters > 5000 && mse > prevMSE)) //|| T < 167 is good for k = 101 of frontiers[k])
			break;				//nIters > 10 is to skip a possible mse > prevMSE case in the beginning phase
		//for k = 127, 128, 129, it was not enough to use 10; maby 20?
	} //end of while (true)

	if (nIters > 2 && mse > prevMSE) //then need to use prevCoords to roll-back to that prevMSE
	{
		cout << "rolling back to coords of mse = " << prevMSE << endl;
		for (int bv = 0; bv < K; bv++)
		{
			pntSet1[bv]->coords[0] = prevCoords1[bv][0]; //roll-back
			pntSet1[bv]->coords[1] = prevCoords1[bv][1]; //roll-back
			pntSet1[bv]->coords[2] = prevCoords1[bv][2]; //roll-back
		}
			
	}

	//FILE* fPtr = fopen("t0.dat", "a"); //see which T0 gave the best/min mse
	float mseUsed = (mse > prevMSE ? prevMSE : mse);
	//fprintf(fPtr, "%d\t%f\t\t%f\t\t%f\t\t%s\n", k, T0, initialMSE, mseUsed, (mseUsed < minMSE ? "min" : ""));
	if (mseUsed < minMSE) minMSE = mseUsed;
	//fclose(fPtr);

	//memo re-capture
	delete[] denomR;
	delete[] denomC;
	delete[] sy;
	delete[] prevCoords1;
	delete mtrx;
	for (int i = 0; i < K; i++) delete[] Y[i]; delete[] Y;
	for (int a = 0; a < K; a++) delete[] V[a]; delete[] V;
	for (int a = 0; a < K; a++) delete[] Vc[a]; delete[] Vc;
	for (int a = 0; a < K; a++) delete[] Q[a]; delete[] Q;
	for (int a = 0; a < K; a++) delete[] RR[a]; delete[] RR;
	for (int a = 0; a < K; a++) delete[] Q1[a]; delete[] Q1;
	for (int a = 0; a < K; a++) delete[] Q2[a]; delete[] Q2;
	for (int a = 0; a < 4; a++) delete[] R[a]; delete[] R;
	for (int a = 0; a < K; a++) delete[] F[a]; delete[] F;
	for (int a = 0; a < K - 4; a++) delete[] tmpK4K[a]; delete[] tmpK4K;
	for (int a = 0; a < K - 4; a++) delete[] tmpK4K4[a]; delete[] tmpK4K4;
	for (int a = 0; a < K - 4; a++) delete[] tmpK4K4a[a]; delete[] tmpK4K4a;
	for (int a = 0; a < K; a++) delete[] tmpKK4[a]; delete[] tmpKK4;
	for (int a = 0; a < 4; a++) delete[] tmp44a[a]; delete[] tmp44a;
	for (int a = 0; a < 4; a++) delete[] tmp44b[a]; delete[] tmp44b;
	for (int a = 0; a < K - 4; a++)	delete[] gamma[a]; delete[] gamma;
	for (int a = 0; a < 4; a++) delete[] Q1t[a]; delete[] Q1t;
	for (int a = 0; a < K - 4; a++) delete[] Q2t[a]; delete[] Q2t;
	for (int a = 0; a < K; a++) delete[] tmpK4[a]; delete tmpK4;
	for (int a = 0; a < 4; a++) delete[] Rt[a]; delete[] Rt;
	for (int a = 0; a < K; a++) delete[] FQ2[a]; delete[] FQ2;
	for (int a = 0; a < 4; a++) delete[] RtQ1t[a]; delete[] RtQ1t;
	for (int a = 0; a < K; a++) delete[] w[a]; delete[] w;
	for (int a = 0; a < 4; a++) delete[] d[a]; delete[] d;
	//////////// TPS-RPM algo ends ////////////

	//cout << "Done! after " << nIters-1 << " TPS-RPM iterations; final mse: " << mse << "\n\n";
	cout << "Done!\n\n";
	return mse;


}







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////CONTOUR FEATURES///////////////////////////////////////////////////////////

vector<vector<int>> completeContoursByGeodesic(Mesh* cF, vector<vector<int>> inFaceParts){

	vector<vector<int>> outCompletedFaceParts;

	//Mesh* cF = new Mesh();
	//cF->loadOff("faces\\face2.off");

	float maxGeoDist = 0;
	Dijsktra::DijsktraSP dsp;
	dsp.setMesh(cF);


	for (int fp = 0; fp < inFaceParts.size(); fp++){

		vector<int> completedPart;
		for (int i = 0; i < inFaceParts[fp].size() - 1; i++){


			int curr = inFaceParts[fp][i];
			int next = inFaceParts[fp][i + 1];

			std::vector<Dijsktra::Node*> dists = dsp.run(curr);

			Dijsktra::Node* result = dists[next]->pred;


			//std::cout << "from " << curr << " to " << next << ":\n";

			while (result){

				//std::cout << result->data << " ";
				completedPart.push_back(result->data);
				result = result->pred;

			}

			//std::cout << "\n----------------------------------\n";

			for (int j = 0; j < dists.size(); j++)
				delete dists[j];

		}

		outCompletedFaceParts.push_back(completedPart);

	}



	return outCompletedFaceParts;


}

vector<vector<int>> findInnerContourVertices(Mesh* cF, vector<vector<int>> inFaceParts){

	vector<vector<int>> innerPoints;

	vector<vector<Edge*>> facePartConturEdges;

	for (int i = 0; i < inFaceParts.size(); i++){

		vector<Edge*> conturEdges;
		for (int p = 0; p < inFaceParts[i].size() - 1; p++){

			int curr = inFaceParts[i][p];
			int next = inFaceParts[i][p + 1];

			vector<int> edges = cF->verts[curr]->edgeList;

			for (int e = 0; e < edges.size(); e++){

				int eid = edges[e];
				if ((cF->edges[eid]->v1i == next) || (cF->edges[eid]->v2i == next)) {

					conturEdges.push_back(cF->edges[eid]);
					break;


				}
			}

		}

		facePartConturEdges.push_back(conturEdges);

	}


	vector<vector<int>> edge2Tri;


	for (int i = 0; i < cF->edges.size(); i++){

		//edgeVisited.push_back(false);

		int v1 = cF->edges[i]->v1i;
		int v2 = cF->edges[i]->v2i;

		vector<int> commonTris;
		vector<int> v1Tris = cF->verts[v1]->triList;
		vector<int> v2Tris = cF->verts[v2]->triList;

		for (int t1 = 0; t1 < v1Tris.size(); t1++){

			for (int t2 = 0; t2 < v2Tris.size(); t2++){
				if (v1Tris[t1] == v2Tris[t2])
					commonTris.push_back(v1Tris[t1]);
			}
		}

		edge2Tri.push_back(commonTris);
	}


	for (int fp = 0; fp < faceParts.size(); fp++){
		vector<int> tempInnerPoints;
		std::cout << "---------------------" << std::endl;
		std::cout << fp << std::endl;
		
		
		bool* edgeVisited = new bool[cF->edges.size()];
		
		for (int jj = 0; jj < cF->edges.size(); jj++)
			edgeVisited[jj] = false;

		std::queue<Edge*> eq;
		edgeVisited[facePartConturEdges[fp][0]->idx] = true;



		int triId0 = edge2Tri[facePartConturEdges[fp][0]->idx][0];
		int triId1 = edge2Tri[facePartConturEdges[fp][0]->idx][1];
		int v0 = cF->tris[triId0]->v1i + cF->tris[triId0]->v2i + cF->tris[triId0]->v3i - facePartConturEdges[fp][0]->v1i - facePartConturEdges[fp][0]->v2i;
		int v1 = cF->tris[triId1]->v1i + cF->tris[triId1]->v2i + cF->tris[triId1]->v3i - facePartConturEdges[fp][0]->v1i - facePartConturEdges[fp][0]->v2i;

		float d0 = 0;
		float d1 = 0;
		cout << "dijsktra" << endl;
		Dijsktra::DijsktraSP dsp2;
		dsp2.setMesh(cF);
		std::vector<Dijsktra::Node*> dists0 = dsp2.run(v0);
		std::vector<Dijsktra::Node*> dists1 = dsp2.run(v1);
		cout << "dijsktra1" << endl;
		for (int i = 0; i < inFaceParts[fp].size(); i++){

			if (d0 < dists0[inFaceParts[fp][i]]->key)
				d0 = dists0[inFaceParts[fp][i]]->key;

			if (d1 < dists1[inFaceParts[fp][i]]->key)
				d1 = dists1[inFaceParts[fp][i]]->key;
		}
		cout << "dijsktra2" << endl;
		int triId = (d0 > d1 ? triId1 : triId0);
		std::cout << triId << std::endl;
		int e1 = cF->tris[triId]->e1;
		int e2 = cF->tris[triId]->e2;
		int e3 = cF->tris[triId]->e3;




		bool isContur1 = false;
		bool isContur2 = false;
		bool isContur3 = false;

		for (int i = 0; i < facePartConturEdges[fp].size(); i++){

			//std::cout << conturEdges[i]->idx << std::endl;

			if (facePartConturEdges[fp][i]->idx == e1)
				isContur1 = true;

			if (facePartConturEdges[fp][i]->idx == e2)
				isContur2 = true;

			if (facePartConturEdges[fp][i]->idx == e3)
				isContur3 = true;
		}

		if (!edgeVisited[e1] && !isContur1)
			eq.push(cF->edges[e1]);

		if (!edgeVisited[e2] && !isContur2)
			eq.push(cF->edges[e2]);

		if (!edgeVisited[e3] && !isContur3)
			eq.push(cF->edges[e3]);



		int ttt = 0;
		while (!eq.empty()){
			//std::cout << "------------" << std::endl;
			//std::cout << "0" << std::endl;
			int count = 0;
			Edge* ed = eq.front();
			eq.pop();


			if (edgeVisited[ed->idx])
				continue;
			ttt++;
			//s->SetTuple1(ed->v1i, 10 + fp * 50); s->SetTuple1(ed->v2i, 10 + fp * 50);

			if (find(tempInnerPoints.begin(), tempInnerPoints.end(), ed->v1i) == tempInnerPoints.end())
				tempInnerPoints.push_back(ed->v1i);

			if (find(tempInnerPoints.begin(), tempInnerPoints.end(), ed->v2i) == tempInnerPoints.end())
				tempInnerPoints.push_back(ed->v2i);

			edgeVisited[ed->idx] = true;
			//std::cout << "1" << std::endl;
			for (int tt = 0; tt < edge2Tri[ed->idx].size(); tt++){

				int triId = edge2Tri[ed->idx][tt];

				int e1 = cF->tris[triId]->e1;
				int e2 = cF->tris[triId]->e2;
				int e3 = cF->tris[triId]->e3;

				if (!edgeVisited[e1]){

					bool t1 = false;
					bool t2 = false;
					for (int p = 0; p < inFaceParts[fp].size(); p++){

						if (inFaceParts[fp][p] == cF->edges[e1]->v1i)
							t1 = true;

						if (inFaceParts[fp][p] == cF->edges[e1]->v2i)
							t2 = true;

					}

					if (!(t1&&t2))
						eq.push(cF->edges[e1]);

				}
				if (!edgeVisited[e2]){

					bool t1 = false;
					bool t2 = false;
					for (int p = 0; p < inFaceParts[fp].size(); p++){

						if (inFaceParts[fp][p] == cF->edges[e2]->v1i)
							t1 = true;

						if (inFaceParts[fp][p] == cF->edges[e2]->v2i)
							t2 = true;

					}

					if (!(t1&&t2))
						eq.push(cF->edges[e2]);

				}
				if (!edgeVisited[e3]){

					bool t1 = false;
					bool t2 = false;
					for (int p = 0; p < inFaceParts[fp].size(); p++){

						if (inFaceParts[fp][p] == cF->edges[e3]->v1i)
							t1 = true;

						if (inFaceParts[fp][p] == cF->edges[e3]->v2i)
							t2 = true;

					}
					if (!(t1&&t2))
						eq.push(cF->edges[e3]);

				}


			}

			//std::cout << "2" << std::endl;
			delete ed;


		}

		innerPoints.push_back(tempInnerPoints);
		std::cout << "done" << std::endl;
	

	}

	return innerPoints;



}

vector<Vertex*> completeContour(vector<Vertex*> contourPoints){

	vector<float> dists;
	for (int j = 1; j < contourPoints.size(); j++){

		float dist = sqrt(dist2Between(contourPoints[j - 1]->coords, contourPoints[j]->coords));
		//cout << dist << endl;
		dists.push_back(dist);
	}

	
	int n = NUM_CONTOUR_POINTS - contourPoints.size();
	
	
	for (int i = 0; i < n; i++){

		float maxDist = 0;
		int maxj = 0;
		for (int j = 0; j < contourPoints.size()-1; j++){

			if (dists[j] > maxDist ){

				maxj = j;
				maxDist = dists[j];

			}
		}

		int curr = maxj;
		int prev = curr - 1;

		float* newCoords = new float[3];

		newCoords[0] = (contourPoints[curr]->coords[0] + contourPoints[curr + 1]->coords[0]) / 2;
		newCoords[1] = (contourPoints[curr]->coords[1] + contourPoints[curr + 1]->coords[1]) / 2;
		newCoords[2] = (contourPoints[curr]->coords[2] + contourPoints[curr + 1]->coords[2]) / 2;

		Vertex* v = new Vertex(-i, newCoords);
		//delete newCoords;
		//contourPoints.push_back(v);
		contourPoints.insert(contourPoints.begin() + curr+1, v);
		float dist1 = sqrt(dist2Between(contourPoints[curr]->coords, contourPoints[curr+1]->coords));
		float dist2 = sqrt(dist2Between(contourPoints[curr+1]->coords, contourPoints[curr+2]->coords));
		dists.erase(dists.begin() + curr);
		dists.insert(dists.begin() + curr, dist1);
		dists.insert(dists.begin() + curr+1, dist2);
		//cout << "--------------------------------------------------\n";
	}
	
	return contourPoints;
}

int* distHist(vector<Vertex*> contourPoints, int index){

	
	int n = contourPoints.size();
	float* dists = new float[n];
	//float* angles = new float[3*n];

	int* hist = new int[20];
	for (int i = 0; i < 20; i++)
		hist[i] = 0;

	
	float mean = 0;
	float min = 999999999999;
	for (int i = 0; i < n; i++){
		dists[i] = sqrt(dist2Between(contourPoints[index]->coords, contourPoints[i]->coords));
		mean += dists[i];

		if (dists[i] < min)
			min = dists[i];

		//cout << dists[i] << endl;

	}

	mean /= n;

	for (int i = 0; i < n; i++)
		dists[i] -= min;

	
	float stdDev = 0;

	for (int i = 0; i < n; i++)
		stdDev += dists[i] * dists[i];


	stdDev = sqrt(stdDev / n);

	for (int i = 0; i < n; i++)
		dists[i] /= stdDev;

	float maxd = 0;
	for (int i = 0; i < n; i++)
		if (dists[i] > maxd)
			maxd = dists[i];

	float h = maxd / 20.;
	//cout << "------------------------------------------\n";
	for (int i = 0; i < n; i++){
		//cout << dists[i] << endl;
		int ind = dists[i] /0.2;
		if (ind < 0)
			ind = 0;
		if (ind>19)
			ind = 19;
		hist[ind]++;

	}
	
	//cout << "******************************************\n";
	return hist;
}

int* angleHist(vector<Vertex*> contourPoints, int index){

	int n = contourPoints.size();
	float* angles = new float[3 * n];

	int* hist = new int[60];
	for (int i = 0; i < 60; i++){
		hist[i] = 0;
	}

	float meanx = 0;
	float meany = 0;
	float meanz = 0;
	float minx = 999999;
	float miny = 999999;
	float minz = 999999;

	for (int i = 0; i < n; i++){

		float dx = contourPoints[index]->coords[0] - contourPoints[i]->coords[0];
		float dy = contourPoints[index]->coords[1] - contourPoints[i]->coords[1];
		float dz = contourPoints[index]->coords[2] - contourPoints[i]->coords[2];
		//cout << dists[i] << endl;

		if (dy == 0)
			angles[3*i] = M_PI / 2;
		else
			angles[3*i] = atan(dx / dy) ;

		if (dz == 0)
			angles[3*i+1] = M_PI / 2;
		else
			angles[3*i + 1] = atan(dy / dz) ;

		if (dx == 0)
			angles[3*i+2] = M_PI / 2;
		else
			angles[3*i + 2] = atan(dz / dx) ;

		//angles[i]   = angles[i]   * 180 / M_PI;
		//angles[i+1] = angles[i+1] * 180 / M_PI;
		//angles[i+2] = angles[i+2] * 180 / M_PI;

		meanx += angles[3*i];
		meany += angles[3*i+1];
		meanz += angles[3*i+2];

		if (minx > angles[3 * i])
			minx = angles[3 * i];

		if (miny > angles[3 * i + 1])
			minx = angles[3 * i + 1];

		if (minz > angles[3 * i + 2])
			minz = angles[3 * i + 2];

		//cout << angles[i] << " " << angles[i+1] << " " << angles[i+2] << " ";

	}
	//cout << endl;

	meanx /= n;
	meany /= n;
	meanz /= n;

	for (int i = 0; i < n; i++){

		angles[ 3*i] -= meanx;
		angles[3*i + 1] -= meany;
		angles[ 3*i + 2] -= meanz;
	}

	float maxx = 0;
	float maxy = 0;
	float maxz = 0;

	for (int i = 0; i < n; i++){

		if (maxx < angles[3 * i])
			maxx = angles[3 * i];

		if (maxy < angles[3 * i + 1])
			maxy = angles[3 * i + 1];

		if (maxz < angles[3 * i + 2])
			maxz = angles[3 * i + 2];
	}


	float hsize = M_PI / 20.;

	/*float hx = maxx / 20.;
	float hy = maxy / 20.;
	float hz = maxz / 20.;*/


	for (int i = 0; i < n; i++){
		//cout << angles[i] << endl;

		int indx = angles[3 * i] / hsize + 10 ;
		if (indx < 0)
			indx = 0;
		if (indx>19)
			indx = 19;

		int indy = angles[3 * i + 1] / hsize +10;
		if (indy < 0)
			indy = 0;
		if (indy>19)
			indy = 19;

		int indz = angles[3 * i + 2] / hsize +10;
		if (indz < 0)
			indz = 0;
		if (indz>19)
			indz = 19;

		
		hist[indx]++;
		hist[20+indy]++;
		hist[40+indz]++;

	}

	return hist;
}

float* consecutivePointFeature(vector<Vertex*> contourPoints){

	float* featureVector = new float[3*contourPoints.size()];

	int i = 0;
	for (; i < contourPoints.size()-1; i++)
	{
		featureVector[3 * i]   = (contourPoints[i]->coords[0] - contourPoints[i + 1]->coords[0]);
		featureVector[3 * i+1] = (contourPoints[i]->coords[1] - contourPoints[i + 1]->coords[1]);
		featureVector[3 * i+2] = (contourPoints[i]->coords[2] - contourPoints[i + 1]->coords[2]);
	}

	featureVector[3 * i]     = (contourPoints[i]->coords[0] - contourPoints[0]->coords[0]);
	featureVector[3 * i + 1] = (contourPoints[i]->coords[1] - contourPoints[0]->coords[1]);
	featureVector[3 * i + 2] = (contourPoints[i]->coords[2] - contourPoints[0]->coords[2]);

	return featureVector;
}

vector<int**> extractContourFeatures(Mesh *face, vector<vector<int>> faceParts){

	vector<int**> partHists;

	for (int i = 0; i < faceParts.size(); i++){


		vector<Vertex*> contourPoints;

		for (int j = 0; j < faceParts[i].size(); j++){

			contourPoints.push_back(face->verts[faceParts[i][j]]);

		}

		vector<Vertex*> newPoints = completeContour(contourPoints);

		



		/*******************************************************************************/
		ofstream f1, f2;
		f1.open("part" + std::to_string(i) + ".xyz");
		f2.open("part" + std::to_string(i) + "_added.xyz");

		for (int j = 0; j < contourPoints.size(); j++){

			f1 << contourPoints[j]->coords[0] << " " << contourPoints[j]->coords[1] << " " << contourPoints[j]->coords[2] << endl;
		}

		for (int j = 0; j < newPoints.size(); j++){

			f2 << newPoints[j]->coords[0] << " " << newPoints[j]->coords[1] << " " << newPoints[j]->coords[2] << endl;
		}

		f1.close();
		f2.close();
		/****************************************************************************/

		int** hists = new int*[NUM_CONTOUR_POINTS];

		for (int j = 0; j < NUM_CONTOUR_POINTS; j++){
			int* hist = new int[80];
			int* hist1 = distHist(newPoints, j);
			int* hist2 = angleHist(newPoints, j);


			for (int h = 0; h < 20; h++)
				hist[h] = hist1[h];
			for (int h = 0; h < 60; h++)
				hist[h+20] = hist2[h];
			

			//angleHist(newPoints, j);
			hists[j] = hist;
		}

		/*
		for (int j = 0; j < NUM_CONTOUR_POINTS; j++){

			for (int k = 0; k < 80; k++)
				cout << hists[j][k] << " ";

			cout << endl;
		}
		cout << "------------------------------------------\n";

		*/
		partHists.push_back(hists);

		//for (int j = 0; j < contourPoints.size(); j++)
			//delete contourPoints[i];


		//for (int j = 0; j < newPoints.size(); j++)
			//delete newPoints[i];
	}

	return partHists;
}


vector<float*> extractContourFeatures2(Mesh *face, vector<vector<int>> faceParts){

	vector<float*> partHists;

	for (int i = 0; i < faceParts.size(); i++){


		vector<Vertex*> contourPoints;

		for (int j = 0; j < faceParts[i].size(); j++){

			contourPoints.push_back(face->verts[faceParts[i][j]]);

		}

		vector<Vertex*> newPoints = completeContour(contourPoints);



		float* f= consecutivePointFeature(newPoints);

		
		partHists.push_back(f);

	}

	return partHists;
}

float distContourFeatures(int* f1, int* f2){

	float cost = 0;

	for (int k = 0; k < 80; k++)
	{
		if (k < 20)
			cost += 1.9*(f1[k] - f2[k]) * (f1[k] - f2[k]);
		else
			cost += 0.1*(f1[k] - f2[k]) * (f1[k] - f2[k]);
	}
	

	cost = sqrt(cost);
	cost /= 80;
	return cost;
}

float findClosestPointCost(int* pHist , int** hists){

	float bestCost=999999999999;
	float besti=-1;
	for (int i = 0; i < NUM_CONTOUR_POINTS; i++){

		float cost = distContourFeatures(pHist, hists[i]);

		if (cost < bestCost){

			bestCost = cost;
			besti = i;
		}

	}

	return bestCost;
}

float matchContour(int** hist1, int** hist2){

	float cost = 0;
	for (int i = 0; i < NUM_CONTOUR_POINTS; i++){
	
		cost += findClosestPointCost(hist1[i], hist2);
	}

	cost /= NUM_CONTOUR_POINTS;

	return cost;

}


float matchFaceParts(vector<int**> f1, vector<int**> f2){

	float cost = 0;
	for (int i = 0; i < f1.size(); i++)
		cost += matchContour(f1[i], f2[i]);

	cost /= f1.size();

	return cost;

}

float matchFaceParts(float* f1, float* f2){

	float temp= 0;

	for (int j = 0; j < NUM_CONTOUR_POINTS; j++)
	{
		temp += (f1[j] - f2[j]) * (f1[j] - f2[j]);
	}

	return sqrt(temp);
}


////////////////////////////////////////graph-based features///////////////////////////////////////////////////////

vector<pair<float, int>> findClosestOut(graph::Vertex* v, Vector* target){

	vector<pair<float, int>> outs;
	float min = 999.;
	int minj = -1;
	for (int j = 0; j < v->outs.size(); j+=2){

		

		float diff = v->outs[j]->v->calculateAngleBetween(*(target));

		if (diff < min){

			min = diff;
			minj = j;
		}

		outs.push_back(make_pair(diff, j));

	}

	//cout<<outs.size()<<endl;
	//getchar();
	return outs;


}

void startSearchFrom(graph::Vertex*  v, vector<Vector*> path, vector<graph::Edge*> &candidates, float &cost, int depth , graph::Edge* last){

	//cout<<"		" <<depth << endl;
	
	//cout << "-------------------------" << std::endl;
	if (path.size() == 0)
		return;

	float min = 999.;
	int minj = -1;
	vector<pair<float, int>> closestOut = findClosestOut(v, path[0]);
	vector<graph::Edge*> mincandidates;
	float minCost = 99999.;
	int mini = -1;

	//std::cout << closestOut.size() << std::endl;
	//getchar();

	/*cout << v->id <<"	:	"<< endl;
	for (int i = 0; i < closestOut.size(); i++){


		cout << v->outs[closestOut[i].second]->to->id << endl;
	}
	getchar();
	*/
	
	for (int i = 0; i < closestOut.size(); i++){

		int thresh1 = 10;
		int thresh2 = 15;

		//std::cout << "		" << i << std::endl;
		min = closestOut[i].first;
		minj = closestOut[i].second;


		vector<Vector*> temppath;
		vector<graph::Edge*> tempcandidates; 
		

		
		for (int j = 1; j < path.size(); j++){
		
			Vector* vv = new Vector(path[j]->X, path[j]->Y, path[j]->Z);
			temppath.push_back(vv);
		}

		if (temppath.size() == 0){
			
			
			//cout << "fucking end" << endl;
			//getchar();
			return;
		}

		if (min <= thresh1 && v->outs[minj]->to->id != v->id){

			//candidates.push_back(v->outs[minj]);
			//cost += min;
			//path.erase(path.begin());

			/*for (int d = 0; d < depth+1; d++)
				cout << " ";
			cout << v->outs[minj]->to->id << endl;*/
		
			float tempcost = 0;
			startSearchFrom(v->outs[minj]->to, temppath, tempcandidates, tempcost, depth + 1, v->outs[minj]);

			if (tempcost < minCost){

				mini = i;
				minCost = tempcost;
				mincandidates.clear();
				for (int c = 0; c < tempcandidates.size();c++)
					mincandidates.push_back( tempcandidates[c]);
			}
		}

		
		for (int j = 0; j < temppath.size(); j++)
			delete temppath[j];
		temppath.clear();

		if (min > thresh1 && v->outs[minj]->to->id != v->id && last){

			 float angle = v->outs[minj]->v->calculateAngleBetween(*(last->v));

			 if (angle < thresh2){

				 for (int j = 0; j < path.size(); j++){

					 Vector* vv = new Vector(path[j]->X, path[j]->Y, path[j]->Z);
					 temppath.push_back(vv);
				 }

				 float tempcost = 0;
				 tempcandidates.clear();
				 startSearchFrom(v->outs[minj]->to, temppath, tempcandidates, tempcost, depth, v->outs[minj]);

				 if (tempcost < minCost){

					 mini = i;
					 minCost = tempcost;
					 mincandidates.clear();
					 for (int c = 0; c < tempcandidates.size(); c++)
						 mincandidates.push_back(tempcandidates[c]);
				 }

				 for (int j = 0; j < temppath.size(); j++)
					 delete temppath[j];
				 temppath.clear();
			 }

		}
		 
		



		/*else{

			if (candidates.size()){

				graph::Edge* last = candidates[candidates.size() - 1];
				//pair<float, int> closestOut = findClosestOut(v->outs[minj]->to, path[0]);

				//float diff = closestOut.first;


				float diff = v->outs[minj]->v->calculateAngleBetween(*(last->v));



				//cout << diff << endl;
				if (diff < 10){

					//candidates.push_back(v->outs[minj]);

					for (int can = 0; can < candidates.size(); can++)
						std::cout << candidates[can]->from->id << " ";

					std::cout << std::endl;

					for (int can = 0; can < candidates.size(); can++)
						std::cout << input->verts[candidates[can]->from->id]->coords[0] << " " << input->verts[candidates[can]->from->id]->coords[1] << " " << input->verts[candidates[can]->from->id]->coords[2] << endl;





					//std::cout << "diff:  " << diff << std::endl;

					//last->v->printVector();
					//v->outs[minj]->v->printVector();
					//std::cout << "*************************" << std::endl;

					startSearchFrom(v->outs[minj]->to, path, candidates, cost);


				}

				else{

					std::cout << input->verts[v->outs[minj]->from->id]->coords[0] << " " << input->verts[v->outs[minj]->from->id]->coords[1] << " " << input->verts[v->outs[minj]->from->id]->coords[2] << endl;

					std::cout << "end -->" << " " << diff << std::endl;
					std::cout << last->from->id << " " << last->to->id << " " << v->outs[minj]->from->id << " " << v->outs[minj]->to->id << std::endl;
					cost = 99999.;

				}

			}

			else
				cost = 99999.;
		}*/

	}

	if (mini != -1){
		candidates.push_back(v->outs[mini]);
		//cout << "added\n";

		for (int i = 0; i < mincandidates.size(); i++)
			candidates.push_back(mincandidates[i]);

		cost += closestOut[mini].first;

	}
	//cost += minCost;

	//cout << "cost updated to: " << cost << endl;
}

void search(vector<graph::Vertex*> verts, vector<Vector*> path, vector<vector<graph::Edge*>> &allCandidates, vector<float> &costs){



	for (unsigned int i = 0; i <verts.size(); i++){
		//std::cout << i << endl;
		vector<graph::Edge*> candidates;
		float cost = 0.;

		graph::Vertex* v = verts[i];

		startSearchFrom(v, path, candidates, cost,0,NULL);

		if (candidates.size()){
			allCandidates.push_back(candidates);
			costs.push_back(cost);
		}

	}

	
}

///////////////////////////////////////////// some show bussines////////////////////////////////////////////////////////////////
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

////////////////////////////////// segmentation evaluation ///////////////////////////////////////////

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




	Mesh* face1 = new Mesh();
	face1->loadOff("faces\\testface1.off");

	Mesh* face2 = new Mesh();
	face2->loadOff("faces\\face1.off");

	
	int i = 0;

	vector<Vertex*> contourPoints1;
	for (int j = 0; j < faceParts[i].size(); j++){

		contourPoints1.push_back(face1->verts[faceParts[i][j]]);
	}

	vector<Vertex*> newPoints1 = completeContour(contourPoints1);


	vector<Vertex*> contourPoints2;
	for (int j = 0; j < faceParts[i].size(); j++){

		contourPoints2.push_back(face2->verts[faceParts[i][j]]);
	}

	vector<Vertex*> newPoints2 = completeContour(contourPoints2);

	
	ofstream f1;
	f1.open("1.xyz");

	for (int j = 0; j < newPoints1.size(); j++){

		f1 << newPoints1[j]->coords[0] << " " << newPoints1[j]->coords[1] << " " << newPoints1[j]->coords[2] << endl;
	}
	f1.close();


	ofstream f2;
	f2.open("2.xyz");

	for (int j = 0; j < newPoints2.size(); j++){

		f2 << newPoints2[j]->coords[0] << " " << newPoints2[j]->coords[1] << " " << newPoints2[j]->coords[2] << endl;
	}
	f2.close();


	nonRigidAlignment(contourPoints1, contourPoints2);


	ofstream f3;
	f3.open("1_modified.xyz");

	for (int j = 0; j < newPoints1.size(); j++){

		f3 << newPoints1[j]->coords[0] << " " << newPoints1[j]->coords[1] << " " << newPoints1[j]->coords[2] << endl;
	}
	f3.close();

	getchar();
	return 5;
	srand(time(NULL));
	faceNo = 3;
	name = "initial";
	initialFace = new Mesh();
	/*
	if (!pointCloudMode)
		initialFace->loadOff("faces\\testface3.off");
	else 
		initialFace->loadxyz("faces\\face1.xyz");
	*/
	//initialFace->scale(0.001);
	// s = initialFace->calculateScale();
	 cube = initialFace->getVTKPolyData(!pointCloudMode);
	
	int nVerts = cube->GetPoints()->GetNumberOfPoints();

	
	isSelected = new bool[nVerts];
	segments.clear();
	for (int i = 0; i < nVerts; i++){
		isSelected[i] = false;
		segments.push_back(-1);

	}

	//initialFace->findNeighborhoodTriangles();
	//initialFace->assignNormalsToTriangles();
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
