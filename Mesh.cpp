#include "Mesh.h"

#include <vtkPointSource.h>
#include <vtkSmartPointer.h>
void Mesh::loadOff(char* meshFile)
{
	
	center.setZero();
	cout << "Mesh initializing (to " << meshFile << ")..\n";

	FILE* fPtr;
	if (! (fPtr = fopen(meshFile, "r")))
	{
		cout << "cannot read " << meshFile << endl;
		exit(0);
	}

	char off[25];
	fscanf(fPtr, "%s\n", &off); //cout << off << " type file\n";
	float a, b, c, d;	//for face lines and the 2nd line (that gives # of verts, faces, and edges) casting these to int will suffice
	fscanf(fPtr, "%f %f %f\n", &a, &b, &c);
	int nVerts = (int) a, v = 0;
/*	minEdgeLen = INF;
	maxEdgeLen = -INF;
	edgeLenTotal = 0.0f;*/
	while (v++ < nVerts) //go till the end of verts coord section
	{
		fscanf(fPtr, "%f %f %f\n", &a, &b, &c);
		

		float* coords = new float[3];
		coords[0] = a;
		coords[1] = b;
		coords[2] = c;

		center.x() += a;
		center.y() += b;
		center.z() += c;
		
		addVertex(coords);
	}
	//verts ready, time to fill triangles
	while (fscanf(fPtr, "%f %f %f %f\n", &d, &a, &b, &c) != EOF) //go till the end of file
	{
		addTriangle((int) a, (int) b, (int) c); //no -1 'cos idxs start from 0 for off files
	}

	center.x() /= nVerts;
	center.y() /= nVerts;
	center.z() /= nVerts;
	//std::cout << center << std::endl;
//	avgEdgeLen = edgeLenTotal / ((int) edges.size());
//	computeBoundingBox();
	fclose(fPtr);

	//cout << "Mesh has " << (int)tris.size() << " tris, " << (int)verts.size() << " verts, " << (int)edges.size() << " edges\nInitialization done\n" << endl;
}

void Mesh::loadxyz(char* fName){
	cout << "Mesh initializing (to " << fName << ")..\n";
	center.setZero();
	ifstream f;
	f.open(fName);
	int numOfVerts;
	f >> numOfVerts;
	cout << numOfVerts << endl;
	for (int i = 0; i < numOfVerts; i++){


		float a, b, c;
		f >> a;
		f >> b;
		f >> c;

		float* coords = new float[3];
		coords[0] = a;
		coords[1] = b;
		coords[2] = c;

		center.x() += a;
		center.y() += b;
		center.z() += c;

		addVertex(coords);

	}

	center.x() /= numOfVerts;
	center.y() /= numOfVerts;
	center.z() /= numOfVerts;

	kdtree* kd = kd_create(3);
	kdres* resultSet = NULL;
	int* voxelIdx;
	float voxelCoord[3];

	for (int i = 0; i < this->verts.size(); i++)
	{

		int* tmp = new int[1];
		*tmp = i;
		kd_insert3f(kd, this->verts[i]->coords[0], this->verts[i]->coords[1], this->verts[i]->coords[2], tmp);
		//delete tmp;
	}


	for (int i = 0; i < this->verts.size(); i++)
	{
		
		resultSet = kd_nearest_range3(kd, this->verts[i]->coords[0], this->verts[i]->coords[1], this->verts[i]->coords[2],150);
		int size = kd_res_size(resultSet);
		
		int n = 0;
		while (n < 10 && n<size){

			voxelIdx = (int*)kd_res_itemf(resultSet, voxelCoord);
			this->verts[i]->vertList.push_back(*voxelIdx);
			n++;
			kd_res_next(resultSet);
		}

		kd_res_free(resultSet);

	}

	
	kd_clear(kd);
	kd_free(kd);

}


void Mesh::toOrigin(){

	center.setZero();

	for (int i = 0; i < verts.size(); i++){

		center.x() += verts[i]->coords[0];
		center.y() += verts[i]->coords[1];
		center.z() += verts[i]->coords[2];
	}


	center.x() /= verts.size();
	center.y() /= verts.size();
	center.z() /= verts.size();

	for (int i = 0; i < verts.size(); i++){
	
		verts[i]->coords[0] -= center.x();
		verts[i]->coords[1] -= center.y();
		verts[i]->coords[2] -= center.z();
	
	
	}


}


void Mesh::calculateAreas(){

	for (int i = 0; i < tris.size(); i++){

		float* a = verts[tris[i]->v1i]->coords;
		float* b = verts[tris[i]->v2i]->coords;
		float* c = verts[tris[i]->v3i]->coords;

		float ab[3];
		float ac[3];

		for (int i = 0; i < 3;i++)
			ab[i] = a[i] - b[i];

		for (int i = 0; i < 3; i++)
			ac[i] = a[i] - c[i];

		float x1 = ab[0];
		float x2 = ab[1];
		float x3 = ab[2];

		float y1 = ac[0];
		float y2 = ac[1];
		float y3 = ac[2];

		tris[i]->area = sqrt(pow(x2*y3 - x3*y2, 2) + pow(x1*y3 - x3*y1, 2) + pow(x1*y2 - x2*y1, 2)) / 2;


	}


}


vtkPolyData* Mesh::getVTKPolyData(bool isMesh){

	int nVerts = this->verts.size();
	int nTris = this->tris.size();
	int *triangles;
	vtkPolyData* cube = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();
	vtkFloatArray *scalars = vtkFloatArray::New();

	// Create the topology of the point (a vertex)
	vtkSmartPointer<vtkCellArray> vertices =
		vtkSmartPointer<vtkCellArray>::New();
	

	// Load the point, cell, and data attributes.
	for (int i = 0; i < nVerts; i++) {

		double p[3];
		p[0] = (double)this->verts[i]->coords[0];
		p[1] = (double)this->verts[i]->coords[1];
		p[2] = (double)this->verts[i]->coords[2];
		points->InsertPoint(i, p);
		vertices->InsertNextCell(1, &i);

	}

	for (int i = 0; i < nTris; i++) {

		triangles = new int[3];
		triangles[0] = this->tris[i]->v1i;
		triangles[1] = this->tris[i]->v2i;
		triangles[2] = this->tris[i]->v3i;
		polys->InsertNextCell(3, triangles);
		delete triangles;

	}
	for (int i = 0; i<nVerts; i++) scalars->InsertTuple1(i, 250);

	
	// We now assign the pieces to the vtkPolyData.
	cube->SetPoints(points);
	points->Delete();
	if (isMesh)
		cube->SetPolys(polys);
	polys->Delete();
	cube->GetPointData()->SetScalars(scalars);
	scalars->Delete();
	if (!isMesh)
		cube->SetVerts(vertices);
	
	return cube;

}



void Mesh::createCube(float sideLength)
{
	float** coords = new float*[8], flbc[3] = {10, 7, 5}, delta[3] = {0, 0, 0};
	for (int v = 0; v < 8; v++)
	{
		coords[v] = new float[3];

		if (v == 1)
			delta[0] += sideLength;
		else if (v == 2)
			delta[1] += sideLength;
		else if (v == 3)
			delta[0] -= sideLength;
		else if (v == 4)
			delta[2] -= sideLength;
		else if (v == 5)
			delta[0] += sideLength;
		else if (v == 6)
			delta[1] -= sideLength;
		else if (v == 7)
			delta[0] -= sideLength;

		for (int c = 0; c < 3; c++)
			coords[v][c] = flbc[c] + delta[c];

		addVertex(coords[v]);
	}

	//connectivity
	addTriangle(0, 1, 2);
	addTriangle(0, 2, 3);

	addTriangle(1, 6, 2);
	addTriangle(6, 5, 2);

	addTriangle(7, 5, 6);
	addTriangle(7, 4, 5);

	addTriangle(0, 3, 7);
	addTriangle(3, 4, 7);

	addTriangle(0, 6, 1);
	addTriangle(0, 7, 6);

	addTriangle(3, 2, 5);
	addTriangle(3, 5, 4);
}

void Mesh::addVertex(float* coords)
{
	int idx = verts.size();
	verts.push_back( new Vertex(idx, coords) );



}

void Mesh::addTriangle(int v1i, int v2i, int v3i)
{
	int idx = tris.size();
	tris.push_back( new Triangle(idx, v1i, v2i, v3i) );

	verts[v1i]->triList.push_back( idx );
	verts[v2i]->triList.push_back( idx );
	verts[v3i]->triList.push_back( idx );

	if (! makeVertsNeighbors(v1i, v2i) )
		addEdge(v1i, v2i);

	if (! makeVertsNeighbors(v1i, v3i) )
		addEdge(v1i, v3i);

	if (! makeVertsNeighbors(v2i, v3i) )
		addEdge(v2i, v3i);
}


void Mesh::calculateDihedrals(){


	for (int i = 0; i < this->tris.size(); i++){


		Triangle* t = this->tris[i];

		float maxAngle = 0;
		for (int j = 0; j < t->neighborhood.size(); j++){


			Triangle* neighbor = t->neighborhood[j];

			float angle = t->tNormal.calculateAngleBetween(neighbor->tNormal);

			if (angle > maxAngle){

				maxAngle = angle;

			}

		}

		t->dihedral = maxAngle;


		         
	}

}

bool Mesh::makeVertsNeighbors(int v, int w)
{
	//try to make v and w neighbor; return true if they already are

	for (int check = 0; check < (int) verts[v]->vertList.size(); check++)
		if (verts[ v ]->vertList[ check ] == w)
			return true;

	verts[v]->vertList.push_back(w);
	verts[w]->vertList.push_back(v);
	return false;
}

inline float distanceBetween(float* a, float* b)
{
	return (float) sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) );
}

void Mesh::addEdge(int a, int b)
{
	int idx = edges.size();
	edges.push_back( new Edge(idx, a, b, distanceBetween(verts[a]->coords, verts[b]->coords) ) );

	verts[a]->edgeList.push_back(idx);
	verts[b]->edgeList.push_back(idx);

	verts[a]->vertList.push_back(b);
	verts[b]->vertList.push_back(a);
}


void Mesh::assignNormalsToTriangles()
{
	for (int i = 0; i < this->tris.size(); i++)
	{
		Vector pva; pva.X = this->verts[this->tris[i]->v1i]->coords[0]; pva.Y = this->verts[this->tris[i]->v1i]->coords[1]; pva.Z = this->verts[this->tris[i]->v1i]->coords[2];
		Vector pvb; pvb.X = this->verts[this->tris[i]->v2i]->coords[0]; pvb.Y = this->verts[this->tris[i]->v2i]->coords[1]; pvb.Z = this->verts[this->tris[i]->v2i]->coords[2];
		Vector pvc; pvc.X = this->verts[this->tris[i]->v3i]->coords[0]; pvc.Y = this->verts[this->tris[i]->v3i]->coords[1]; pvc.Z = this->verts[this->tris[i]->v3i]->coords[2];

		Vector normal;
		Vector lhs;
		Vector rhs;
		lhs.X = pvb.X - pva.X;
		lhs.Y = pvb.Y - pva.Y;
		lhs.Z = pvb.Z - pva.Z;

		rhs.X = pvc.X - pva.X;
		rhs.Y = pvc.Y - pva.Y;
		rhs.Z = pvc.Z - pva.Z;

		lhs.crossProduct(rhs, normal);
		normal.Normalize();

		this->tris[i]->tNormal = normal;
		this->tris[i]->normalId = -1;

	}
}



void Mesh::findNeighborhoodTriangles()
{
	assignEdgesToTriangles();

	for (unsigned int i = 0; i < tris.size(); i++) {
		int re1 = tris[i]->e1; int re2 = tris[i]->e2; int re3 = tris[i]->e3;
		for (unsigned int j = 0; j < tris.size(); j++) {
			int se1 = tris[j]->e1; int se2 = tris[j]->e2; int se3 = tris[j]->e3;
			if (i != j) {
				int count = 0;
				if (re1 == se1) count++;
				if (re1 == se2) count++;
				if (re1 == se3) count++;
				if (re2 == se1) count++;
				if (re2 == se2) count++;
				if (re2 == se3) count++;
				if (re3 == se1) count++;
				if (re3 == se2) count++;
				if (re3 == se3) count++;


				if (count >0){
					tris[i]->neighborhood.push_back(tris[j]);

					
				}
			}
		}
	}
}

int Mesh::findEdgeBetween(int v, int w) {
	for (unsigned int i = 0; i < edges.size(); i++) {
		if ((edges[i]->v1i == v && edges[i]->v2i == w) || (edges[i]->v1i == w && edges[i]->v2i == v))
		{
			return edges[i]->idx;
		}
	}
	return -1;
}

void Mesh::assignEdgesToTriangles()
{
	for (unsigned int i = 0; i < tris.size(); i++) {
		int a = tris[i]->v1i;
		int b = tris[i]->v2i;
		int c = tris[i]->v3i;

		int eid1 = findEdgeBetween(a, b);
		int eid2 = findEdgeBetween(a, c);
		int eid3 = findEdgeBetween(b, c);

		tris[i]->e1 = eid1;
		tris[i]->e2 = eid2;
		tris[i]->e3 = eid3;
	}
}



void Mesh::addExtraEdgesToTheMesh( float threshold)
{
	
	std::vector<std::pair<unsigned int, unsigned int> > newEdges;

	for (unsigned int i = 0; i < this->tris.size(); i++) {
		Triangle *refTri = this->tris[i];
		std::vector<Triangle*> planar_neighborhood;
		planar_neighborhood.push_back(refTri);
		for (unsigned int j = 0; j < refTri->neighborhood.size(); j++) {
			Triangle *srcTri = refTri->neighborhood[j];

			Vector n1 = refTri->tNormal;
			Vector n2 = srcTri->tNormal;

			float dot = n1.dotProduct(n2);
			float angle = acos(dot);
			if (angle <threshold) //Means that they are planar
			{
			
				planar_neighborhood.push_back(srcTri);
			}

		}


		
			for (unsigned int j = 0; j < planar_neighborhood.size(); j++){

				int uncommonVertex = -1;
				int uncommonVertex2 =- 1;
				Triangle *neighbor = planar_neighborhood[j];
				int re1 = refTri->v1i; int re2 = refTri->v2i; int re3 = refTri->v3i;
				int se1 = neighbor->v1i; int se2 = neighbor->v2i; int se3 = neighbor->v3i;

				if (se1 != re1 && se1 != re2 && se1 != re3)
					uncommonVertex = se1;
				else if (se2 != re1 && se2 != re2 && se2 != re3)
					uncommonVertex = se2;
				else if (se3 != re1 && se3 != re2 && se3 != re3)
					uncommonVertex = se3;

				if (re1 != se1 && re1 != se2 && re1 != se3)
					uncommonVertex2 = re1;
				else if (re2 != se1 && re2 != se2 && re2 != se3)
					uncommonVertex2 = re2;
				else if (re3 != se1 && re3 != se2 && re3 != se3)
					uncommonVertex2 = re3;

				if (uncommonVertex != -1 && uncommonVertex2 != -1)
					this->addEdge(uncommonVertex, uncommonVertex2);

			}


		

	}

	
}






double cotangentWeight(Mesh* md, int iv, int jv) {
	
	if (md->tris.size() == 0)
		return 1;

	int noTri = 0;
	double cots = 0.0;
	for (size_t i = 0; i < md->verts[iv]->triList.size(); i++) {
		Triangle * ivTri = md->tris[md->verts[iv]->triList[i]];
		int b, n1 = iv, n2 = jv;
		bool shareTriangle = false;

		if (ivTri->v1i == jv) {
			noTri++;
			shareTriangle = true;
			if (ivTri->v2i == iv) {
				b = ivTri->v3i;
			}
			else {
				b = ivTri->v2i;
			}
		}
		else if (ivTri->v2i == jv) {
			noTri++;
			shareTriangle = true;
			if (ivTri->v1i == iv) {
				b = ivTri->v3i;
			}
			else {
				b = ivTri->v1i;
			}
		}
		else if (ivTri->v3i == jv) {
			noTri++;
			shareTriangle = true;
			if (ivTri->v1i == iv) {
				b = ivTri->v2i;
			}
			else {
				b = ivTri->v1i;
			}
		}
		if (shareTriangle) {
			Vertex *vn1 = md->verts[n1];
			Vertex *vn2 = md->verts[n2];
			Vector base; base.X = md->verts[b]->coords[0]; base.Y = md->verts[b]->coords[1]; base.Z = md->verts[b]->coords[2];
			Vector ne1; ne1.X = vn1->coords[0]; ne1.Y = vn1->coords[1]; ne1.Z = vn1->coords[2];
			Vector ne2; ne2.X = vn2->coords[0]; ne2.Y = vn2->coords[1]; ne2.Z = vn2->coords[2];

			Vector v1 = ne1.Subtract(base); v1.Normalize();
			Vector v2 = ne2.Subtract(base); v2.Normalize();

			float dot = v1.dotProduct(v2);

			float angle = acos(dot);

			cots += 1.0 / tan(angle);
		}

	}
	return cots / noTri;
}



void Mesh::initialize(){
	tripletList.clear();
	for (int i = 0; i < (int)this->verts.size(); i++)
	{
		double sum_w = 0.0;
		for (int j = 0; j < (int)this->verts[i]->vertList.size(); j++) {
			double w_ij = cotangentWeight(this, i, this->verts[i]->vertList[j]);
			//w_ij = 1.0; //Compute cotanjant weights here
			sum_w += w_ij;
		}
		for (int j = 0; j < (int)this->verts[i]->vertList.size(); j++) {
			double w_ij =  cotangentWeight(this, i, this->verts[i]->vertList[j]);
			//w_ij = 1.0; //Compute cotanjant weights here

			if (i != this->verts[i]->vertList[j]) {
				double entry = double((-1.0)*w_ij / sum_w);
				tripletList.push_back(Trip(i, this->verts[i]->vertList[j], entry));
			}
		}
		tripletList.push_back(Trip(i, i, 1.0));
	}
	int LapSize = this->verts.size();
	L = SpMat(LapSize, LapSize);
	L.setFromTriplets(tripletList.begin(), tripletList.end());
	cx = Eigen::VectorXd(LapSize);
	cy = Eigen::VectorXd(LapSize);
	cz = Eigen::VectorXd(LapSize);
	vx = Eigen::VectorXd(LapSize);
	vy = Eigen::VectorXd(LapSize);
	vz = Eigen::VectorXd(LapSize);
	dx = Eigen::VectorXd(LapSize);
	dy = Eigen::VectorXd(LapSize);
	dz = Eigen::VectorXd(LapSize);
	LTLx = Eigen::VectorXd(LapSize);
	LTLy = Eigen::VectorXd(LapSize);
	LTLz = Eigen::VectorXd(LapSize);
}


void Mesh::deform(Mesh* target,float weigth)
{

	target->initialize();
	
	int LapSize = this->verts.size();
	

	SpMat W = SpMat(LapSize, LapSize);
	W.setIdentity();
	
	kdres* resultSet = NULL;
	kdtree* kd = kd_create(3);
	
	int* voxelIdx;
	float voxelCoord[3];
	
	for (int i = 0; i < target->verts.size(); i++)
	{

		int* tmp = new int[1];
		*tmp = i;
		kd_insert3f(kd, target->verts[i]->coords[0], target->verts[i]->coords[1], target->verts[i]->coords[2], tmp);
		delete tmp;
	}


	for (int i = 0; i < this->verts.size(); i++)
	{

	

		resultSet = kd_nearest3f(kd, this->verts[i]->coords[0], this->verts[i]->coords[1], this->verts[i]->coords[2]);
		voxelIdx = (int*)kd_res_itemf(resultSet, voxelCoord);
		
		cx[i] = voxelCoord[0];
		cy[i] = voxelCoord[1];
		cz[i] = voxelCoord[2];


		//cx[i] = target->verts[i]->coords[0];
		//cy[i] = target->verts[i]->coords[1];
		//cz[i] = target->verts[i]->coords[2];

		vx[i] = this->verts[i]->coords[0];
		vy[i] = this->verts[i]->coords[1];
		vz[i] = this->verts[i]->coords[2];
		kd_res_free(resultSet);
	
		

		
	}

	kd_clear(kd);
	kd_free(kd);
	
	
	//SpMat L_ = target->L;


	SpMat A = W + weigth*(L.transpose()*L);

	LTLx = (L.transpose()*L)*vx;
	LTLy = (L.transpose()*L)*vy;
	LTLz = (L.transpose()*L)*vz;

	bx = W*cx + weigth*LTLx;
	by = W*cy + weigth*LTLy;
	bz = W*cz + weigth*LTLz;

	Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of LpTLp
	vx = chol.solve(bx);
	vy = chol.solve(by);
	vz = chol.solve(bz);
	
	
	for (int i = 0; i < LapSize; i++){

		 this->verts[i]->coords[0] = vx[i];
		 this->verts[i]->coords[1] = vy[i];
		 this->verts[i]->coords[2] = vz[i];

	}

	
	
	
}


void Mesh::write(char* filename, bool isMesh){

	if (isMesh){
		std::ofstream of;
		of.open(filename);
		of << "OFF" << endl;
		of << this->verts.size() << " " << this->tris.size() << " 0" << endl;

		for (int i = 0; i < this->verts.size(); i++)
			of << this->verts[i]->coords[0] << " " << this->verts[i]->coords[1] << " " << this->verts[i]->coords[2] << endl;

		for (int i = 0; i < this->tris.size(); i++)
			of << "3 " << this->tris[i]->v1i << " " << this->tris[i]->v2i << " " << this->tris[i]->v3i << endl;


		of.close();
	}
	else {

		std::ofstream of;
		of.open(filename);
		
		of << this->verts.size() << endl;

		for (int i = 0; i < this->verts.size(); i++)
			of << this->verts[i]->coords[0] << "	" << this->verts[i]->coords[1] << "		" << this->verts[i]->coords[2] << endl;

		of.close();


	}



}




void Mesh::deform(Mesh* target, float weigth, int iterations){



	for (int i = 0; i < iterations; i++){
		cout << i << endl;
		deform(target, weigth);
	}


}


float Mesh::calculateScale(){

	float dist = 0;
	for (int i = 0; i < verts.size(); i++){

		dist += sqrt(pow(verts[i]->coords[0] - center.x(), 2) + pow(verts[i]->coords[1] - center.y(), 2) + pow(verts[i]->coords[2] - center.z(), 2));


	}

	return dist / verts.size();




}


void Mesh::scale(float s){


	for (int i = 0; i < verts.size(); i++){
		
		verts[i]->coords[0] *= s;
		verts[i]->coords[1] *= s;
		verts[i]->coords[2] *= s;
	
	}

}


int Mesh::smoothNeighborsNormal(int tid, int numOfProcessed){


	Triangle* t = this->tris[tid];
	Vector n1 = t->tNormal;
	vector<int> availableNeighbors;
	for (int i = 0; i < t->neighborhood.size(); i++){


		Triangle* t2 = t->neighborhood[i];
		if (t2->normalId == -1){

			Vector n2 = t2->tNormal;

			float angle = n2.calculateAngleBetween(n1);


			if (angle < 8){

				//t2->tNormal = n1;
				t2->normalId = t->normalId;
				numOfProcessed++;
				availableNeighbors.push_back(i);
				
			}

		}

	}

	for (int i = 0; i < availableNeighbors.size(); i++){
		
		
		numOfProcessed = smoothNeighborsNormal(t->neighborhood[availableNeighbors[i]]->idx, numOfProcessed);


	}

	return numOfProcessed;

}


int Mesh::smoothAllNormals(int tid){


	int lastId = 0;
	this->tris[tid]->normalId = lastId;

	int numOfProcessed = 1;
	
	
	
	while (true) {


		bool found = false;

		smoothNeighborsNormal(tid, numOfProcessed);

		std::cout << numOfProcessed << std::endl;

		lastId++;


		for (int i = 0; i < this->tris.size(); i++){

			if (this->tris[i]->normalId == -1){
				
				this->tris[i]->normalId = lastId;
				numOfProcessed++;
				numOfProcessed = smoothNeighborsNormal(i, numOfProcessed);
				found = true;
				break;

			}

		}


		if (!found)
			break;

	}


	std::cout << lastId << std::endl;

	return lastId;

}

void Mesh::angleXY(){

	float minx = 999;
	float miny = 999;
	float minz = 999;
	float maxx = -999;
	float maxy = -999;
	float maxz = -999;
	for (int i = 0; i < this->tris.size(); i++){

		Triangle* t = this->tris[i];


		float x  = t->tNormal.calculateAngleBetween(Vector(1, 0, 0));
		float y  = t->tNormal.calculateAngleBetween(Vector(0, 1, 0));
		float z  = t->tNormal.calculateAngleBetween(Vector(0, 0, 1));


		if (x < minx)
			minx = x;
		if (x > maxx)
			maxx = x;

		
		if (y < miny)
			miny = y;
		if (y > maxy)
			maxy = y;


		if (z < minz)
			minz = z;
		if (z > maxz)
			maxz = z;




		t->axisAngles = Vector(x, y, z);




	}


	std::cout << minx <<"	"<< miny<<"	" << minz << std::endl;
	std::cout << maxx << "	" << maxy << "	" << maxz << std::endl;

}
