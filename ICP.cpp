#include "ICP.hpp"

void ICP::align(Mesh *r,Mesh *im){
	
	reference = new Mesh();
	input = new Mesh();
	
	for (int v = 0; v < r->verts.size(); v++){
		float* vertex = new float[3];
		vertex[0] = r->verts[v]->coords[0];
		vertex[1] = r->verts[v]->coords[1];
		vertex[2] = r->verts[v]->coords[2];
		
		//cout << vertex[0] << "   " << vertex[1] << "    " << vertex[2] << endl;
		reference->addVertex(vertex);
	}


	for (int v = 0; v < im->verts.size(); v++){
		float* vertex = new float[3];
		vertex[0] = im->verts[v]->coords[0];
		vertex[1] = im->verts[v]->coords[1];
		vertex[2] = im->verts[v]->coords[2];
		//cout << vertex[0] << "   " << vertex[1] << "    " << vertex[2] << endl;
		input->addVertex(vertex);
	}

	Transformations.clear();
	globalTrans = Eigen::Matrix4f::Identity();
	globalTrans1 = Eigen::Matrix4f::Identity();
	globalTrans2= Eigen::Matrix4f::Identity();

	Vector comR(0,0,0);
	Vector comT(0, 0, 0);
	for (int i = 0; i < reference->verts.size(); i++){
	
		comR.X += reference->verts[i]->coords[0];
		comR.Y += reference->verts[i]->coords[1];
		comR.Z += reference->verts[i]->coords[2];
	
	}
		
	
	comR = comR.Multiply(1.0/ reference->verts.size());
	
	
	for (int i = 0; i < input->verts.size(); i++){
	
		comT.X += input->verts[i]->coords[0];
		comT.Y += input->verts[i]->coords[1];
		comT.Z += input->verts[i]->coords[2];
	
	}
		

	comT = comT.Multiply(1.0 / input->verts.size());
	
	
	Eigen::MatrixXf ms1(3, input->verts.size()), ms2(3, input->verts.size());
	Eigen::Matrix3f cov12; 
	Eigen::Matrix4f Q; 

	kdtree* kd = kd_create(3);
	kdres* resultSet = NULL;
	int* voxelIdx;
	float voxelCoord[3];

	for (int i = 0; i < reference->verts.size(); i++)
	{
		
		int* tmp = new int[1];
		*tmp = i;
		kd_insert3f(kd, reference->verts[i]->coords[0], reference->verts[i]->coords[1], reference->verts[i]->coords[2], tmp);
	}
	
	int nIters = 0, maxIters = 100;
	float mse = 999999999.0f, prevMse = 0.0f, 
		  alignErr = 100.0f, closeness = 0.0000001f;

	while (++nIters <= maxIters)
	{

		comT =Vector(0,0,0);
		for (int i = 0; i < input->verts.size(); i++){

			comT.X += input->verts[i]->coords[0];
			comT.Y += input->verts[i]->coords[1];
			comT.Z += input->verts[i]->coords[2];

	

		}

		comT = comT.Multiply(1.0 / input->verts.size());
		
		
		std::vector< int > matchedVoxels; 
		std::vector<int> matchedVertices;
		for (int i = 0; i<reference->verts.size(); i++)
			matchedVertices.push_back(-1);

		for (int i = 0; i < input->verts.size(); i++) 
		{

			resultSet = kd_nearest3f(kd, input->verts[i]->coords[0], input->verts[i]->coords[1], input->verts[i]->coords[2]); 
			voxelIdx = (int*) kd_res_itemf(resultSet, voxelCoord); 

			for (int k = 0; k < reference->verts.size(); k++){
				if (reference->verts[k]->coords[0] == voxelCoord[0] && reference->verts[k]->coords[1] == voxelCoord[1] && reference->verts[k]->coords[2] == voxelCoord[2]){
					*voxelIdx = k;
					break;
				}
			}


			matchedVertices[*voxelIdx] = i;
			matchedVoxels.push_back(*voxelIdx);
		}

		
		for (int i = 0; i < 3; i++) 
			for (int j = 0; j < input->verts.size(); j++) 
				ms1(i, j) = input->verts[matchedVertices[matchedVoxels[j]]]->coords[i] - comT[i];


		comR = Vector(0,0,0);
		for (int v = 0; v < input->verts.size(); v++){


			comR.X += reference->verts[matchedVoxels[v]]->coords[0];
			comR.Y += reference->verts[matchedVoxels[v]]->coords[1];
			comR.Z += reference->verts[matchedVoxels[v]]->coords[2];

		}
			
		comR = comR.Multiply(1.0/ input->verts.size());
		
		
		for (int i = 0; i < 3; i++) 
			for (int j = 0; j < input->verts.size(); j++) 
				ms2(i, j) = reference->verts[matchedVoxels[j]]->coords[i] - comR[i];
		
		Eigen::MatrixXf ms2t = ms2.transpose(); 
		cov12 = ms1 * ms2t; 

		Eigen::Matrix4f Transformation;

		if (!svdBased){
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					cov12(i, j) /= (float)input->verts.size();

			Q(0, 0) = cov12(0, 0) + cov12(1, 1) + cov12(2, 2);
			Q(1, 0) = Q(0, 1) = cov12(1, 2) - cov12(2, 1); 
			Q(2, 0) = Q(0, 2) = cov12(2, 0) - cov12(0, 2); 
			Q(3, 0) = Q(0, 3) = cov12(0, 1) - cov12(1, 0); 
			Q(1, 1) = 2 * cov12(0, 0) - Q(0, 0);		Q(1, 2) = cov12(0, 1) + cov12(1, 0);	Q(1, 3) = cov12(0, 2) + cov12(2, 0);
			Q(2, 1) = cov12(1, 0) + cov12(0, 1);	Q(2, 2) = 2 * cov12(1, 1) - Q(0, 0);		Q(2, 3) = cov12(1, 2) + cov12(2, 1);
			Q(3, 1) = cov12(2, 0) + cov12(0, 2);	Q(3, 2) = cov12(2, 1) + cov12(1, 2);	Q(3, 3) = 2 * cov12(2, 2) - Q(0, 0);
			
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> eigensolver(Q); 
			if (eigensolver.info() != Eigen::Success) abort();
			Eigen::Vector4f eig_vals = eigensolver.eigenvalues(); 
			Eigen::Matrix4f eig_vecs = eigensolver.eigenvectors(); 

			float q0 = eig_vecs(0, 3), q1 = eig_vecs(1, 3), q2 = eig_vecs(2, 3), q3 = eig_vecs(3, 3);

			// rotation
			Transformation(0, 0) = q0*q0 + q1*q1 - q2*q2 - q3*q3;	Transformation(0, 1) = 2 * (q1*q2 - q0*q3);				Transformation(0, 2) = 2 * (q1*q3 + q0*q2);
			Transformation(1, 0) = 2 * (q1*q2 + q0*q3);				Transformation(1, 1) = q0*q0 + q2*q2 - q1*q1 - q3*q3;	Transformation(1, 2) = 2 * (q2*q3 - q0*q1);
			Transformation(2, 0) = 2 * (q1*q3 - q0*q2);				Transformation(2, 1) = 2 * (q2*q3 + q0*q1);				Transformation(2, 2) = q0*q0 + q3*q3 - q1*q1 - q2*q2;
			Transformation(3, 0) = 0;								Transformation(3, 1) = 0;								Transformation(3, 2) = 0;
			

		}
		else{

			Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov12, Eigen::ComputeThinU | Eigen::ComputeThinV);
			Eigen::MatrixXf U = svd.matrixU();
			Eigen::MatrixXf V = svd.matrixV();
			Eigen::MatrixXf R = V * U.transpose();

			// rotation
			Transformation(0, 0) = R(0, 0);	Transformation(0, 1) = R(0, 1);	Transformation(0, 2) = R(0, 2);
			Transformation(1, 0) = R(1, 0);	Transformation(1, 1) = R(1, 1);	Transformation(1, 2) = R(1, 2);
			Transformation(2, 0) = R(2, 0);	Transformation(2, 1) = R(2, 1);	Transformation(2, 2) = R(2, 2);
			Transformation(3, 0) = 0;		Transformation(3, 1) = 0;		Transformation(3, 2) = 0;

		}

		// translation
		Transformation(0, 3) = comR[0] - (Transformation(0, 0)*comT[0] + Transformation(0, 1)*comT[1] + Transformation(0, 2)*comT[2]);
		Transformation(1, 3) = comR[1] - (Transformation(1, 0)*comT[0] + Transformation(1, 1)*comT[1] + Transformation(1, 2)*comT[2]);
		Transformation(2, 3) = comR[2] - (Transformation(2, 0)*comT[0] + Transformation(2, 1)*comT[1] + Transformation(2, 2)*comT[2]);
		Transformation(3, 3) = 1;

		globalTrans = Transformation*globalTrans;

		if (nIters == 0)
			globalTrans1 = globalTrans;
		else if (nIters == 1)
			globalTrans2 = globalTrans;
		else if (nIters == 2)
			globalTrans3 = globalTrans;
		else if (nIters == 3)
			globalTrans4 = globalTrans;

		
		for (int v = 0; v < input->verts.size(); v++)
			{

				
				float x = input->verts[v]->coords[0];
				float y = input->verts[v]->coords[1];
				float z = input->verts[v]->coords[2];
					
				//rotation
				input->verts[v]->coords[0] = Transformation(0, 0)*x + Transformation(0, 1)*y + Transformation(0, 2)*z;
				input->verts[v]->coords[1] = Transformation(1, 0)*x + Transformation(1, 1)*y + Transformation(1, 2)*z;
				input->verts[v]->coords[2] = Transformation(2, 0)*x + Transformation(2, 1)*y + Transformation(2, 2)*z;
				//translation
				input->verts[v]->coords[0] += Transformation(0, 3);
				input->verts[v]->coords[1] += Transformation(1, 3);
				input->verts[v]->coords[2] += Transformation(2, 3);

			}			
		
		mse = 0.0f;
		for (int bv2 = 0; bv2 < input->verts.size(); bv2++)
		{
			mse += pow(reference->verts[matchedVoxels[bv2]]->coords[0] - input->verts[matchedVertices[ matchedVoxels[bv2] ]]->coords[0], 2.0f) +
				pow(reference->verts[matchedVoxels[bv2]]->coords[1] - input->verts[matchedVertices[matchedVoxels[bv2]]]->coords[1], 2.0f) +
				pow(reference->verts[matchedVoxels[bv2]]->coords[2] - input->verts[matchedVertices[matchedVoxels[bv2]]]->coords[2], 2.0f);
		}
		mse /= input->verts.size();
	
		//std::cout << "mse of ICP iter" << nIters-1 << ": " << mse << "\tmse diff\t" << fabs(mse - prevMse) << std::endl;
		if(fabs(mse - prevMse) < closeness){
			break;
		}
		prevMse = mse;

	}
	//std::cout << nIters-1 << "'th ICP iteration w/ final mse: " << mse << "\nICP done!\n\n";

	

	
	ofstream  of;
	of.open("out.xyz");
	for (int i = 0; i < input->verts.size(); i++){
		of << input->verts[i]->coords[0] << " " << input->verts[i]->coords[1] << " " << input->verts[i]->coords[2] << endl;
		
		im->verts[i]->coords[0] = input->verts[i]->coords[0];
		im->verts[i]->coords[1] = input->verts[i]->coords[1];
		im->verts[i]->coords[2] = input->verts[i]->coords[2];
	}
	of.close();

	delete reference;
	delete input;
	
}