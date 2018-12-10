#include "diff_react.h"

void SetVariables(vector<double>& var)
{
	double var0[5] = {0.2,1.,0.2,10.,0.033};//DA, DX, DB, k_on, k_off
	var.clear();
	var.resize(5);
	var.assign(var0, var0 + 5);
}

void ReadMesh2D(string fn, vector<array<double, 2>>& pts, vector<Element2D>& mesh, vector<int>& pid_loc)//need vtk file with point label
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3];
			for (int j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]][0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]][1];
			}

		}
		fin.close();
		cout << "Mesh Loaded!" << endl;
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

//void ReadMesh(string fn, vector<array<double, 3>>& pts, vector<Element3D>& mesh, vector<array<double, 3>>& pts_b, vector<Element2D>& mesh_b, vector<int>& pid_loc)
//{
//	string fname(fn), stmp;
//	int npts, neles, itmp;
//	ifstream fin;
//	fin.open(fname);
//	if (fin.is_open())
//	{
//		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
//		fin >> stmp >> npts >> stmp;
//		pts.resize(npts);
//		for (int i = 0; i<npts; i++)
//		{
//			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
//		}
//		getline(fin, stmp);
//		fin >> stmp >> neles >> itmp;
//		mesh.resize(neles);
//		for (int i = 0; i<neles; i++)
//		{
//			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3] >>
//				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7];
//			for (int j = 0; j < 8; j++)
//			{
//				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]][0];
//				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]][1];
//				mesh[i].pts[j][2] = pts[mesh[i].IEN[j]][2];
//			}
//		}
//		fin.close();
//	}
//	else
//	{
//		cerr << "Cannot open " << fname << "!\n";
//	}
//
//	vector<array<int, 4>> face;
//	vector<array<int, 4>> fc_sort;
//	vector<array<int, 6>> mesh_fc(mesh.size());
//	int fcloc[6][4] = { { 0, 3, 2, 1 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 0, 4, 7, 3 }, { 4, 5, 6, 7 } };
//	for (unsigned int i = 0; i < mesh.size(); i++)
//	{
//		for (int j = 0; j < 6; j++)
//		{
//			array<int, 4> tmp1 = { mesh[i].IEN[fcloc[j][0]], mesh[i].IEN[fcloc[j][1]], mesh[i].IEN[fcloc[j][2]], mesh[i].IEN[fcloc[j][3]] };
//			array<int, 4> tmp2(tmp1);
//			sort(tmp2.begin(), tmp2.end());
//			vector<array<int, 4>>::iterator it = find(fc_sort.begin(), fc_sort.end(), tmp2);
//			mesh_fc[i][j] = it - fc_sort.begin();
//			if (it == fc_sort.end())
//			{
//				face.push_back(tmp1);
//				fc_sort.push_back(tmp2);
//			}
//		}
//	}
//	vector<vector<int>> fc2mesh(face.size());
//	for (unsigned int i = 0; i < mesh.size(); i++)
//	{
//		for (int j = 0; j < 6; j++)
//		{
//			fc2mesh[mesh_fc[i][j]].push_back(i);
//		}
//	}
//	vector<int> pt_flag(pts.size(), 0);
//	for (unsigned int i = 0; i < face.size(); i++)
//	{
//		if (fc2mesh[i].size() == 1)
//		{
//			Element2D etmp;
//			for (int j = 0; j < 4; j++)
//			{
//				pt_flag[face[i][j]] = 1;
//				etmp.IEN[j] = face[i][j];
//				etmp.pts[j][0] = pts[face[i][j]][0];
//				etmp.pts[j][1] = pts[face[i][j]][1];
//				etmp.pts[j][2] = pts[face[i][j]][2];
//			}
//			mesh_b.push_back(etmp);
//		}
//	}
//	int count(0);
//	vector<int> pt_loc(pts.size(),-1);
//	for (unsigned int i = 0; i < pts.size(); i++)
//	{
//		if (pt_flag[i] == 1)
//		{
//			pts_b.push_back(pts[i]);
//			pt_loc[i] = count++;
//		}
//	}
//	for (unsigned int i = 0; i < mesh_b.size(); i++)
//	{
//		for (int j = 0; j < 4; j++)
//		{
//			if (pt_loc[mesh_b[i].IEN[j]] == -1)
//			{
//				cerr << "wrong index!\n"; getchar();
//			}
//			mesh_b[i].cnct[j] = pt_loc[mesh_b[i].IEN[j]];
//		}
//	}
//	pid_loc = pt_loc;
//
//	//cout << "# pts_b: " << pts_b.size() << "\n";
//	//cout << "# mesh_b: " << mesh_b.size() << "\n";
//	//getchar();
//}

void SetInitialCondition(int ndof, int ndof_b, vector<double>& CA0, vector<double>& CA0_b, vector<double>& NX0, vector<double>& NB0)
{
	double val_ini[3] = {100.,1.e-6,0.};//CA, NX, NB
	CA0.clear();
	CA0_b.clear();
	NX0.clear();
	NB0.clear();
	CA0.resize(ndof, val_ini[0]);
	CA0_b.resize(ndof_b, val_ini[0]);
	NX0.resize(ndof_b, val_ini[1]);
	NB0.resize(ndof_b, val_ini[2]);
}
