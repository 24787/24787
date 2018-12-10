#include <iostream>
#include <vector>
#include <array>
#include "PDE.h"
#include <sstream>
#include <iomanip>
#include <random>
#include <ctime>

using namespace std;

int main()
{
	double tstep(0.01);
	double D_array[5] = { 1.,3.,6.,8.,10. };
	double K_array[5] = { 1.,3.,6.,8.,10. };
	double bc_array[5] = { 0.0, 0.25, 0.67, 0.8, 1.0 };
	double bc_set[4] = { 2. / 3.,0.,0.,0. };
	double D_diff;
	double K_reaction;
	int nstep(1000);
	//nstep=21;
	PDE diffuse_react;
	diffuse_react.ReadMesh2D("../io/ML_project/mesh_21_21.vtk");
		
	
	string fn_time_bc("../io/ML_project/time_bc.txt");
	ofstream fout_time_bc;
	fout_time_bc.open(fn_time_bc.c_str(),ios::app);

	for (int aa = 0; aa < 4000; aa++)
	{
		random_device rd;
		mt19937 mt(rd());
		uniform_real_distribution<double> dist(0.0, 1.0);

		for (int i = 0; i < 4; ++i)
			bc_set[i] = dist(mt);
		D_diff = 1.0;
		K_reaction = 1.0;

		stringstream ss_fn;
		ss_fn << "mesh_" << bc_set[0] << "_" << bc_set[1] << "_" << bc_set[2] << "_" << bc_set[3] << "_";
		string fn_para("");
		string fn_out1("../io/ML_project/physical_bc/" + ss_fn.str());
		string fn_out2("../io/ML_project/parametric_bc/" + ss_fn.str());
		string fn_out3("../io/ML_project/traindata_bc/" + ss_fn.str());

		clock_t start = clock();
		
		diffuse_react.InitializeProblem(D_diff, K_reaction, tstep, bc_set);

		for (int i = 0; i < nstep; i++)
		{
			/*if ((i != 0 && (i + 1) % 10 == 0))
			{
				cout << "Step: " << i << "\n";
			}*/

			if ((i != 0 && (i + 1) % 10 == 0))
			{
				stringstream ss;
				if (i != 0)
					ss << (i + 1) * tstep;
				else
					ss << i;
				diffuse_react.VisualizeVTK(fn_out1 + ss.str());
				diffuse_react.VisualizeVTKParametric(fn_out2 + ss.str());
				diffuse_react.OutputSolution_ML(fn_out3 + ss.str(), 21, 21);
				//cout << "Step " << i << " done!\n";
				//getchar();
			}
			diffuse_react.Run();
		}
		cout << "BC: " << aa << " Done!\n";
		clock_t end = clock();
		fout_time_bc << ss_fn.str() << "     " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	}
	fout_time_bc.close();
	
	///Scalar feature D, K, t
	string fn_time("../io/ML_project/time_DKt.txt");
	ofstream fout_time;
	fout_time.open(fn_time.c_str(), ios::app);
	for (int aa = 0; aa < 5; aa++)
	{
		for (int bb = 0; bb < 5; bb++)
		{
			D_diff = D_array[aa];
			K_reaction = K_array[bb];

			stringstream ss_fn;
			ss_fn << "mesh_" << D_diff << "_" << K_reaction << "_";
			string fn_para("");
			string fn_out1("../io/ML_project/physical/" + ss_fn.str());
			string fn_out2("../io/ML_project/parametric/" + ss_fn.str());
			string fn_out3("../io/ML_project/traindata/" + ss_fn.str());

			clock_t start = clock();
			diffuse_react.InitializeProblem(D_diff, K_reaction, tstep, bc_set);			

			for (int i = 0; i < nstep; i++)
			{
				if (i == 0 || (i != 0 && (i + 1) % 10 == 0))
				{
					cout << "Step: " << i << "\n";
				}

				if (i == 0 || (i != 0 && (i + 1) % 10 == 0))
				{
					stringstream ss;
					if (i != 0)
						ss << (i + 1) * tstep;
					else
						ss << i;
					diffuse_react.VisualizeVTK(fn_out1 + ss.str());
					diffuse_react.VisualizeVTKParametric(fn_out2 + ss.str());
					diffuse_react.OutputSolution_ML(fn_out3 + ss.str(), 21, 21);
					cout << "Step " << i << " done!\n";
					//getchar();
				}
				diffuse_react.Run();
			}
			cout << "D, K: " << D_diff << "_" << K_reaction << " Done!\n";
			clock_t end = clock();
			fout_time << ss_fn.str() << "     " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
		}
	}
	fout_time.close();


	cout<<"DONE!\n";
	getchar();
	return 0;
}

