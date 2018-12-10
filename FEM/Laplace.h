#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"

using namespace std;
using namespace Eigen;

class Laplace
{
public:
	Laplace();
	void GaussInfo(int ng=2);
	void InitializeProblem(const vector<double>& CA0, const vector<double>& var, double tstep);
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[8], double dNdx[8][3], double& detJ);
	void ElementMassMatrix(double Nx[8], double detJ, double EM[8][8]);
	void ElementStiffMatrix(double dNdx[8][3], double detJ, double EK[8][8]);
	void Assembly(double EM[8][8], double EK[8][8], const vector<int>& IEN, SparseMatrix<double>& GM, SparseMatrix<double>& GK);
	void BuildLinearSystem(const vector<Element3D>& mesh, SparseMatrix<double>& GM, SparseMatrix<double>& GK);
	void Solver(SparseMatrix<double>& GM, SparseMatrix<double>& GK, const vector<double>& Bv);
	void VisualizeVTK(const vector<array<double, 3>>& pts, const vector<Element3D>& mesh, string fn);
	void Run(const vector<array<double, 3>>& pts, const vector<Element3D>& bzmesh, const vector<double>& Bv, vector<double>& CA1, string fn);

private:
	vector<double> Gpt;
	vector<double> wght;
	vector<double> CA;
	double DA;
	double dt;
	int ndt;
};

#endif