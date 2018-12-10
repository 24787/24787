#include "Laplace.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned int uint;

Laplace::Laplace()
{
}

void Laplace::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch(ng)
	{
	case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.1127016653792583;			Gpt[1]=0.5;			Gpt[2]=0.8872983346207417;
			wght[0]=0.5555555555555556;			wght[1]=0.8888888888888889;			wght[2]=0.5555555555555556;
			break;
		}
	case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.06943184420297371;			Gpt[1]=0.33000947820757187;			Gpt[2]=0.6699905217924281;			Gpt[3]=0.9305681557970262;
			wght[0]=0.3478548451374539;			wght[1]=0.6521451548625461;			wght[2]=0.6521451548625461;			wght[3]=0.3478548451374539;
			break;
		}
	case 5:
	{
			  Gpt.resize(ng);
			  wght.resize(ng);
			  Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;			Gpt[2] = 0.5;			Gpt[3] = 0.7692346550528415;  Gpt[4] = 0.953089922969332;
			  wght[0] = 0.2369268850561891;			wght[1] = 0.4786286704993665;			wght[2] = 0.5688888888888889;			wght[3] = 0.4786286704993665; wght[4] = 0.2369268850561891;
			  break;
	}
	default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	}
}

void Laplace::InitializeProblem(const vector<double>& CA0, const vector<double>& var, double tstep)
{
	CA.resize(CA0.size());
	CA = CA0;
	if (var.size() != 0) DA = var[0];
	else
	{
		cerr << "0 variables!\n"; getchar();
	}
	dt = tstep;
}

void Laplace::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ)
{
	double Nu[2] = { (1. - u), u };
	double Nv[2] = { (1. - v), v };
	double dNdu[2] = { -1., 1. };
	double dNdv[2] = { -1., 1. };
	double dNdt[4][2];
	double X(0), Y(0), Z(0);

	int i, j, k, a, b, loc;
	loc = 0;

	for (j = 0; j < 2; j++)
	{
		X += pt[loc][0] * Nu[j] * Nv[0];
		Y += pt[loc][1] * Nu[j] * Nv[0];
		loc++;
	}
	for (j = 1; j >= 0; j--)
	{
		X += pt[loc][0] * Nu[j] * Nv[1];
		Y += pt[loc][1] * Nu[j] * Nv[1];
		loc++;
	}

	loc = 0;
	for (k = 0; k < 2; k++)
	{
		Nx[loc] = Nu[k] * Nv[0];
		dNdt[loc][0] = dNdu[k] * Nv[0];
		dNdt[loc][1] = Nu[k] * dNdv[0];
		loc++;
	}
	for (k = 1; k >= 0; k--)
	{
		Nx[loc] = Nu[k] * Nv[1];
		dNdt[loc][0] = dNdu[k] * Nv[1];
		dNdt[loc][1] = Nu[k] * dNdv[1];
		loc++;
	}

	Matrix2d dxdt = Matrix2d::Zero();
	loc = 0;
	for (i = 0; i<2; i++)
	{
		for (j = 0; j<2; j++)
		{
			for (a = 0; a<2; a++)
			{
				for (b = 0; b<2; b++)
				{
					dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
				}
			}
			loc++;
		}
	}

	Matrix2d dtdx = dxdt.inverse();
	for (i = 0; i < 4; i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0);
		dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1);
	}
	detJ = dxdt.determinant();
	detJ = 0.25*detJ;
}

void Laplace::ElementMassMatrix(double Nx[4], double detJ, double EK[4][4])
{
	int i, j;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			EK[i][j] += Nx[i] * Nx[j] * detJ;
		}
	}
}

void Laplace::ElementStiffMatrix(double dNdx[8][3], double detJ, double EK[8][8])
{
	int i, j;
	for (i = 0; i<8; i++)
	{
		for (j = 0; j<8; j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::Assembly(double EM[8][8], double EK[8][8], const vector<int>& IEN, SparseMatrix<double>& GM, SparseMatrix<double>& GK)
{
	unsigned int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			GM.coeffRef(A, B) += EM[i][j];
			GK.coeffRef(A, B) += EK[i][j];
		}
	}
}

void Laplace::BuildLinearSystem(const vector<Element3D>& mesh, SparseMatrix<double>& GM, SparseMatrix<double>& GK)
{
	unsigned int e, i, j, k, a, b;
	double EM[8][8], EK[8][8];
	double Nx[8];
	double dNdx[8][3];
	double detJ;
	for (e = 0; e<mesh.size(); e++)
	{
		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				EM[i][j] = 0.;
				EK[i][j] = 0.;
			}
		}
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], mesh[e].pts, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					ElementMassMatrix(Nx, detJ, EM);
					ElementStiffMatrix(dNdx, detJ, EK);
				}
			}
		}
		for (i = 0; i<8; i++)
		{
			for (j = 0; j<8; j++)
			{
				EK[i][j] = EM[i][j] - dt*DA*EK[i][j];
			}
		}
		Assembly(EM, EK, mesh[e].IEN, GM, GK);
	}
}

void Laplace::Solver(SparseMatrix<double>& GM, SparseMatrix<double>& GK, const vector<double>& Bv)
{
	VectorXd tmp(CA.size());
	for (uint i = 0; i < CA.size(); i++) tmp[i] = CA[i];
	VectorXd GF = GK*tmp;
	for (uint i = 0; i < CA.size(); i++) GF[i] += dt*Bv[i];

	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GM).solve(GF);
	for(uint i=0; i<CA.size(); i++)
	{
		CA[i] = sol[i];
	}
}

void Laplace::VisualizeVTK(const vector<array<double, 3>>& spt, const vector<Element3D>& mesh, string fn)
{
	string fname = fn + "_CA.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 9 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "8 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
				<< " " << mesh[i].IEN[4] << " " << mesh[i].IEN[5] << " " << mesh[i].IEN[6] << " " << mesh[i].IEN[7] << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "12\n";
		}
		fout<<"POINT_DATA "<<CA.size()<<"\nSCALARS CA float 1\nLOOKUP_TABLE default\n";
		for(uint i=0; i<CA.size(); i++)
		{
			fout<<CA[i]<<"\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void Laplace::Run(const vector<array<double, 3>>& pts, const vector<Element3D>& bzmesh, const vector<double>& Bv, vector<double>& CA1, string fn)
{
	cout << "Diffusion...\n";
	GaussInfo();
	SparseMatrix<double> GM(CA.size(), CA.size());
	SparseMatrix<double> GK(CA.size(), CA.size());
	GM.setZero();
	GK.setZero();
	cout<<"Building linear system...\n";
	BuildLinearSystem(bzmesh, GM, GK);
	Solver(GM, GK, Bv);
	cout << "Done diffusion!\n";

	CA1.clear();
	CA1.resize(CA.size());
	CA1 = CA;
}