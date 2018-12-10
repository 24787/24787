#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip> 

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//const int phys_dim(3);

/////////////////////////////////

class Element2D
{
public:
	vector<int> cnct;
	vector<int> IEN;
	vector<array<double, 3>> pts;//tmp

	Element2D();
};

class Element3D
{
public:
	vector<int> IEN;
	vector<array<double, 3>> pts;//tmp

	Element3D();
};


//mesh format conversion

void Raw2Vtk_hex(string fn);

void Rawn2Vtk_hex(string fn);

#endif