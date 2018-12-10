#ifndef DIFF_REACT_H
#define DIFF_REACT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
using namespace std;

//problem setting

void SetVariables(vector<double>& var);

void ReadMesh2D(string fn, vector<array<double, 2>>& pts, vector<Element2D>& mesh, vector<int>& pid_loc);

void SetInitialCondition(int ndof, int ndof_b, vector<double>& CA0, vector<double>& CA0_b, vector<double>& NX0, vector<double>& NB0);


#endif