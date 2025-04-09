#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <Windows.h>
using namespace std;

const double EPS = 1e-12;

vector<double> F(const vector<double>& x);
vector<vector<double>> J_analytical(const vector<double>& x);
vector<vector<double>> J_numerical(const vector<double>& x, double M);
vector<double> GaussMethod(vector<vector<double>> A, vector<double> b);
void NewtonMethod(const vector<double>& x0, double eps1, double eps2, int max_iter, double M = -1.0);