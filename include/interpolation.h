#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#define FUNC(x) exp(x)
#define NUM 10
#define LEFT -2
#define RIGHT 2

double  **get_uniform_grid(double a, double b, int n);
void    print_grid(double **grid, int n);
double  *spline_inter(double **grid, double x, const int n);
double  *sweep_method(double *a, double *b, double *c, double *d, const int n);