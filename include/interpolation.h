#pragma once

#include <iostream>
#include <iomanip>

using namespace std;

#define FUNC(x) x*x
#define NUM 4
#define LEFT -1
#define RIGHT 1

double  **get_uniform_grid(double a, double b, int n);
void    print_grid(double **grid, int n);
double  spline_inter(double **grid, double x, const int n);
double  *sweep_method(double *a, double *b, double *c, double *d, const int n);