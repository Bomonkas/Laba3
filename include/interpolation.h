#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

#define FUNC(x) x*x
#define NUM 11
#define LEFT -1
#define RIGHT 1
#define EPS 10e-7
#define H 1

double      **get_uniform_grid(double a, double b, int n);
void        print_grid(double **grid, int );
double      lagrang_inter(double **grid, double x, int n);
double		**get_chebysh_grid(double a, double b, int n);
double		**get_lagr_points(double ** grid, double a, double b, int n, int h);
double		**get_chebysh_points(double ** grid, double a, double b, int n, int h);
void		put_zero(double *x);
double  	*spline_inter(double **grid, double x, const int n);
double  	*sweep_method(double *a, double *b, double *c, double *d, const int n);
