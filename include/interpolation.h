#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstring>

using namespace std;

#define FUNC(x) x*x //pow(4*x*x*x+2*x*x-4*x+2,sqrt(2))+asin(1/(5+x-x*x))-5
#define NUM 3
#define LEFT -1
#define RIGHT 1
#define EPS 10e-7
#define H 100

double      **get_uniform_grid(double a, double b, int n);
void        print_grid(double **grid, int );
double      lagrang_inter(double **grid, double x, int n);
double		**get_chebysh_grid(double a, double b, int n);
double		**get_lagr_uniform_points(double ** grid, double a, double b, int n, int h);
double		**get_lagr_chebysh_points(double ** grid, double a, double b, int n, int h);
void		put_zero(double *x);
void	  	spline_inter(string name_file, double **grid, const int n, double a_, double b_);
double  	*sweep_method(double *a, double *b, double *c, double *d, const int n);
void 		print_v(double *vec, const int size);
void		solveMatrix(int n, double *a, double *c, double *b, double *f, double *x);
