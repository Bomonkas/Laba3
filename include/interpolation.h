#pragma once

#include <iostream>
#include <iomanip>

using namespace std;

#define FUNC(x) x*x
#define NUM 10
#define LEFT -1
#define RIGHT 1

double      **get_uniform_grid(double a, double b, int n);
void        print_grid(double **grid, int n);