#include "interpolation.h"

double   *get_uniform_grid(double a, double b, int n)
{
    double *grid = new double[n + 1];
    double  x;

    for (int i = 0; i <= n; i++)
    {
        x = a + (b - a) / n * i;
        grid[i] = FUNC(x);
    }
    return (grid);
}

void      print_uniform_grid(double *grid, double a, double b, int n)
{
    cout <<"uniform grid :\n" << " x =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) <<  a + (b - a) / n * i<< " | ";
    cout << endl << " y =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) << setprecision(4) <<  grid[i] << " | ";
    cout << endl;
}
