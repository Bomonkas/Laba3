#include "interpolation.h"

double   *get_uniform_grid(double a, double b, int n)
{
    double *grid = new double[n + 1];
    double  x;

    for (int i = 0; i <= n; i++)
    {
        x = a + (b - a) / n * i;
        cout << "x = " << x << " FUNC(x) = " << FUNC(x) << endl;
        grid[i] = FUNC(x);
        n--;
    }
    return (grid);
}

void      print_uniform_grid(double *grid, double a, double b, int n)
{
    for (int i = 0; i <= n; i++)
        cout << setw(6) <<  a + (b - a) / n * i << " | ";
    cout << endl;
    for (int i = 0; i <= n; i++)
        cout << setw(6) << setprecision(4) <<  grid[i] << " | ";
    cout << endl;
}
