#include "interpolation.h"

double   **get_uniform_grid(double a, double b, int n)
{
    double **grid = new double*[2];
    double x;

    grid[0] = new double[n + 1];
    grid[1] = new double[n + 1];
    for (int i = 0; i <= n; i++)
    {
        x = a + (b - a) / n * i;
        grid[0][i] = x;
        grid[1][i] = FUNC(x);
    }
    return (grid);
}

void      print_uniform_grid(double **grid, int n)
{
    cout <<"uniform grid :\n" << " x =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) <<  grid[0][i]<< " | ";
    cout << endl << " y =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) << setprecision(4) <<  grid[1][i] << " | ";
    cout << endl;
}
