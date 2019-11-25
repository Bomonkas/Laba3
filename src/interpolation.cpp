#include "interpolation.h"

double      **get_uniform_grid(double a, double b, int n)
{
    double **grid = new double* [2];
    double x;

    grid[0] = new double[n];
    grid[1] = new double[n];
    for (int i = 0; i < n; i++)
    {
        x = a + (b - a) / (n - 1) * i;
        grid[0][i] = x;
        grid[1][i] = FUNC(x);
    }
    return (grid);
}

double **get_chebysh_grid(double a, double b, int n)
{
    double **grid = new double* [2];
    double x;

    grid[0] = new double [n];
    grid[1] = new double [n];
    for (int i = 0; i < n; i++)
    {
        x = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n)));
        grid[0][n - i - 1] = x;
        grid[1][n - i - 1] = FUNC(x);
    }
    return (grid);
}

double lagrang_inter(double **grid, double x, int n)
{
    double res = 0.;
    double prod;

    for (int i = 0; i <= n; i++)
    {
        prod = 1.;
        for (int j = 0; j <= n; j++)
        {
            if (i != j)
                prod *= (x - grid[0][j]) / (grid[0][i] - grid[0][j]);
        }
        res += prod * grid[1][i]; 
    }
    return (res);
}

void         print_grid(double **grid, int n)
{
    cout << " x = | ";
    for (int i = 0; i < n; i++)
    {
        cout << setw(10) << setprecision(4) <<  grid[0][i] << " | ";

    }
    cout << endl << " y = | ";
    for (int i = 0; i < n; i++)
    {
        cout << setw(10) << setprecision(4) <<  grid[1][i] << " | ";
    }
    cout << endl;
}

double **get_lagr_points(double ** grid, double a, double b, int n, int h)
{
    double **lagr_points;
    fstream out;

    lagr_points = new double*[2];
    lagr_points[0] = new double[n * h - 1];
    lagr_points[1] = new double[n * h - 1];
    out.open("grid.txt", ofstream::out);
    for (int i = 0; i < n * h - 1; i++)
    {
        lagr_points[0][i] = a + ((b - a) / (n - 1) * i / h);
        out << lagr_points[0][i] << " ";
    }
    out << endl;
    for (int j = 0; j < n * h - 1; j++)
    {
        lagr_points[1][j] = lagrang_inter(grid, lagr_points[0][j], n);
        put_zero(&(lagr_points[1][j]));
        out << lagr_points[1][j] << " ";
    }
    out << endl;
    out << endl;
    out.close();
    return (lagr_points);
}

double **get_chebysh_points(double ** grid, double a, double b, int n, int h)
{
    double **cheb_points;
    fstream out;

    cheb_points = new double*[2];
    cheb_points[0] = new double[n * h - 1];
    cheb_points[1] = new double[n * h - 1];
    out.open("grid.txt", ofstream::app);
    for (int i = 0; i < n * h - 1; i++)
    {
        cheb_points[0][n * h - i - 1] = (a + b) / 2 + (b - a) / h /  2 * cos((2 * i + 1) * 3.14 / (2 * (n)));
        out << cheb_points[0][n * h - i - 1] << " ";
    }
    out << endl;
    for (int j = 0; j < n * h - 1; j++)
    {
        cheb_points[1][j] = lagrang_inter(grid, cheb_points[0][j], n);
        put_zero(&(cheb_points[1][j]));
        out << cheb_points[1][j] << " ";
    }
    out << endl;
    out.close();
    return (cheb_points);
}

void     put_zero(double *x)
{
    if (abs(*x) < EPS)
        *x = 0;
}