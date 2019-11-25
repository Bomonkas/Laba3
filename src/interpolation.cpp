#include "interpolation.h"

double  **get_uniform_grid(double a, double b, int n)
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
double  *spline_inter(double **grid, double x, const int n) {
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *d = new double[n];

    double *aa = new double[n];
    double *bb = new double[n];
    double *cc = new double[n];
    double *dd = new double[n];

    double *h = new double[n];
    double *g = new double[n];

    double *res = new double[n];

    for (int i = 1; i <= n; i++)
    {
        a[i] = grid[1][i - 1];
        h[i] = grid[0][i] - grid[0][i - 1];
        g[i] = (grid[1][i] - grid[1][i - 1]) / h[i];
    }
    for (int i = 2; i <= n; i++)
    {
        aa[i - 2] = h[i - 1];
        bb[i - 2] = 2 * (h[i - 1] + h[i]);
        cc[i - 2] = h[i];
        dd[i - 2] = 3 * (g[i] - g[i - 1]);
    }
    c = sweep_method(aa, bb, cc, dd, n - 2);
    for (int i = 1; i <= n; i++)
    {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    for (int i = 1; i <= n; i++)
        res[i - 1] =  a[i] + b[i] * (x - grid[0][i - 1]) + //
        c[i] * (x - grid[0][i - 1]) * (x - grid[0][i - 1]) +//
        d[i] * (x - grid[0][i - 1]) * (x - grid[0][i - 1]) * (x - grid[0][i - 1]);
        
    return (res);
}

double  *sweep_method(double *a, double *b, double *c, double *d, const int n)
{
    double *alpha = new double[n];
    double *beta = new double[n];
    double *y = new double[n];
    double *x = new double[n];

    y[0] = b[0];
    alpha[0] = -c[0]/y[0];
    beta[0] = d[0]/y[0];
    for (int i = 1; i <= n - 1; i++)
    {
        y[i] = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i]/y[i];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
    }
    y[n] = b[n] + a[n] * alpha[n - 1];
    beta[n] = (d[n] - a[n] * beta[n - 1]) / y[n];
    x[n] = beta[n];
    for (int i = n - 1; i >= 1; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    delete[] alpha;
    delete[] beta;
    delete[] y;

    return (x);
}
