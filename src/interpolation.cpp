#include "interpolation.h"

double  **get_uniform_grid(double a, double b, int n)
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

void    print_grid(double **grid, int n)
{
    cout <<"uniform grid :\n" << " x =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) <<  grid[0][i]<< " | ";
    cout << endl << " y =";
    for (int i = 0; i <= n; i++)
        cout << setw(6) << setprecision(4) <<  grid[1][i] << " | ";
    cout << endl;
}

double  spline_inter(double **grid, double x, const int n) {
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

    for (int i = 0; i < n; i++)
    {
        a[i] = grid[1][i];
        h[i] = grid[0][i + 1] - grid[0][i];
        g[i] = (grid[1][i + 1] - grid[1][i])/h[i];
        c[0] = 0.;
        aa[i] = h[i - 1];
        bb[i] = 2 * (h[i - 1] + h[i]);
        cc[i] = h[i];
        dd[i] = 3 * (g[i] - g[i - 1]);
    }
    c = sweep_method(aa, bb, cc, dd, n);
    for (int i = 0; i < n; i++)
    {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    double sum = 0.;
    for (int i = 1; i < n; i++)
        sum += a[i] + b[i] * (x - grid[0][i - 1]) + //
        c[i] * (x - grid[0][i - 1]) * (x - grid[0][i - 1]) +//
        d[i] * (x - grid[0][i - 1]) * (x - grid[0][i - 1]) * (x - grid[0][i - 1]);
        
    return (sum);
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
    for (int i = 1; i < n - 1; i++)
    {
        y[i] = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i]/y[i];
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y[i];
    }
    y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 1];
    beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 1]) / y[n];
    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = alpha[i] * x[i + 1] + beta[i];
    delete[] alpha;
    delete[] beta;
    delete[] y;

    return (x);
}
