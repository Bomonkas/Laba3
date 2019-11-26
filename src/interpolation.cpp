#include "interpolation.h"

double  **get_uniform_grid(double a, double b, int n)
{
    double **grid = new double* [2];
    ofstream out("uniform.txt");
    grid[0] = new double[n];
    grid[1] = new double[n];
    for (int i = 0; i < n; i++)
    {
        grid[0][i] = a + (b - a) / (n - 1) * i;
        out << grid[0][i] << " ";
    }
    out << endl;
    for (int i = 0; i < n; i++)
    {
        grid[1][i] = FUNC(grid[0][i]);
        out << grid[1][i] << " ";
    }
    out << endl;
    out.close();
    return (grid);
}

double **get_chebysh_grid(double a, double b, int n)
{
    double **grid = new double* [2];
    ofstream out("chebysh.txt");
    grid[0] = new double [n];
    grid[1] = new double [n];
    for (int i = 0; i < n; i++)
    {
        grid[0][n - i - 1] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n)));
        out << grid[0][n - i - 1] << " ";
    }
    out << endl;
    for (int i = 0; i < n; i++)
    {
        grid[1][i] = FUNC(grid[0][i]);
        out << grid[1][i] << " ";
    }
    out << endl;
    out.close();
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

double **get_lagr_uniform_points(double ** grid, double a, double b, int n, int h)
{
    double **lagr_points;
    fstream out("points.txt");

    lagr_points = new double*[2];
    lagr_points[0] = new double[n * h - 1];
    lagr_points[1] = new double[n * h - 1];
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

double **get_lagr_chebysh_points(double ** grid, double a, double b, int n, int h)
{
    double **cheb_points;
    fstream out("points.txt");

    cheb_points = new double*[2];
    cheb_points[0] = new double[n * h - 1];
    cheb_points[1] = new double[n * h - 1];
    for (int i = 0; i < n * h - 1; i++)
    {
        cheb_points[0][n * h - i - 2] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n)));
        out << cheb_points[0][n * h - i - 2] << " ";
    }
    out << endl;
    for (int j = 0; j < n * h - 2; j++)
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

void    spline_inter(double **grid, const int m)
{
    ofstream fout("spline.txt");
    fout << m << endl;
    int n = m - 1;
    double *a = new double[n + 1];
    double *b = new double[n + 1];
    double *c = new double[n + 1];
    double *d = new double[n + 1];

    double *aa = new double[n];
    double *bb = new double[n];
    double *cc = new double[n];
    double *dd = new double[n];

    double *h = new double[n + 1];
    double *g = new double[n + 1];
    
    for (int i = 1; i <= n; i++)
    {
        a[i] = grid[1][i - 1];
        h[i] = grid[0][i] - grid[0][i - 1];
        g[i] = (grid[1][i] - grid[1][i - 1]) / h[i];
    }
    for (int i = 2; i <= n; i++)
    {
        aa[i - 1] = h[i - 1];
        bb[i - 1] = 2 * (h[i - 1] + h[i]);
        cc[i - 1] = h[i];
        dd[i - 1] = 3 * (g[i] - g[i - 1]);
    }
    c = sweep_method(aa, bb, cc, dd, n);
    for (int i = 1; i <= n; i++)
    {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    // cout << "\t" <<"a" << "\t|\t" << "b" << "\t|\t" << "c" << "\t|\t" << "d" //
    //     << "\t|\t" << "(x[i-1], x[i])\n" //
    //     << "----------------------------------------------------------------------------------------\n";
    for (int i = 1; i <= n; i++)
    {
        put_zero(&a[i]);
        put_zero(&b[i]);
        put_zero(&c[i]);
        put_zero(&d[i]);
        // cout << setw(10) << a[i] << "\t|" << setw(10) << b[i] << "\t|" << setw(10) //
        //     << c[i] << "\t|" << setw(10) << d[i] << "\t|" << setw(10) << grid[0][i-1] //
        //     << " < x < " << grid[0][i] << endl;
        fout << a[i] << " " << b[i] << " " << c[i] << " " << d[i] << endl;
    }
    fout.close();
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] aa;
    delete[] bb;
    delete[] cc;
    delete[] dd;
    delete[] h;
    delete[] g;
}

double  *sweep_method(double *a, double *b, double *c, double *d, const int n)
{
    double *alpha = new double[n + 1];
    double *beta = new double[n + 1];
    double y = 0.;
    double *x = new double[n + 1];

    alpha[2] = c[1] / b[1];
    beta[2] = d[1] / b[1];
    for (int i = 2; i <= n - 1; i++)
    {
        y = b[i] - a[i] * alpha[i];
        alpha[i + 1] = c[i] / y;
        beta[i + 1] = (d[i] + a[i] * beta[i]) / y;
    }
    
    x[n] = (d[n-1] + a[n-1] * beta[n-1]) / (b[n-1] - a[n-1] * alpha[n-1]);
    for (int i = n - 1; i >= 1; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    delete[] alpha;
    delete[] beta;
    
    return (x);
}
