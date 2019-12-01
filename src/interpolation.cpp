#include "interpolation.h"

double **get_uniform_grid(double a, double b, int n)
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
        grid[0][n - i - 1] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n)));
    for (int i = 0; i < n; i++)
        out << grid[0][i] << " ";
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

    for (int i = 0; i < n; i++)
    {
        prod = 1.;
        for (int j = 0; j < n; j++)
        {
            if (i != j)
                prod *= (x - grid[0][j]) / (grid[0][i] - grid[0][j]);
        }
        res += prod * grid[1][i]; 
    }
    return (res);
}

void print_grid(double **grid, int n)
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
    fstream out("uni_points.txt");

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
    fstream out("cheb_points.txt");

    cheb_points = new double*[2];
    cheb_points[0] = new double[n * h];
    cheb_points[1] = new double[n * h];
    for (int i = 0; i < n * h; i++)
        cheb_points[0][n * h - 1 - i] = (a + b) /  2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n * h)));
    for (int i = 0; i < n * h; i++)
        out << cheb_points[0][i] << " ";
    out << endl;
    for (int j = 0; j < n * h; j++)
    {
        cheb_points[1][j] = lagrang_inter(grid, cheb_points[0][j], n);
        put_zero(&(cheb_points[1][j]));
        out << cheb_points[1][j] << " ";
    }   
    out << endl;
    out.close();
    return (cheb_points);
}

void put_zero(double *x)
{
    if (abs(*x) < EPS)
        *x = 0.;
}

void print_v(double *vec, const int size)
{
	if (!vec) {
		cout << "Vector is not exist\n";
		return;
	}
	for (int i = 0; i < size; i++) {
		cout << vec[i] << " ";
	}
	cout << endl;
}

void spline_inter(string name_file, double **grid, const int m)
{
    ofstream fout(name_file);
    fout << m << endl;
    int n = m - 1;
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n + 1];
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
        g[i] = (grid[1][i + 1] - grid[1][i]) / h[i];
    }
    for (int i = 0; i < n; i++)
    {
        aa[i] = h[i];
        bb[i] = 2 * (h[i] + h[i + 1]);
        cc[i] = h[i + 1];
        dd[i] = 3 * (g[i + 1] - g[i]);
        put_zero(&aa[i]);
        put_zero(&bb[i]);
        put_zero(&cc[i]);
        put_zero(&dd[i]);
    }

    solveMatrix(n, aa, bb, cc, dd, c);
    // c = sweep_method(aa, bb, cc, dd, n);
    c[0] = 0;
    c[n] = 0;
   
    for (int i = 0; i < n; i++)
    {
        b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
   
    // cout << "\t" <<"a" << "\t|\t" << "b" << "\t|\t" << "c" << "\t|\t" << "d" //
    //     << "\t|\t" << "(x[i-1], x[i])\n" //
    //     << "----------------------------------------------------------------------------------------\n";
    for (int i = 0; i < n; i++)
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

void solveMatrix(int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = 0;
    
	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}

// double *sweep_method(double *a, double *b, double *c, double *d, const int n)
// {
//     double *alpha = new double[n + 1];
//     double *beta = new double[n + 1];
//     double y = 0.;
//     double *x = new double[n + 2];

//     alpha[0] = -c[0] / b[0];
//     beta[0] = d[0] / b[0];
//     for (int i = 1; i < n; i++)
//     {
//         y = b[i] + a[i] * alpha[i - 1];
//         alpha[i] = -c[i] / y;
//         beta[i] = (d[i] - a[i] * beta[i - 1]) / y;
//     }
//     x[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 2]) / (b[n - 1] + a[n - 1] * alpha[n - 2]);
//     for (int i = n - 2; i >= 0; i--)
//     {
//         x[i] = alpha[i] * x[i + 1] + beta[i];
//     }

//     delete[] alpha;
//     delete[] beta;
    
//     return (x);
// }
