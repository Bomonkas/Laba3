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

double **get_lagr_uniform_points(double **grid, double a, double b, int n, int h)
{
    double **uni_points;
    fstream out("uni_points.txt");

    uni_points = new double*[2];
    int m = n * h;
    uni_points[0] = new double[m];
    uni_points[1] = new double[m];
    for (int i = 0; i < m; i++)
    {
        uni_points[0][i] = a + (b - a) / (m - 1) * i;
        out << uni_points[0][i] << " ";
    }
    out << endl;
    for (int j = 0; j < m; j++)
    {
        uni_points[1][j] = lagrang_inter(grid, uni_points[0][j], n);
        put_zero(&(uni_points[1][j]));
        out << uni_points[1][j] << " ";
    }
    out << endl;
    out << endl;
    out.close();
    return (uni_points);
}

double **get_lagr_chebysh_points(double **grid, double a, double b, int n, int h)
{
    double **cheb_points;
    fstream out("cheb_points.txt");

    int m = n * h;
    cheb_points = new double*[2];
    cheb_points[0] = new double[m];
    cheb_points[1] = new double[m];
    for (int i = 0; i < m; i++)
        cheb_points[0][m - 1 - i] = (a + b) /  2 + (b - a) / 2 * cos((2 * i + 1) * 3.14 / (2 * m));
    for (int i = 0; i < m; i++)
        out << cheb_points[0][i] << " ";
    out << endl;
    for (int j = 0; j < m; j++)
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

double get_lagr_error(double **grid, double a, double b, int n)
{
    double error = 0.;
    double x = 0.;
    double dif = 0.;
    while (x < b)
    {
        x += (b - a) / 1000;
        dif = fabs(FUNC(x) - lagrang_inter(grid, x, n));
        if (error < dif)
        error = dif;
    }
    return (error);
}

void spline_inter(string name_file, double **grid, const int m, double a_, double b_)
{
    ofstream fout(name_file);
    fout << m << endl;
    int n = m - 1;
    double *a = new double[n];
    double *b = new double[n];
    double *res = new double[n];
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
    aa[0] = 0;
    for (int i = 0; i < n; i++)
    {
        if (i != 0)
        aa[i] = h[i];
        bb[i] = 2 * (h[i] + h[i + 1]);
        if (i != n)
        cc[i] = h[i + 1];
        dd[i] = 3 * (g[i + 1] - g[i]);
        put_zero(&aa[i]);
        put_zero(&bb[i]);
        put_zero(&cc[i]);
        put_zero(&dd[i]);
    }
    cc[n - 1] = 0;

    solveMatrix(n, aa, bb, cc, dd, res);
    double *c = new double[n + 1];
    for (int i = 0; i < n; i++)
        c[i + 1] = res[i]; 
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

    if (name_file == "unispline.txt")
    {
        double **uni_points;
        fstream out("spline_uni_points.txt");

        int l = m * H;
        uni_points = new double*[2];
        uni_points[0] = new double[l];
        uni_points[1] = new double[l];
        for (int i = 0; i < l; i++)
        {
            uni_points[0][i] = a_ + (b_ - a_) / (l - 1) * i;
            out << uni_points[0][i] << " ";
        }
        out << endl;
        int k = 0;
        print_v(grid[0], m);
        for (int j = 0; j < l; j++)
        {
            k = 0;
            while ((uni_points[0][j] > grid[0][k + 1]) && (k + 1 < m - 1))
            k++;
            uni_points[1][j] = a[k] + b[k] * (uni_points[0][j] - grid[0][k]) + //
                                c[k] * pow((uni_points[0][j] - grid[0][k]), 2) + //
                                d[k] * pow((uni_points[0][j] - grid[0][k]), 3);
            put_zero(&(uni_points[1][j]));
            out << uni_points[1][j] << " ";
        }
        out << endl;
        out << endl;
        out.close();
    } else
    {
        double **cheb_points;
        fstream out1("spline_cheb_points.txt");

        int l = m * H;
        cheb_points = new double*[2];
        cheb_points[0] = new double[l];
        cheb_points[1] = new double[l];
        print_v(grid[0], m);
        a_ = grid[0][0];
        b_ = grid[0][m - 1];
        cout << a_ << " " << b_ << endl;
        for (int i = 0; i < l; i++)
            cheb_points[0][l - 1 - i] = (a_ + b_) /  2 + (b_ - a_) / 2 * cos((2 * i + 1) * 3.14 / (2 * l));
        for (int i = 0; i < l; i++)
            out1 << cheb_points[0][i] << " ";
        out1 << endl;
        int k = 0;
        for (int j = 0; j < l; j++)
        {
            k = 0;
            while ((cheb_points[0][j] > grid[0][k + 1]) && (k + 1 < m - 1))
                k++;
            cheb_points[1][j] = a[k] + b[k] * (cheb_points[0][j] - grid[0][k]) + //
                                c[k] * pow((cheb_points[0][j] - grid[0][k]), 2) + //
                                d[k] * pow((cheb_points[0][j] - grid[0][k]), 3);
            put_zero(&(cheb_points[1][j]));
            out1 << cheb_points[1][j] << " ";
        }
        out1 << endl;
        out1 << endl;
        out1.close();
    }

    delete[] a;
    delete[] res;
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

	x[n-1] = f[n-1] / c[n-1];
    
	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}

double *sweep_method(double *a, double *b, double *c, double *d, const int n)
{
    double *alpha = new double[n + 1];
    double *beta = new double[n + 1];
    double y = 0.;
    double *x = new double[n + 2];

    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        y = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / y;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y;
    }
    x[n - 1] = (d[n - 1] + a[n - 1] * beta[n - 2]) / (b[n - 1] + a[n - 1] * alpha[n - 2]);
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    delete[] alpha;
    delete[] beta;
    
    return (x);
}
