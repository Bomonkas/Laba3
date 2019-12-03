#include "interpolation.h"

int main()
{
    double  a = LEFT;
    double  b = RIGHT;
    int     n = NUM;

    double **uniform_grid = get_uniform_grid(a, b, n);
    double **chebysh_grid = get_chebysh_grid(a, b, n);
    get_lagr_uniform_points(uniform_grid, a, b, n, H);
    get_lagr_chebysh_points(chebysh_grid, a, b, n, H);
    cout << "Uniform grid : " << endl;
    print_grid(uniform_grid, n);
    cout << "Chebysh grid : " << endl;
    print_grid(chebysh_grid, n);

    double x = 0.5;
    double res = lagrang_inter(uniform_grid, x, n);
    double f = FUNC(x);
    cout << setprecision(10) << "Ln(" << x << ") = " << res << ",   f(" << x <<") = " << f << endl;
    cout << setprecision(10) << "Ln(" << x << ") - f(" << x <<") = " << fabs(res - f) << endl;
    spline_inter("unispline.txt", uniform_grid, n, a, b);
    spline_inter("chebspline.txt", chebysh_grid, n, a, b);

    delete[] uniform_grid[0];
    delete[] chebysh_grid[0];
    delete[] uniform_grid[1];
    delete[] chebysh_grid[1];
    delete[] uniform_grid;
    delete[] chebysh_grid;
    return (0);
}