#include "interpolation.h"

int     main()
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

    double x = 1.5;
    cout << "Ln(" << x << ") = " << lagrang_inter(uniform_grid, x, n) << ",   f(" << x <<") = " << FUNC(x) << endl;

    spline_inter(uniform_grid, n);

    delete[] uniform_grid[0];
    delete[] chebysh_grid[0];
    delete[] uniform_grid[1];
    delete[] chebysh_grid[1];
    delete[] uniform_grid;
    delete[] chebysh_grid;
    return (0);
}