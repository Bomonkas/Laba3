#include "interpolation.h"

int     main()
{
    double  a = LEFT;
    double  b = RIGHT;
    int     n = NUM;

    double **uniform_grid = get_uniform_grid(a, b, n);
    double **chebysh_grid = get_chebysh_grid(a, b, n);
    get_lagr_points(uniform_grid, a, b, n, H);
    get_chebysh_points(chebysh_grid, a, b, n, H);
    cout << "Uniform grid : " << endl;
    print_grid(uniform_grid, n);
    cout << "Chebysh grid : " << endl;
    print_grid(chebysh_grid, n);

    lagrang_inter(uniform_grid, 0.5, n);
    spline_inter(uniform_grid, n);
    return (0);
}