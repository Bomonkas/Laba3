#include "interpolation.h"

int     main()
{
    double  a = LEFT;
    double  b = RIGHT;
    int     n = NUM;

    double **uniform_grid = get_uniform_grid(a, b, n);
    print_grid(uniform_grid, n);

    double *x;
    x = spline_inter(uniform_grid, -2, n);
    for (int i = 0; i < n; i++)
        cout << "\t" << setw(7) << x[i];
    cout << endl;
    cout << FUNC(-2) << endl;
    return (0);
}