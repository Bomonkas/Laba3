#include "interpolation.h"

int     main()
{
    double  a = LEFT;
    double  b = RIGHT;
    int     n = NUM;

    double *uniform_grid = get_uniform_grid(a, b, n);
    print_uniform_grid(uniform_grid, a, b, n);
    cout << FUNC(4) << endl;
    return (0);
}