#include "interpolation.h"

int     main()
{
    double *uniform_grid = get_uniform_grid();
    print_uniform_grid(uniform_grid);
    cout << FUNC(4) << endl;
    return (0);
}