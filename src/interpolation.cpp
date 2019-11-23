#include "interpolation.h"

double   *get_uniform_grid()
{
    int     n = NUM;
    double *grid = new double[n];
    int     i = 0;

    while (n-- >= 0)
    {
        grid[i] = LEFT + (LEFT - RIGHT) / NUM * i;
        i++;
    }
    return (grid);
}

void      print_uniform_grid(double *grid)
{
    for (int i = 0; i < NUM; i++)
        cout << setw(3) <<  LEFT + (LEFT - RIGHT) / NUM * i << " ";
    cout << endl;
    for (int i = 0; i < NUM; i++)
        cout << setw(3) <<  grid[i];
    cout << endl;
}
