#include <iostream>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>
#include <nvector/nvector_serial.h>


int main()
{

    int N = 5;
    int retval;

    SUNContext ctx;
    sunindextype mu, ml;
    SUNMatrix A;
    double *kthcol;

    mu = 1;
    ml = 1;

    printf("Testing banded matrix\n\n");

    retval = SUNContext_Create(NULL, &ctx);
    A = SUNBandMatrix(N, mu, ml, ctx);

    for (int j = 1; j < N-1; j++)
    {
        kthcol = SUNBandMatrix_Column(A,j);
        kthcol[-1] = 3;
        kthcol[0] = -2;
        kthcol[1] = 1;
    }

    kthcol = SUNBandMatrix_Column(A,0);
    kthcol[0] = -2;
    kthcol[1] = 1;
    kthcol = SUNBandMatrix_Column(A,N-1);
    kthcol[0] = -2;
    kthcol[-1] = 3;

    N_Vector x, y;
    SUNLinearSolver LS;
    x = N_VNew_Serial(N, ctx);
    y = N_VNew_Serial(N, ctx);

    double *xdata, *ydata;

    //  Create a vector of our known values.  ydaa is a pointer to the
    //  Sundials array y
    ydata = N_VGetArrayPointer(y);
    ydata[0] = 1;
    ydata[1] = 2;
    ydata[2] = 3;
    ydata[3] = 4;
    ydata[4] = 5;

    //  Set up our banded solver
    LS = SUNLinSol_Band(y, A, ctx);

    //  Set up the solver.  This isn't technically needed, but the solver
    //  gives the wrong answer otherwise.  I don't know why that's the
    //  case.
    retval = SUNLinSolSetup(LS, A);
    retval = SUNLinSolSolve(LS, A, x, y, 1e-9);

    //  Grab a pointer to the vector x and iterate over it to print the
    //  solution
    xdata = N_VGetArrayPointer(x);
    for (int i = 0; i < N; i++)
        std::cout << xdata[i] << std::endl;

    return 0;
}

