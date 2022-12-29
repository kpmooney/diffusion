#include <iostream>
#include <iomanip>
#include <math.h>
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>


using namespace std;

int equations(realtype t, N_Vector T, N_Vector Tprime, N_Vector res,
		void *userdata);

int main()
{

    int N = 101;

    int retval;
    SUNContext ctx;
    N_Vector T, Tprime, res, avtol, id;

    void* mem;
    double *Tval, *Tpval, *atval;
    double rtol, t0, tout, tret;
    int iout, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;

    mem = NULL;
    T = Tprime = avtol = id = NULL;
    Tval = Tpval = atval = NULL;

    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);

    /* Allocate N-vectors. */
    T = N_VNew_Serial(N, ctx);
    Tprime = N_VClone(T);
    avtol = N_VClone(T);
    id = N_VClone(T);

    Tval = N_VGetArrayPointer(T);
    Tval[0]  = 0;
    for (int i = 1; i < N-1; i++)
        Tval[i] = 1;
    Tval[N-1]  = 30;

    rtol = 1.0e-4;

    t0 = 0;
    tout = 1.;

    atval = N_VGetArrayPointer(avtol);
    for (int i = 0; i < N; i++)
        atval[i] = 1.0e-8;


    mem = IDACreate(ctx);
    retval = IDASetId(mem, id);

    A = SUNDenseMatrix(N, N, ctx);

    retval = IDAInit(mem, equations, t0, T, Tprime);
    retval = IDASVtolerances(mem, rtol, avtol);
    LS = SUNLinSol_Dense(T, A, ctx);
    retval = IDASetLinearSolver(mem, LS, A);

    //retval = IDAGetConsistentIC(mem, T, Tprime);
    retval = IDACalcIC(mem, IDA_YA_YDP_INIT, 0.0);

    while (1)
    {
        retval = IDASolve(mem, tout, &tret, T, Tprime, IDA_NORMAL);
        if (tout >= 6000)
            break;

        tout += 1;

    }


    for (int i = 0; i < N; i++)
        cout << setprecision(6) << Tval[i] << endl;

    IDAFree(&mem);

    return 0;
}

int equations(realtype t, N_Vector Y, N_Vector Yprime, N_Vector Yres,
		void *userdata)
{

    int N = 101;
    double D = 1.0e-4;
    double Delta = 1.0 / (N - 1);

    double *T, *Tprime, *resval;

    T = N_VGetArrayPointer(Y);
    Tprime = N_VGetArrayPointer(Yprime);
    resval = N_VGetArrayPointer(Yres);

    resval[0] = T[0] - 0;
    resval[N-1] = T[N-1] - 30;
        for (int i = 1; i <= N - 2; i++)
            resval[i] = Tprime[i] - D /( Delta*Delta) * (T[i+1] - 2*T[i] + T[i-1]);


    return 0;
}
