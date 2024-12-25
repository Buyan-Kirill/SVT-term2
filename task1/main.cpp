#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

int main(int argc, char *argv[])
{
    // Get number of nodes
    unsigned N = 100;

    // Create sparse matrix, RHS vector and solution vector
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;
    // Set their size
    A.SetInterval(0, N);
    b.SetInterval(0, N);
    sol.SetInterval(0, N);
    // Make A identity matrix
    // Make b: b_i = i
    for(unsigned i = 0; i < N; i++){
        A[i][i] = 1.0;
        b[i] = i;
    }

    // Get solver
    // All inner INMOST solvers are BiCGStab
    // with different preconditioners, let's use ILU2
    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");

    // Set matrix in the solver;
    // this also computes preconditioner
    S.SetMatrix(A);

    // Solve
    bool solved = S.Solve(b, sol);
    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;
    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }
	return 0;
}
