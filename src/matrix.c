#include "matrix.h"

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <stdio.h>
#include <stdlib.h>

void tridiag(double *A_diag, double *A_upper, double *A_lower, double *B, double *u, int n) {
    HYPRE_IJMatrix A;
    HYPRE_IJVector x, y;
    HYPRE_Solver solver;

    // Création de la matrice HYPRE
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, n - 1, 0, n - 1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    // Remplissage avec les valeurs de A_diag, A_upper, et A_lower
    for (int i = 0; i < n; i++) {
        HYPRE_BigInt row = i;
        int num_cols = 0;
        HYPRE_BigInt cols[3];
        double values[3];

        // Diagonale
        cols[num_cols] = i;
        values[num_cols] = A_diag[i];
        num_cols++;

        // Sup
        if (i < n - 1) {
            cols[num_cols] = i + 1;
            values[num_cols] = A_upper[i];
            num_cols++;
        }

        // Inf
        if (i > 0) {
            cols[num_cols] = i - 1;
            values[num_cols] = A_lower[i];
            num_cols++;
        }

        HYPRE_IJMatrixSetValues(A, 1, &num_cols, &row, cols, values);
    }

    HYPRE_IJMatrixAssemble(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n - 1, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n - 1, &y);
    HYPRE_IJVectorSetObjectType(y, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(y);

    for (int i = 0; i < n; i++) {
        HYPRE_BigInt index = i;
        double value = B[i];
        HYPRE_IJVectorSetValues(x, 1, &index, &value);
    }

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorAssemble(y);

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);

    HYPRE_ParCSRGMRESSetMaxIter(solver, 1000);
    HYPRE_ParCSRGMRESSetTol(solver, 1e-6);
    // HYPRE_ParCSRGMRESSetPrintLevel(solver, 2);

    HYPRE_ParCSRMatrix parcsr_matrix;
    HYPRE_IJMatrixGetObject(A, (void**)&parcsr_matrix);

    HYPRE_ParVector par_rhs;
    HYPRE_IJVectorGetObject(x, (void**)&par_rhs);

    HYPRE_ParVector par_solution;
    HYPRE_IJVectorGetObject(y, (void**)&par_solution);

    HYPRE_ParCSRGMRESSetup(solver, parcsr_matrix, par_rhs, par_solution);
    HYPRE_ParCSRGMRESSolve(solver, parcsr_matrix, par_rhs, par_solution);

    for (int i = 0; i < n; i++) {
        HYPRE_BigInt index = i;
        HYPRE_IJVectorGetValues(y, 1, &index, &u[i]);
    }

    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_IJVectorDestroy(x);
    HYPRE_IJVectorDestroy(y);
    HYPRE_IJMatrixDestroy(A);
}
