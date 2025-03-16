#include "solver.h"

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "matrix.h"

double kappa(double u) {
    return u * u;
}

void init_solution(double *u, int n) {
    for (int i = 0; i <= n; i++) {
        u[i] = 1.0;
    }
}

void solve_diffusion_implicit_linear(int rank, double *u, int n, double dx, int max_iter) {
    double *Fn = malloc((n + 1) * sizeof(double));
    double *A_diag = malloc((n + 1) * sizeof(double));
    double *A_upper = malloc((n + 1) * sizeof(double));
    double *A_lower = malloc((n + 1) * sizeof(double));
    double *Un1 = malloc((n + 1) * sizeof(double));
    double *Kn12 = malloc(n * sizeof(double));

    if (!Fn || !A_diag || !A_upper || !A_lower || !Un1 || !Kn12) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        exit(EXIT_FAILURE);
    }

    int iter = 0;

    for (iter = 0; iter < max_iter; iter++) {
        // calcul du dt
        double Umax = u[0];
        for (int i = 1; i <= n; i++) {
            if (u[i] > Umax) {
                Umax = u[i];
            }
        }

        double Umax_sq = kappa(Umax);
        double dx_sq = kappa(dx);
        double dt = GAMMA * 2.0 / (4 * SIGMA * (Umax_sq * Umax) + 4 * KAPPA0 * Umax_sq / dx_sq);

        // % calcul du 2nd membre
        for (int i = 1; i < n; i++) {
            double Q = (i * dx < 0.2) ? BETA : 0.0;
            Fn[i] = u[i] + dt * (Q + SIGMA);
        }

        // C.L.
        Fn[0] = u[0] + dt * (BETA + SIGMA);
        Fn[n] = 1.0;

        // calcul de la matrice 3Diag
        for (int i = 1; i < n; i++) {
            Kn12[i - 1] = 0.5 * (kappa(u[i - 1]) + kappa(u[i]));
        }

        A_diag[0] = 1 + dt * (KAPPA0 / dx_sq * (2 * Kn12[0]) + SIGMA * pow(u[0], 3));
        A_upper[0] = -dt * KAPPA0 / dx_sq * (2 * Kn12[0]);
        A_lower[0] = 0;

        for (int i = 1; i < n - 1; i++) {
            A_diag[i] = 1 + dt * (KAPPA0 / dx_sq * (Kn12[i - 1] + Kn12[i]) + SIGMA * pow(u[i], 3));
            A_upper[i] = -dt * KAPPA0 / dx_sq * Kn12[i];
            A_lower[i] = -dt * KAPPA0 / dx_sq * Kn12[i - 1];
        }

        A_diag[n - 1] = 1;
        A_upper[n - 1] = 0;
        A_lower[n - 2] = 0;

        // resolution
        tridiag(A_diag, A_upper, A_lower, Fn, Un1, n);

        // u[0] = u[1];
        // u[n] = 1.0;

        // erreur
        double err = 0.0;
        for (int i = 0; i <= n; i++) {
            err += (Un1[i] - u[i]) * (Un1[i] - u[i]);
        }
        err = sqrt(err) / n;

        // fin
        for (int i = 0; i <= n; i++) {
            u[i] = Un1[i];
        }

        // printf("Itération %d : Umax = %f, dt = %e, erreur = %e\n", iter, Umax, dt, err);

        if (err < EPSILON) {
            if (rank == 0) {
                printf("Convergence atteinte en %d itérations\n", iter);
            }
            break;
        }
    }

    if (iter == max_iter && rank == 0) {
        printf("Convergence non atteinte\n");
    }

    free(Fn);
    free(A_diag);
    free(A_upper);
    free(A_lower);
    free(Un1);
    free(Kn12);
}

void solve_diffusion_newton(int rank, double *u, int n, double dx, int max_iter) {
    double *F = malloc(n * sizeof(double));
    double *J_diag = malloc(n * sizeof(double));
    double *J_upper = malloc(n * sizeof(double));
    double *J_lower = malloc(n * sizeof(double));
    double *du = malloc(n * sizeof(double));
    double *Kn12 = malloc((n - 1) * sizeof(double));
    double *Q = malloc(n * sizeof(double));

    if (!F || !J_diag || !J_upper || !J_lower || !du || !Kn12 || !Q) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        exit(EXIT_FAILURE);
    }

    double delta = 0.2;

    for (int i = 0; i < n; i++) {
        double x = i * dx;
        Q[i] = (x < delta) ? BETA : 0;
    }

    double dx_sq = kappa(dx);
    int iter = 0;

    for (iter = 0; iter < max_iter; iter++) {
        // calcul du 2nd membre -Fi
        // nds internes
        for (int i = 1; i < n; i++) {
            Kn12[i - 1] = 0.5 * (kappa(u[i - 1]) + kappa(u[i]));
        }

        for (int i = 1; i < n - 1; i++) {
            F[i] = KAPPA0 / dx_sq * (Kn12[i - 1] * (u[i - 1] - u[i]) + Kn12[i] * (u[i + 1] - u[i])) - SIGMA * (pow(u[i], 4) - 1) + Q[i];
        }

        // C.L en 0
        F[0] = KAPPA0 / dx_sq * (Kn12[0] * (u[1] - u[0]) + Kn12[0] * (u[1] - u[0])) - SIGMA * (pow(u[0], 4) - 1) + Q[0];
        F[n] = 0.0;

        // calcul de la matrice Jacobienne
        J_diag[0] = KAPPA0 / dx_sq * (2 * Kn12[0]) + 4 * SIGMA * pow(u[0], 3);
        J_upper[0] = -KAPPA0 / dx_sq * (2 * Kn12[0]);
        J_lower[0] = 0;

        for (int i = 1; i < n - 1; i++) {
            J_diag[i] = KAPPA0 / dx_sq * (Kn12[i - 1] + Kn12[i]) + 4 * SIGMA * pow(u[i], 3);
            J_upper[i] = -KAPPA0 / dx_sq * (Kn12[i]);
            J_lower[i - 1] = -KAPPA0 / dx_sq * (Kn12[i - 1]);
        }

        J_diag[n - 1] = 1;
        J_upper[n - 1] = 0;
        J_lower[n - 2] = 0;

        // resolution
        tridiag(J_diag, J_upper, J_lower, F, du, n);

        // erreur
        double norm_F = 0.0;
        for (int i = 0; i < n; i++) {
            norm_F += F[i] * F[i];
        }

        double err = sqrt(norm_F / n);

        // fin
        for (int i = 0; i <= n; i++) {
            u[i] += du[i];
        }

        // printf("Itération %d : erreur = %e\n", iter, err);

        if (err < EPSILON) {
            if (rank == 0) {
                printf("Convergence atteinte en %d itérations\n", iter);
            }
            break;
        }
    }

    if (iter == max_iter && rank == 0) {
        printf("Convergence non atteinte\n");
    }

    free(F);
    free(J_diag);
    free(J_upper);
    free(J_lower);
    free(du);
    free(Kn12);
    free(Q);
}

void save_results(const char *filename, double *u, int n) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Erreur lors de l'ouverture du fichier");
        return;
    }
    for (int i = 0; i <= n; i++) {
        fprintf(file, "%lf\n", u[i]);
    }
    fclose(file);
}
