#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "solver.h"

int main(int argc, char *argv[]) {
    char *path_ils;
    char *path_newton;

    if (argc == 3) {
        path_ils = argv[1];
        path_newton = argv[2];
    } else {
        path_ils = "../solution_ils.dat";
        path_newton = "../solution_newton.dat";
    }

    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double dx = 1.0 / (N - 1);
    double *u = (double *)malloc((N + 1) * sizeof(double));

    if (!u) {
        fprintf(stderr, "Erreur d'allocation mémoire pour u\n");
        MPI_Finalize();
        return 0;
    }

    if (rank == 0) {
        printf(
            "===============================================\n"
            "   Solver équation de diffusion non linéaire   \n"
            "                    Flamme                     \n"
            "===============================================\n\n");

        printf(
            "Paramètre\n"
            "=====================\n"
            "Dimension: %d\n"
            "Itération Max: %d\n"
            "Tolérance: %.e\n"
            "kappa0: %.4f\n"
            "sigma: %.4f\n"
            "beta: %.4f\n"
            "gamma: %.4f\n"
            "=====================\n\n",
            N, MAX_ITER, EPSILON, KAPPA0, SIGMA, BETA, GAMMA);

        printf(
            "Fichier de sauvegarde\n"
            "Schéma implicite linéarisé: %s\n"
            "Méthode de Newton: %s\n\n",
            path_ils, path_newton);
    }

    if (rank == 0) {
        printf(
            "Résolution avec le schéma implicite linéarisé\n"
            "========================\n");
        printf("Initialisation de la solution ...\n");
    }
    init_solution(u, N);

    if (rank == 0) {
        printf("Résolution avec schéma implicite linéarisé ...\n");
    }
    solve_diffusion_implicit_linear(rank, u, N, dx, MAX_ITER);

    if (rank == 0) {
        printf(
            "Sauvegarde de la solution dans %s ...\n"
            "========================\n\n",
            path_ils);
    }
    save_results(path_ils, u, N);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf(
            "Résolution avec la méthode de Newton\n"
            "========================\n");
        printf("Initialisation de la solution ...\n");
    }
    init_solution(u, N);

    if (rank == 0) {
        printf("Résolution avec la méthode de Newton ...\n");
    }
    solve_diffusion_newton(rank, u, N, dx, MAX_ITER);

    if (rank == 0) {
        printf(
            "Sauvegarde de la solution dans %s ...\n"
            "========================\n",
            path_newton);
    }
    save_results(path_newton, u, N);

    MPI_Barrier(MPI_COMM_WORLD);

    free(u);
    MPI_Finalize();

    return 0;
}
