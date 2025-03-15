#pragma once

double kappa(double u);
void init_solution(double *u, int n);
void solve_diffusion_implicit_linear(int rank, double *u, int n, double dx, int max_iter);
void solve_diffusion_newton(int rank, double *u, int n, double dx, int max_iter);
void save_results(const char *filename, double *u, int n);
