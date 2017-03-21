#include <string>
using namespace std;

void Static_Solver(double A, double E, double I, double L, double l, double qx, double qy, double Fy, int Nx, int ku, int kl, int N);

void Dynamic_Solver_1(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N);

void Dynamic_Solver_2(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N);

void Dynamic_Solver_MPI_1(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N, int size, int rank);

void Dynamic_Solver_MPI_2(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N, int size, int rank);