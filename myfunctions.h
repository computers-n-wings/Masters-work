#include <string>
using namespace std;

void Print_Matrix(double A[], int n, int m);

void Print_Matrix_Parallel(double A[], int n, int m, int size, int rank);

void Print_Vector(double b[], int n);

void Print_Vector_Parallel(double b[], int n, int size, int rank);

void Zero_Vector(double a[], int N);

void Copy_Vector(double a[], double b[], int N);

void Print_Norm(double a[], int N);

void Build_K_elemental(double Ke[], int n, double A, double E, double l, double I);

void Build_K_global_banded(double Kb[], double A, double E, double l, double L, double I, int Nx, int N, int kl, int ku, int lda, int bfr);

void Build_K_global_banded_altered(double Kb[], double Mb[], double A, double E, double l, double L, double I, double del_t, int Nx, int N, int kl, int ku, int lda, int bfr);

void Build_F_global(double F[], double Fy, double qx, double qy, double time, double l, int Nx, int N);

void Build_M_global_banded(double Mb[], double rho, double A, double l, int Nx, int N);

void Banded_Matrix_Solver(double Ab[], double b[], int N,int lda, int kl, int ku);

void Write_Vector(double F[], int N, double l, double L, string mystring);

void Write_Point_Displacement(double F[], int N, string mystring);

void Build_Multiplier1(double F[], double Kb[], double M[], double u0[], double u1[], double del_t, int N, int lda, int kl, int ku);

double RMS_error(double u1[], double S[], double M[], int N); 

void Build_Keff(double Keff[], double Mb[], double coeff1, double A, double E, double l, double L, double I, int Nx, int N, int kl, int ku, int lda, int bfr);

void Build_Multiplier2(double S[], double Mb[], double F[], double u0[], double udot0[], double udotdot0[], double coeff1, double coeff2, double coeff3, int N);

void Build_udotdot(double udotdot1[], double u1[], double u0[], double udot0[], double udotdot0[], double coeff1, double coeff2, double coeff3, int N);

void Build_udot(double udot1[], double udot0[], double udotdot0[], double udotdot1[], double coeff4, double coeff5, int N);

void Build_M_global_banded_MPI(double Mbloc[], double rho, double A, double l, int Nx, int Nxloc, int Nloc, int size, int rank);

void Build_K_global_banded_MPI(double Kbloc[], double A, double E, double l, double L, double I, int Nx, int Nxloc, int N, int kl, int ku, int lda, int bfr, int size, int rank);

void Build_K_global_banded_altered_MPI(double Kbloc[], double Mbloc[], double A, double E, double l, double L, double I, double del_t, int Nx, int Nxloc, int N, int Nloc, int kl, int ku, int lda, int bfr, int size, int rank);

void Build_F_global_MPI(double Floc[], double Fy, double qx, double qy, double time, double l, int Nx, int Nxloc, int N, int Nloc, int size, int rank);

void Build_Multiplier1_MPI(double F[], double Kbloc[], double Mbloc[], double u0loc[], double u1loc[], double del_t, int N, int Nloc, int lda, int kl, int ku, int size, int rank);