#include <string>
using namespace std;

void Print_Matrix(double A[], int n, int m);

void Print_Vector(double b[], int n);

void Zero_Vector(double a[], int N);

void Copy_Vector(double a[], double b[], int N);

void Build_K_elemental(double Ke[], int n, double A, double E, double l, double I);

void Build_K_global_banded(double Kb[], double A, double E, double l, double L, double I, int Nx, int N, int kl, int ku, int lda, int bfr);

void Build_F_elemental(double Fe[], double qx, double qy, double l, double time);

void Build_F_global(double F[], double Fy, double qx, double qy, double time, double l, int Nx, int N);

void Build_M_elemental(double Me[], int Ne, double rho, double A, double l);

void Build_M_global_banded(double Mb[], double rho, double A, double l, int Nx, int N);

void Matrix_System_Solver(double A[], double b[], int N);

void Banded_Matrix_Solver(double Ab[], double b[], int N,int lda, int kl, int ku);

void Write_Vector(double F[], int N, double l, double L, string mystring);

void Write_Point_Displacement(double F[], int N, string mystring);

void Build_Fn(double F[], double Fn[], double del_t, int N);

void Build_Un1_Multiplier(double Un1[], double Kb[], double Mb[], double u1[], int N, int lda, double del_t, int kl, int ku);

void Build_Un0_Multiplier(double Un0[], double Mb[], double u0[], int N);

void Build_Multiplier1(double S[], double F[], double Kb[], double M[], double u0[], double u1[], double del_t, int N, int lda, int kl, int ku);

double RMS_error(double u1[], double S[], double M[], int N); 

void Build_K_eff(double Kb[], double Mb[], double del_t, double beta, int N, int ku, int kl, int lda);

void Build_Multiplier2(double S[], double F[], double Mb[], double u0[], double udot[], double uddot0[], double del_t, double beta, int N);

void Build_uddot(double uddot1[], double uddot0[], double u1[], double u0[], double udot[], double beta, double del_t, int N);

void Build_udot(double udot[], double u0[], double uddot0[] , double uddot1[], double del_t, double gamma, int N);