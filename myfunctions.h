#include <string>
using namespace std;

void Print_Matrix(double A[], int n, int m);

void Print_Vector(double b[], int n);

void Zero_Vector(double a[], int N);

void Copy_Vector(double a[], double b[], int N);

void Build_K_elemental(double Ke[], int n, double A, double E, double l, double I);

void Build_F_elemental(double Fe[], double qx, double qy, double l, double time);

void Build_global_matrix(double K[], double Ke[], int Nx, int N);

void Build_K_global_banded(double Kb[], double Ke[], int Nx, int N, int kl, int ku);

void Banded_Storage(double Kb[], double K[], int N, int kl, int ku);

void Build_F_global(double F[], double Fe[], int Nx, double Fy, double time, int N);

void Matrix_System_Solver(double A[], double b[], int N);

void Write_Vector(double F[], int N, double l, double L, string mystring);

void Write_Point_Displacement(double F[], int N, string mystring);

void Build_M_elemental(double Me[], int Ne, double rho, double A, double l);

void Build_Fn(double F[], double Fn[], double del_t, int N);

void Build_Un1_Multiplier(double Un1[], double K[], double M[], double u1[], int N, double del_t);

void Build_Un0_Multiplier(double Un0[], double M[], double u0[], int N);

void Build_Multiplier(double S[], double F[], double K[], double M[], double u0[], double u1[], double del_t, int N);

double RMS_error(double u1[], double S[], double M[], int N); 

