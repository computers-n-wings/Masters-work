#include <iomanip>
#include <iostream>
#include <fstream>
using namespace std;

#define F77NAME(x) x##_
extern "C" 
{
	void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B, const int& ldb,
                    int& info);
	void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, 
                    const int& nrhs, const double * A, const int& ldab, 
                    int * ipiv, double * B, const int& ldb, int* info);
	void F77NAME(dcopy)(const int& n, const double * x, const int& incx, 
						const double * y, const int& incy);
	void F77NAME(dscal) (const int& n, const double& alpha, const double * x,
						const int& incx);
	void F77NAME(daxpy) (const int& n, const double& alpha, const double * x,
						const int& incx, const double * y, const int& incy);
	void F77NAME(dgemv) (const char& trans,  const int& m,
         				const int& n,       const double& alpha,
         				const double* a,    const int& lda,
         				const double* x,    const int& incx,
         				const double& beta, double* y, const int& incy);
}	


void Print_Matrix(double A[], int m, int n) 
{
	for (int i = 0; i<m; ++i)
	{
		for (int j = 0; j<n; ++j)
		{
			cout << setprecision(3) << setw(15) << A[i*n+j];
		} 
		cout << endl;
	}
}

void Print_Vector(double b[], int n)
{
	for (int i = 0; i<n; ++i)
	{
		cout << setw(15) << b[i] << endl;
	}
}

void Build_K_elemental(double Ke[], int n, double A, double E, double l, double I)
{
	double K1 = A*E/l; 							
	double K2 = 12.0*E*I/(l*l*l); 				
	double K3 = 6.0*E*I/(l*l); 					
	double K4 = E*I/l;
	for (int i = 0; i<n; ++i)
	{
		for (int j = 0; j<n; ++j)
		{
			if (i <= j)
			{
				Ke[0*n+0] = K1;
				Ke[0*n+3] = -K1;
				Ke[1*n+1] = K2;
				Ke[1*n+2] = K3;
				Ke[1*n+4] = -K2;
				Ke[1*n+5] = K3;
				Ke[2*n+2] = 4.0*K4;
				Ke[2*n+4] = -K3;
				Ke[2*n+5] = 2.0*K4;
				Ke[3*n+3] = K1;
				Ke[4*n+4] = K2;
				Ke[4*n+5] = -K3;
				Ke[5*n+5] = 4.0*K4;
			}
			else
			{
				Ke[i*n+j] = Ke[j*n+i];
			}
		}
	}
}

void Build_F_elemental(double Fe[], double qx, double qy, double l, double time)
{
	const double F1 = time*qx*l/2.0; 							
	const double F2 = time*qy*l/2.0; 							
	const double F3 = time*qy*l*l/12.0; 	
	Fe[0] = F1;
	Fe[1] = F2;
	Fe[2] = F3;
	Fe[3] = F1;
	Fe[4] = F2;
	Fe[5] = -F3;
}

void Build_global_matrix(double K[], double Ke[], int Nx, int N)
{
	for (int k = 0; k<Nx; ++k)
	{
		if (k == 0)
		{
			for (int i = 0; i<3; ++i)
			{
				for (int j = 0; j<3; ++j)
				{
					K[i*N + j] += Ke[(i+3)*6+j+3];
				}
			}
		}
		else if (k == Nx-1)
		{
			for (int i = 0; i<3; ++i)
			{
				for (int j = 0; j<3; ++j)
				{
					K[(3*(k-1) + i)*N + 3*(k-1) + j] += Ke[i*6+j];
				}
			}
		}
		else
		{
			for (int i = 0; i<6; ++i)
			{
				for (int j = 0; j<6; ++j)
				{
					K[(3*(k-1) + i)*N + 3*(k-1) + j] += Ke[i*6+j];
				}
			}
		}
	}
}

void Build_K_global_banded(double Kb[], double Ke[], int Nx, int N, int kl, int ku)
{
	for (int k = 0; k<Nx; ++k)
	{
		if (k == 0)
		{
			for (int i = 0; i<3; ++i)
			{
				for (int j = 0; j<3; ++j)
				{
					Kb[(ku+i-j)*N + 3*(k-1) + j] += Ke[(i+3)*6+j+3];
				}
			}
		}
		else if (k == Nx-1)
		{
			for (int i = 0; i<3; ++i)
			{
				for (int j = 0; j<3; ++j)
				{
					Kb[(3*(k-1) + (ku+i-(3*(k-1) + j)))*N + 3*(k-1) + j] += Ke[i*6+j];
				}
			}
		}
		else
		{
			for (int i = 0; i<6; ++i)
			{
				for (int j = 0; j<6; ++j)
				{
					Kb[(3*(k-1) + (ku+i-(3*(k-1) + j)))*N + 3*(k-1) + j] += Ke[i*6+j];
				}
			}
		}
	}
}

void Build_F_global(double F[], double Fe[], int Nx, double Fy, double time)
{
	for (int k = 0; k<Nx; ++k)
	{
		if (k == 0)
		{
			F[0] = Fe[3];
			F[1] = Fe[4];
			F[2] = Fe[5];
		}
		else if (k == Nx-1)
		{
			for (int i = 0; i<3; ++i)
			{
				F[3*(k-1)+i] += Fe[i];
			}
		}
		else
		{
			for (int i = 0; i<6; ++i)
			{
				F[3*(k-1) +i] += Fe[i];
			}
		}
		if (k == Nx/2)
		{
			F[3*(k-1) +1] += time*Fy;
		}
	}
}

void Banded_Storage(double Kb[], double K[], int N, int kl, int ku)
{
	for (int i = 0; i<N; ++i)
	{
		for (int j = 0; j<N; ++j)
		{
			Kb[(ku+i-j)*N+j] = K[i*N+j];
		}
	}
}

void Matrix_System_Solver(int N, double K[], double F[])
{
	const int nrhs = 1;
    int info = 0;
    int * ipiv = new int[N];
    F77NAME(dgesv) (N, nrhs, K, N, ipiv, F, N, info);
}

void Write_Text_File(double F[], int N, double l, double L)
{
	ofstream File;
	File.open("displacements.txt");

	if (File.good())
	{
		File << 0 << ' ' << setw(5) << setprecision(5) << 0 << '\n';
		for (int i = 0; i < N/3; ++i)
		{
			File << l*(i+1) << ' ' << setw(5) << setprecision(5) << F[3*i+1] << '\n';
		}
		File << L << ' ' << setw(5) << setprecision(5) << 0;
	}
	else
	{
		cout << "The file is not good." << endl;
	}

	File.close();
}

void Build_M_elemental(double Me[], int Ne, double rho, double A, double l)
{
	double alpha = 0.0416667;
	double M1 = 0.5*rho*A*l;
   	double M2 = rho*A*alpha*l*l*l;

	Me[0] = M1;
	Me[Ne + 1] = M1;
	Me[2*Ne + 2] = M2;
	Me[3*Ne + 3] = M1;
	Me[4*Ne + 4] = M1;
	Me[5*Ne + 5] = M2;
}

void Build_Un1_Multiplier(double K_temp[], double K[], double M[], int N, double del_t)
{
	F77NAME(dcopy) (N*N, K, 1, K_temp, 1);
	F77NAME(dscal) (N*N, -del_t*del_t, K_temp, 1);
	F77NAME(daxpy) (N*N, 2.0, M, 1, K_temp, 1);
}

void Build_Multiplier(double S[], double F[], double K_temp[], double M[], double del_t, double u0[], double u1[], int N)
{
	double * F_temp = new double[N];
	double * U_n1 = new double[N];
	double * U_n0 = new double[N];

	F77NAME(dcopy) (N, F, 1, F_temp, 1);
	F77NAME(dscal) (N, del_t*del_t, F_temp, 1);
	F77NAME(dgemv) ('N', N, N, 1.0, K_temp, N, u1, 1, 0.0, U_n1, 1);
	F77NAME(dgemv) ('N', N, N, -1.0, M, N, u0, 1, 0.0, U_n0, 1);
	F77NAME(daxpy) (N, 1.0, U_n0, 1, U_n1, 1);
	F77NAME(daxpy) (N, 1.0, U_n1, 1, F_temp, 1);
	F77NAME(dcopy) (N, F_temp, 1, S, 1);

	delete[] F_temp;
	delete[] U_n1;
	delete[] U_n0;
}

void Copy_Vector(int N, double a[], double b[])
{
	F77NAME(dcopy) (N, a, 1, b, 1);
}