#include <iostream>
#include <string>
#include "myfunctions.h"
using namespace std;

#include <mpi.h>

#define F77NAME(x) x##_
extern "C" 
{
  void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, 
                const int& nrhs, const double * A, const int& ldab, 
                int * ipiv, double * B, const int& ldb, int& info);
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
  double F77NAME(dnrm2) (const int& n, const double* a, const int& incx);
  void F77NAME(dgbmv)(const char& trans, const int& m, const int& n,
                const int& kl, const int& ku,
                const double& alpha, const double* a, const int& lda,
                const double* x, const int& incx, const double& beta,
                double* y, const int& incy);
} 

void Static_Solver(double A, double E, double I, double L, double l, double qx, double qy, double Fy, int Nx, int ku, int kl, int N)
{
	const int lda = 1 + 2*kl + ku;

	double * Kb = new double[lda*N]();
	double * F = new double[N]();

	Build_K_global_banded(Kb, A, E, l, L, I, Nx, N, kl, ku, lda, kl);
  Build_F_global(F, Fy, qx, qy, 1.0, l, Nx, N);
  Banded_Matrix_Solver(Kb, F, N, lda, kl, ku);
  Write_Vector(F,N, l, L, "Task1");

  delete[] Kb;
  delete[] F;
}

void Dynamic_Solver_1(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N)
{
 	const int lda = 1 + kl + ku;
  const double del_t = T/(double)Nt;
 		
  double * Mb = new double[N]();
  double * u0 = new double[N]();
  double * u1 = new double[N]();
  double * F = new double[N]();
 	double * Kb = new double[lda*N]();
 	double * d = new double[Nt]();

 	Build_M_global_banded(Mb, rho, A, l, Nx, N);
 	Build_K_global_banded_altered(Kb, Mb, A, E, l, L, I, del_t, Nx, N, kl, ku, lda, 0);
 			
 	for (int i = 0; i<Nt; ++i)
 	{
 		Build_F_global(F, Fy, qx, qy, (double)i*del_t/T, l, Nx, N);
 		Build_Multiplier1(F, Kb, Mb, u0, u1, del_t, N, lda, kl, ku);
    F77NAME(dcopy)(N, u1, 1, u0, 1);
 		for (int j = 0; j<N; ++j)
 		{
 			u1[j] = F[j]/Mb[j];
 		}
 		d[i] = u1[(N-1)/2];					   				
 	}
 	Write_Point_Displacement(d, Nt, "Task2");

 	delete[] Mb;
 	delete[] u1;
 	delete[] F;
  delete[] Kb;
  delete[] d;
}

void Dynamic_Solver_2(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N)
{
  const int lda = 1 + 2*kl + ku;
  const double del_t = T/(double)Nt;
 	const double beta = 0.25;
 	const double gamma = 0.5;
 	const double coeff1 = 1.0/(beta*del_t*del_t);
 	const double coeff2 = 1.0/(beta*del_t);
 	const double coeff3 = (1.0-2.0*beta)/(2.0*beta);
 	const double coeff4 = (1.0-gamma)*del_t;
 	const double coeff5 = gamma*del_t;

 	double * Mb = new double[N]();
 	double * Keff = new double[lda*N]();
 	double * F = new double[N]();
 	double * u0 = new double[N]();
 	double * u1 = new double[N]();
 	double * udot = new double[N]();
 	double * udotdot0 = new double[N]();
 	double * udotdot1 = new double[N]();

  Build_M_global_banded(Mb, rho, A, l, Nx, N);
 	Build_Keff(Keff, Mb, coeff1, A, E, l, L, I, Nx, N, kl, ku, lda, kl);

  for (int i = 0; i<Nt; ++i)
  {
    Build_F_global(F, Fy, qx, qy, (double)(i+1)*del_t/T, l, Nx, N);
    Build_Multiplier2(u1, F, Mb, u0, udot, udotdot0, coeff1, coeff2, coeff3, N);
    Banded_Matrix_Solver(Keff, u1, N, lda, kl, ku);
    Build_udotdot(udotdot1, u1, u0, udot, udotdot0, coeff1, coeff2, coeff3, N);
    Build_udot(udot, udotdot0, udotdot1, coeff4, coeff5, N);
    F77NAME(dcopy) (N, u1, 1, u0, 1);
    F77NAME(dcopy) (N, udotdot1, 1, udotdot0, 1);
  }
  Write_Vector(u0, N, l, L, "Task1");

 	delete[] Keff;
 	delete[] u0;
 	delete[] u1;
 	delete[] udot;
 	delete[] udotdot0;
 	delete[] udotdot1;
}

void Dynamic_Solver_MPI_1(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N, int size, int rank)
{
	const int lda = 1 + kl + ku;
	const double del_t = T/(double)Nt;
	const int Nxloc = Nx/2 + 2;
	const int Nloc = (Nxloc-1)*3 ;

	double * Mbloc = new double[Nloc]();
	double * u0loc = new double[Nloc]();
  double * u1loc = new double[Nloc]();
	double * Kbloc = new double[Nloc*lda]();
	double * Floc = new double[Nloc]();

  Build_M_global_banded_MPI(Mbloc, rho, A, l, Nx, Nxloc, Nloc, size, rank);
  Build_K_global_banded_altered_MPI(Kbloc, Mbloc, A, E, l, L, I, del_t, Nx, Nxloc, N, Nloc, kl, ku, lda, 0, size, rank);

  for (int i = 0; i<Nt; ++i)
  {
    Build_F_global_MPI(Floc, Fy, qx, qy, (double)i*del_t/T, l, Nx, Nxloc, N, Nloc, size, rank);
    Build_Multiplier1_MPI(Floc, Kbloc, Mbloc, u0loc, u1loc, del_t, N, Nloc, lda, kl, ku, size, rank);
    F77NAME(dcopy)(Nloc, u1loc, 1, u0loc, 1);
    for (int j = 0; j<Nloc; ++j)
    {
      u1loc[j] = Floc[j]/Mbloc[j];
    }  
    MPI_Barrier(MPI_COMM_WORLD); 
    if (rank == 0)
    {
      MPI_Send(&u1loc[(Nx/2 - 2)*3], 3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&u1loc[3*Nx/2], 3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (rank == 1)
    {
      MPI_Recv(&u1loc[0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&u1loc[6], 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }  
    MPI_Barrier(MPI_COMM_WORLD);     
  }
  Print_Vector_Parallel(u1loc, Nloc, size, rank);

  if (rank == 0)
  {
    double * Mb = new double[N]();
    double * u0 = new double[N]();
    double * u1 = new double[N]();
    double * F = new double[N]();
    double * Kb = new double[lda*N]();

    Build_M_global_banded(Mb, rho, A, l, Nx, N);
    Build_K_global_banded_altered(Kb, Mb, A, E, l, L, I, del_t, Nx, N, kl, ku, lda, 0);

    for (int i = 0; i<Nt; ++i)
    {
      Build_F_global(F, Fy, qx, qy, (double)i*del_t/T, l, Nx, N);
      Build_Multiplier1(F, Kb, Mb, u0, u1, del_t, N, lda, kl, ku);
      F77NAME(dcopy)(N, u1, 1, u0, 1);
      for (int j = 0; j<N; ++j)
      {
        u1[j] = F[j]/Mb[j];
      }
    }
    cout << endl;
    cout << "True Solution" << endl;
    Print_Vector(u1, N);                
 	}
  delete[] Mbloc;
  delete[] u0loc;
  delete[] u1loc;
  delete[] Kbloc;
  delete[] Floc;
}

void Dynamic_Solver_MPI_2(double A, double E, double I, double L, double l, double qx, double qy, double Fy, double rho, double T, int Nt, int Nx, int ku, int kl, int N, int size, int rank)
{
  const int nb = (int)(N/size) + 1;
  const int lda = 1 + 2*kl + 2*ku;
  const double del_t = T/(double)Nt;
  const double beta = 0.25;
  const double gamma = 0.5;
  const double coeff1 = 1.0/(beta*del_t*del_t);
  const double coeff2 = 1.0/(beta*del_t);
  const double coeff3 = (1.0-2.0*beta)/(2.0*beta);
  const double coeff4 = (1.0-gamma)*del_t;
  const double coeff5 = gamma*del_t;

  double * Mbloc = new double[nb]();
  double * Kbloc = new double [nb*lda]();
  double * Floc = new double[nb]();
  double * u0loc = new double[nb]();
  double * u1loc = new double[nb]();
  double * udotloc = new double[nb]();
  double * udotdot0loc = new double[nb]();
  double * udotdot1loc = new double[nb]();
  double * Mb = new double[N]();
  double * Keff = new double[lda*N]();
  double * F = new double[N]();

  Build_M_global_banded(Mb, rho, A, l, Nx, N);
  Build_Keff(Keff, Mb, coeff1, A, E, l, L, I, Nx, N, kl, ku, lda, kl+ku);
  Build_Block_Array(Mbloc, Mb, N, N, nb, nb, 1, 0, 0, 1.0, size, rank);
  Build_Block_Array(Kbloc, Keff, N, lda*N, nb, nb*lda, lda, kl+ku, ku, 1.0, size, rank);

  for (int i = 0; i < Nt; ++i)
  {
    Build_F_global(F, Fy, qx, qy, (double)(i+1)*del_t/T, l, Nx, N);
    Build_Block_Array(Floc, F, N, N, nb, nb, 1, 0, 0, 0.0, size, rank);
    Build_Multiplier2(u1loc, Floc, Mbloc, u0loc, udotloc, udotdot0loc, coeff1, coeff2, coeff3, nb);
    Banded_Matrix_Solver_Parallel(Kbloc, u1loc, size*nb, nb, lda, kl, ku, size);
    Build_udotdot(udotdot1loc, u1loc, u0loc, udotloc, udotdot0loc, coeff1, coeff2, coeff3, nb);
    Build_udot(udotloc, udotdot0loc, udotdot1loc, coeff4, coeff5, nb);
    F77NAME(dcopy) (nb, u1loc, 1, u0loc, 1);
    F77NAME(dcopy) (nb, udotdot1loc, 1, udotdot0loc, 1);
  }
  
  Print_Vector_Parallel(u0loc, nb, size, rank);

}