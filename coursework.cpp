#include <iostream>
#include <cmath>
#include "myfunctions.h"
using namespace std;

int main()
{
	// ########################################## Parameters ################################################
	// ######################################################################################################

	const double L = 10000.0; 
	const int Nx = 24;
	const double A = 12000.0;
	const double I = 14400000.0; 
	const double E = 2.1e+06;
	const double qx = 0.0;
	const double qy = 1.0;
	const double l = L/(double)Nx;
	const double Fy = 1000.0;
	int time_dependence = 0;
	int scheme = 0;

	const double T = 1.0;
	const int Nt = 1000;
	const double rho = 7850.0e-12;
	const int Ne = 6; 	
	const int N = (Nx-1)*(Ne/2);											
	
	// ################################### Find nodal displacements #########################################
	// ######################################################################################################
    
   if (time_dependence == 0)
   {
   		double * Ke = new double[Ne*Ne]; 				
		double * Fe = new double[Ne]; 
		double * K = new double[N*N];					
		double * F = new double[N];	

   		Build_K_elemental(Ke, Ne, A, E, l, I);
   		Build_global_matrix(K, Ke, Nx, N);	
   		Build_F_elemental(Fe, qx, qy, l, 1.0);
   		Build_F_global(F, Fe, Nx, Fy, 1.0);

   		Matrix_System_Solver(N, K, F);
   		Write_Text_File(F,N, l, L);
   		Print_Vector(F, N);

   		delete[] Ke;
   		delete[] Fe;
   		delete[] K;
   		delete[] F;
   }
   else if (time_dependence == 1)
   {
   		if (scheme == 0)
   		{
   			const double del_t = T/(double)Nt;

   			double * Ke = new double[Ne*Ne]; 				
			double * Fe = new double[Ne]; 
			double * Me = new double[Ne*Ne];
			double * K = new double[N*N];
			double * F = new double[N];						
   			double * M = new double[N*N];
   			double * S = new double[N];
   			double * u0 = new double[N];
   			double * u1 = new double[N];
   			double * K_temp = new double[N*N];

   			Build_K_elemental(Ke, Ne, A, E, l, I);
   			Build_global_matrix(K, Ke, Nx, N);
   			Build_M_elemental(Me, Ne, rho, A, l);
   			Build_global_matrix(M, Me, Nx, N);

   			Build_Un1_Multiplier(K_temp, K, M, N, del_t);

   			for (int i = 0; i<Nt; ++i)
   			{
   				double * F = new double[N];	
   				double * S = new double[N];
   				Build_F_elemental(Fe, qx, qy, l, del_t*(i+1));
   				Build_F_global(F, Fe, Nx, Fy, del_t*(i+1));
   				Build_Multiplier(S, F, K_temp, M, del_t, u0, u1, N);
   				Copy_Vector(N, u1, u0);
   				Matrix_System_Solver(N, M, S);
   				Copy_Vector(N, S, u1);
   				Print_Vector(u0, N);
   				cout << endl;
   				Print_Vector(u1, N);
   				cout << endl;
   			}

   			delete[] Fe;
   			delete[] F;
   			delete[] Ke;
   			delete[] Me;
   			delete[] K;
   			delete[] M;
   			delete[] u0;
   			delete[] u1;
   			delete[] K_temp;
   			delete[] S;
   		}
   		else if (scheme == 1)
   		{
   			cout << "The scheme has not been defined" << endl;
   		}
   		else
   		{
   			cout << "The integration scheme has not been correctly specified" << endl;
   		}
   }
   else 
   {
   		cout << "The time dependence has not been correctly specified." << endl;
   }


	

	return 0;

}