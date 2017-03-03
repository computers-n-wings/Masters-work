#include <iostream>
#include <cmath>
#include <string>
#include "myfunctions.h"
using namespace std;

int main()
{
	// ########################################## Parameters ################################################
	// ######################################################################################################
	int time_dependence = 1;
	int scheme = 0;

	const double L = 10000.0; 
	const int Nx = 24;
	const double A = 12000.0;
	const double I = 14400000.0; 
	const double E = 2.1e+06;
	const double qx = 0.0;
	const double qy = 1.0;
	const double l = L/(double)Nx;
	const double Fy = 1000.0;
	const double T = 1.0;
	const int Nt = 1000;
	const double rho = 7850.0e-09;
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
   		Build_F_global(F, Fe, Nx, Fy, 1.0, N);

   		Matrix_System_Solver(K, F, N);
   		Write_Text_File(F,N, l, L, "Task1");
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
   			double * u0 = new double[N];
   			double * u1 = new double[N];
   			double * S = new double[N];
   			
   			Build_K_elemental(Ke, Ne, A, E, l, I);
   			Build_global_matrix(K, Ke, Nx, N);
   			Build_M_elemental(Me, Ne, rho, A, l);
   			Build_global_matrix(M, Me, Nx, N);

   			for (int i = 0; i<Nt; ++i)
   			{
   				Build_F_elemental(Fe, qx, qy, l, (double)(i+1)*del_t/T);
   				Build_F_global(F, Fe, Nx, Fy, (double)(i+1)*del_t/T, N);
   				Build_Multiplier(S, F, K, M, u0, u1, del_t, N);
   				Copy_Vector(u1, u0, N);
   				Matrix_System_Solver(M, S, N);
   				Copy_Vector(S, u1, N);
   				if (i % (Nt/10) == 0)
   				{
   					Write_Text_File(u1,N, l, L, "Task2_Time_Step" + to_string(i));
   				}					   				
   			}
   			Print_Vector(u1, N);

   			
   			delete[] Fe;
   			delete[] F;
   			delete[] Ke;
   			delete[] Me;
   			delete[] K;
   			delete[] M;
   			delete[] u0;
   			delete[] u1;
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