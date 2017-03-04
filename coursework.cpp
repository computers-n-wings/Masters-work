#include <iostream>
#include <cmath>
#include <string>
#include "myfunctions.h"
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
// int main()
int main(int argc, char* argv[])
{
	po::options_description desc("Solving for displacement along a beam");
    desc.add_options()
        ("L",  po::value<double>(),    "Length of the beam")
        ("Nx", po::value<int>(), "Number of elements")
        ("A", po::value<double>(), "Area of cross-section")
        ("I", po::value<double>(), "Second moment of inertia")
        ("E", po::value<double>(), "Young's modulus")
        ("time_dependence", po::value<int>(), "Time dependence: 0 for static, 1 for dynamic")
        ("scheme", po::value<int>(), "Integration scheme: 0 for explicit, 1 for implicit")
        ("T", po::value<double>(), "Number of seconds if dynamic problem")
        ("Nt", po::value<int>(), "Number of time steps if dynamic problem")
        ("help", 					"Produce help message");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) 
        {
        	cout << "Solving for displacement along a beam." << endl;
        	cout << desc << endl;
        	return 0;
    	}

    	const double L = vm.count("L") ? vm["L"].as<double>() : 10000.0;
    	const int Nx = vm.count("Nx") ? vm["Nx"].as<int>() : 24;
    	const double A = vm.count("A") ? vm["A"].as<double>() : 12000.0;
    	const double I = vm.count("I") ? vm["I"].as<double>() : 14400000.0;
    	const double E = vm.count("E") ? vm["E"].as<double>() : 2.1e+06;
    	const int time_dependence = vm.count("time_dependence") ? vm["time_dependence"].as
    	<int>() : 0;
    	const int scheme = vm.count("scheme") ? vm["scheme"].as<int>() : 0;
    	const double T = vm.count("T") ? vm["T"].as<double>() : 1.0;
    	const int Nt = vm.count("Nt") ? vm["Nt"].as<int>() : 10000;

    // ########################################## Parameters ################################################
	// ######################################################################################################
	const double qx = 0.0;
	const double qy = -1.0;
	const double l = L/(double)Nx;
	const double Fy = -1000.0;
	const double rho = 7850.0e-09;
	const int Ne = 6; 	
	const int N = (Nx-1)*(Ne/2);											
	
	// ################################### Find nodal displacements #########################################
	// ######################################################################################################
    
   if (time_dependence == 0)
   {
   		double * Ke = new double[Ne*Ne](); 				
		double * Fe = new double[Ne](); 
		double * K = new double[N*N]();					
		double * F = new double[N]();

		int ku = 4; int kl = 4;
		int m = ku+2*kl+1;
		double * Kb = new double[m*N]();	
   		Build_K_elemental(Ke, Ne, A, E, l, I);
   		Build_global_matrix(K, Ke, Nx, N);	
   		Print_Matrix(K, N, N);
		cout << endl;
		Banded_Storage(Kb, K, N, kl, ku);
		Print_Matrix(Kb, m, N);
		cout << endl;
   		// Build_F_elemental(Fe, qx, qy, l, 1.0);
   		// Build_F_global(F, Fe, Nx, Fy, 1.0, N);

   		// Matrix_System_Solver(K, F, N);
   		// Write_Vector(F,N, l, L, "Task1");

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

   			double * Ke = new double[Ne*Ne](); 				
			double * Fe = new double[Ne](); 
			double * Me = new double[Ne*Ne]();
			double * K = new double[N*N]();
			double * F = new double[N]();						
   			double * M = new double[N*N]();
   			double * u0 = new double[N]();
   			double * u1 = new double[N]();
   			double * S = new double[N]();
   			double * d = new double[Nt]();
   			
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
   				// if (i % (Nt/10) == 0)
   				// {
   				// 	Write_Vector(u1,N, l, L, "Task2_Time_Step" + to_string(i));
   				// }
   				d[i] = u1[(N-1)/2];					   				
   			}
   			Write_Point_Displacement(d, Nt, "Task2");
   			// Print_Vector(u1, N);

   			
   			delete[] Fe;
   			delete[] F;
   			delete[] Ke;
   			delete[] Me;
   			delete[] K;
   			delete[] M;
   			delete[] u0;
   			delete[] u1;
   			delete[] S;
   			delete[] d;

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