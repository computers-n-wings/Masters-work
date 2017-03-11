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
	const int N = (Nx-1)*3;
	const int ku = 4; 
	const int kl = 4;
							
	double * F = new double[N]();
	
	// ################################### Find nodal displacements #########################################
	// ######################################################################################################
    
   if (time_dependence == 0)
   {
		const int lda = 1 + 2*kl + ku;

		double * Kb = new double[lda*N]();											

		Build_K_global_banded(Kb, A, E, l, L, I, Nx, N, kl, ku, lda, kl);
		Print_Matrix(Kb, N, lda);
   		Build_F_global(F, Fy, qx, qy, 1.0, l, Nx, N);
   		Banded_Matrix_Solver(Kb, F, N, lda, kl, ku);
   		// Write_Vector(F,N, l, L, "Task1");

   }
   else if (time_dependence == 1)
   {
   		const double del_t = T/(double)Nt;
   		const double rho = 7850.0e-09;
   		
     	double * Mb = new double[N]();
     	double * u0 = new double[N]();
        double * u1 = new double[N]();

   		Build_M_global_banded(Mb, rho, A, l, Nx, N);

   		if (scheme == 0)
   		{	
   			const int lda = 1 + kl + ku;

   			double * Kb = new double[lda*N]();
        	double * S = new double[N]();
   			double * d = new double[Nt]();

   			Build_K_global_banded(Kb, A, E, l, L, I, Nx, N, kl, ku, lda, 0);
   			
   			for (int i = 0; i<Nt; ++i)
   			{
   				Build_F_global(F, Fy, qx, qy, (double)i*del_t/T, l, Nx, N);
   				Build_Multiplier1(S, F, Kb, Mb, u0, u1, del_t, N, lda, kl, ku);
   				Copy_Vector(u1, u0, N);
   				Banded_Matrix_Solver(Mb, S, N, 1, 0, 0);
   				Copy_Vector(S, u1, N);
   				d[i] = u1[(N-1)/2];					   				
   			}
   			Write_Point_Displacement(d, Nt, "Task2");

   		}
   		else if (scheme == 1)
   		{
   			const int lda = 1 + 2*kl + ku;
   			const double beta = 0.25;
   			const double gamma = 0.5;

   			double * udot = new double[N]();
   			double * uddot0 = new double[N]();
   			double * uddot1 = new double[N]();
   			double * S = new double[N]();
   			double * Kb = new double[lda*N]();

   			Build_K_global_banded(Kb, A, E, l, L, I, Nx, N, kl, ku, lda, kl);
   			Build_K_eff(Kb, Mb, del_t, beta, N, ku, kl, lda);
   			
   			for (int i = 0; i<Nt; ++i)
   			{
   				Build_F_global(F, Fy, qx, qy, (double)(i+1)*del_t/T, l, Nx, N);
   				Build_Multiplier2(S, F, Mb, u0, udot, uddot0, del_t, beta, N);
   				Banded_Matrix_Solver(Kb, S, N, lda, kl, ku);
   				Copy_Vector(S, u1, N);
   				Build_uddot(uddot1, uddot0, u1, u0, udot, beta, del_t, N);
				Build_udot(udot, u0, uddot0, uddot1, del_t, gamma, N);
				Copy_Vector(u1, u0, N);
				Copy_Vector(uddot1, uddot0, N);	
   			}
   			Write_Vector(u1,N, l, L, "Task1");


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