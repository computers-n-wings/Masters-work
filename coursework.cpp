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
        ("rho", po::value<double>(), "Density")
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

    	const double L = vm.count("L") ? vm["L"].as<double>() : 10.0;
    	const int Nx = vm.count("Nx") ? vm["Nx"].as<int>() : 24;
    	const double A = vm.count("A") ? vm["A"].as<double>() : 0.012;
    	const double I = vm.count("I") ? vm["I"].as<double>() : 1.44e-05;
    	const double E = vm.count("E") ? vm["E"].as<double>() : 2.1e+11;
      const double rho = vm.count("rho") ? vm["rho"].as<double>() : 7850.0;
    	const int time_dependence = vm.count("time_dependence") ? vm["time_dependence"].as
    	<int>() : 0;
    	const int scheme = vm.count("scheme") ? vm["scheme"].as<int>() : 0;
    	const double T = vm.count("T") ? vm["T"].as<double>() : 1.0;
    	const int Nt = vm.count("Nt") ? vm["Nt"].as<int>() : 10000;

    // ########################################## Parameters ################################################
	// ######################################################################################################
	const double qx = 0.0;
	const double qy = -1000.0;
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
   	Build_F_global(F, Fy, qx, qy, 1.0, l, Nx, N);
   	Banded_Matrix_Solver(Kb, F, N, lda, kl, ku);
   	Write_Vector(F,N, l, L, "Task1");

   }
   else if (time_dependence == 1)
   {
   		const double del_t = T/(double)Nt;
   		
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
   			Write_Vector(u1,N, l, L, "Task1");
   			// Write_Point_Displacement(d, Nt, "Task2");

   		}
   		else if (scheme == 1)
   		{
   			const int lda = 1 + 2*kl + ku;
   			const double beta = 0.25;
   			const double gamma = 0.5;
   			const double coeff1 = 1/(beta*del_t*del_t);
   			const double coeff2 = 1/(beta*del_t);
   			const double coeff3 = (1-2*beta)/(2*beta);
   			const double coeff4 = (1-gamma)*del_t;
   			const double coeff5 = gamma*del_t;

   			double * Keff = new double[lda*N]();
   			double * u0 = new double[N]();
   			double * u1 = new double[N]();
   			double * udot0 = new double[N]();
        double * udot1 = new double[N]();
   			double * udotdot0 = new double[N]();
   			double * udotdot1 = new double[N]();
        double * S = new double[N]();

   			Build_Keff(Keff, Mb, coeff1, A, E, l, L, I, Nx, N, kl, ku, lda, kl);

        for (int i = 0; i<Nt; ++i)
        {
            Build_F_global(F, Fy, qx, qy, (double)(i+1)*del_t/T, l, Nx, N);
            Build_Multiplier2(S, Mb, F, u0, udot0, udotdot0, coeff1, coeff2, coeff3, N);
            Banded_Matrix_Solver(Keff, S, N, lda, kl, ku);
            Copy_Vector(S, u1, N);
            Build_udotdot(udotdot1, u1, u0, udot0, udotdot0, coeff1, coeff2, coeff3, N);
            Build_udot(udot1, udot0, udotdot0, udotdot1, coeff4, coeff5, N);
            Copy_Vector(u1, u0, N);
            Copy_Vector(udot1, udot0, N);
            Copy_Vector(udotdot1, udotdot0, N);
        }
        Write_Vector(u0, N, l, L, "Task1");

   			delete[] Keff;
   			delete[] u0;
   			delete[] u1;
   			delete[] udot0;
        delete[] udot1;
   			delete[] udotdot0;
   			delete[] udotdot1;
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