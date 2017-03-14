#include <iostream>
#include <cmath>
#include <string>
#include "myfunctions.h"
#include "solvers.h"
using namespace std;

#include <mpi.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
        ("parallel", po::value<int>(), "Parallelization: 0 for serial, 1 for parallel")
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
      const int parallel = vm.count("parallel") ? vm["parallel"].as<int>() : 0;

    // ########################################## Parameters ################################################
	// ######################################################################################################
	const double qx = 0.0;
	const double qy = -1000.0;
	const double l = L/(double)Nx;
	const double Fy = -1000.0; 	
	const int N = (Nx-1)*3;
	const int ku = 4; 
	const int kl = 4;
	
	// ################################### Find nodal displacements #########################################
	// ######################################################################################################
    
   if (time_dependence == 0)
   {
      Static_Solver(A, E, I, L, l, qx, qy, Fy, Nx, ku, kl, N);
   }
   else if (time_dependence == 1)
   {
      if (scheme == 0)
      {
        if (parallel == 0)
        {
          Dynamic_Solver_1(A, E, I, L, l, qx, qy, Fy, rho, T, Nt, Nx, ku, kl, N);
        }
        else if (parallel == 1)
        {
          // Problem with sending the number of processors to solver
          // Cannot convert from char to int
          int np = (int)argv[0];
          Dynamic_Solver_MPI_1(A, E, I, L, l, qx, qy, Fy, rho, T, Nt, Nx, ku, kl, N, np);
        }
        else
        {
          cout << "The parallelization has not been correctly specified" << endl;
        }          
      }
   		else if (scheme == 1)
   		{
        if (parallel == 0)
        {
          Dynamic_Solver_2(A, E, I, L, l, qx, qy, Fy, rho, T, Nt, Nx, ku, kl, N);
        }
        else if (parallel == 1)
        {
          cout << endl;
        }
        else
        {
          cout << "The parallelization has not been correctly specified" << endl;
        }                
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