#include "DiscreteProblem.hpp"
#include "Intel_Pardiso.hpp"
#include "NewtonSolver.hpp"
#include "R_Precond_Solver.hpp"
#include "Intel_Prec_GMRES.hpp"
#include "Intel_ILU0.hpp"
#include "read_write.h"
//typedef GENSOL::Intel_Pardiso                       LINEARSOLVER;
//typedef NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;
typedef GENSOL::R_Precond_Solver< GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 > >  LINEARSOLVER;
//typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;
#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>

void
read_from_file(const char * filename, std::vector<double> &v)
{
	double tmp;
	std::ifstream strm(filename);
	if (strm.good())
	{
		for (std::size_t i = 0; i < v.size(); ++i)
		if (strm >> tmp) v[i] = tmp;
	}
}

void
dump_solution(const char * filename, const DiscreteProblem::StateVector & v, const double& t)
{
	std::ofstream strm(filename, std::ios_base::app);
	/*strm << "Po\tSw\tSg\tSo\tRso\n";*/
	strm << "time = " << t << std::endl;
	for (std::size_t i = 0; i < v.size(); ++i)
		strm << v[i].Po.value() << "\t"
		<< v[i].Sw.value() << "\t"
		<< v[i].Sg.value() << "\t"
		<< 1.0 - v[i].Sw.value() - v[i].Sg.value() << "\t"
		<< v[i].Rso.value() << std::endl;
	strm.close();
}

void
dump_field(const char * filename,
const std::vector<double> &phi,
const std::vector<double> &kx,
const std::vector<double> &ky,
const std::vector<double> &kz)
{
	std::ofstream strm(filename);
	for (std::size_t i = 0; i < phi.size(); ++i)
		strm << phi[i] << "\t"
		<< kx[i] << "\t"
		<< ky[i] << "\t"
		<< kz[i] << std::endl;
	strm.close();
}

int main()
{

	const std::size_t MAX_NLNITER = 12;
	const double      DT_INIT = 1.0;
	const double      DT_CUT = 0.5;
	const double      DT_GROW = 2.0;
	//const double      T_FINAL     = 365.0;

	const std::size_t NX = 250, NY = 1, NZ = 1;
	const double      LX = 25.0, LY = 50.0, LZ = 200;
	std::vector<double> vPORO = readSingleValue("./input/_poro.dat");
	std::vector<double> vKX = readSingleValue("./input/_permX.dat");
	std::vector<double> vKY = vKX;
	std::vector<double> vKZ = vKX;
	//std::vector<double> vPORO ( NX*NY*NZ, 0.2 );
	//std::vector<double> vKX ( NX*NY*NZ, 20 );
	//std::vector<double> vKY ( NX*NY*NZ, 20 );
	//std::vector<double> vKZ ( NX*NY*NZ, 4 );
	//   read_from_file( "./input/permx.dat", vKX  );   
	//   read_from_file( "./input/permy.dat", vKY  );   
	//   read_from_file( "./input/permz.dat", vKZ  );   
	//   read_from_file( "./input/spe_phi.dat", vPORO  );   
	//   dump_field( "./output/field.out", vPORO, vKX, vKY, vKZ );

	const double OPT_NLNITER = (3 < MAX_NLNITER ? MAX_NLNITER : 3);
	DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
	LINEARSOLVER lnsolver(model.max_num_eqns(),
		model.max_num_nnz());

	STDN newton(model, MAX_NLNITER, 1);

	DiscreteProblem::StateVector uOld, uNew;
	model.initialize_state(uOld);
	double DT = DT_INIT;
	double time = 0.0;

	// read time step and calculate acculmulate time steps
	std::vector<double> tStep = readSingleValue("./input/_tstep.dat");
	std::vector<double> cumTime;
	cumTime.resize(tStep.size());
	std::partial_sum(tStep.begin(), tStep.end(), cumTime.begin());


	for (size_t i = 0; i < cumTime.size(); ++i){
		double T_TARGET = cumTime[i];
		while (time < T_TARGET)
		{
			// check if DT go over the T_TARGET(the target time step we want to output
			if ((time + DT)>T_TARGET){
				DT = T_TARGET - time;
			}

			// run simulation
			uNew = uOld;
			STDN::report_t stdsmry = newton.solve_timestep(uNew, uOld, DT, model, lnsolver);
			if (stdsmry.is_converged)
			{
				uOld = uNew;
				time += DT;
				if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
				std::cout << "CONVERGED t = " << time << " days" << std::endl;
			}
			else
			{
				DT *= DT_CUT;
				std::cout << "FAILED " << std::endl;
			}
		}
		dump_solution("./output/results.out", uNew, time);
	}


	return -1;
}
