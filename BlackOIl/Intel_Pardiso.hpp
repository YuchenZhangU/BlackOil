#ifndef __INTEL_PARDISO_HPP__
#define __INTEL_PARDISO_HPP__

#include <cstring>
#include <vector>
#include "mkl_pardiso.h"
#include "CSR_Matrix.hpp"

namespace GENSOL {

	class Intel_Pardiso
	{
	public:
		typedef double                                     double_type;
		typedef MKL_INT                                    int_type;

		typedef CSR_Matrix< double_type, int_type >        A_type;
		typedef std::vector< double_type >                 x_type;

	public:

		Intel_Pardiso(std::size_t _N, std::size_t _NNZ_MAX) :
			pt(),
			maxfct(1),
			mnum(1),
			mtype(11),
			phase(13),
			n(_N),
			idummy(0),
			nrhs(1),
			iparm(),
			msglvl(0),
			error(0)
		{
			std::memset(pt, 0, sizeof(pt));
			std::memset(iparm, 0, sizeof(iparm));
			iparm[0] = 1;        /* No default values for solver */
			iparm[1] = 3;        /* Fill-in reducing ordering */
			iparm[3] = 0;        /* No iterative-direct algorithm */
			iparm[4] = 0;        /* No user fill-in reducing permutation */
			iparm[5] = 0;        /* Write solution on x */
			iparm[6] = 0;        /* No output of iterative refinement progress */
			iparm[7] = 2;        /* Maximum number of iterative refinement steps */
			iparm[9] = 13;       /* Perturb the pivot elements with 1E-13 */
			iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
			iparm[11] = 0;        /* Do not solve transposed matrix*/
			iparm[12] = 1;        /* Maximum weighted matching (default for nonsymmetric) */
			iparm[17] = 1;        /* No Output: Number of non-zero values in the factor LU */
			iparm[18] = 1;        /* No Output: Mflops for LU factorization */
			iparm[19] = 0;        /* No Output: Numbers of CG Iterations */
			iparm[34] = 1;        /* Zero based indexing */
		}

		~Intel_Pardiso()
		{
			phase = -1;
			double dsentinel;
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
				&n, &dsentinel, &idummy, &idummy, &idummy,
				&nrhs, iparm, &msglvl, &dsentinel, &dsentinel, &error);
		}

		static int_type offset() { return 0; }

		int solve(A_type & _A, x_type & _x_, x_type & _b)
		{
			n = _A.N();
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
				&n, _A.value().data(), _A.rowptr().data(), _A.colind().data(),
				&idummy, &nrhs, iparm, &msglvl, _b.data(), _x_.data(), &error);
			return error;
		}

	private:
		void     * pt[64];
		int_type   maxfct;
		int_type   mnum;
		int_type   mtype;
		int_type   phase;
		int_type   n;
		int_type   idummy;
		int_type   nrhs;
		int_type   iparm[64];
		int_type   msglvl;
		int_type   error;
	};

};


#endif
