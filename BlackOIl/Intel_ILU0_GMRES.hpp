#pragma once

#include <limits>
#include <iostream>
#include <cmath>
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_rci.h"
#include "mkl_service.h"
#include "fastl/containers/pod_vector_unbounded.hpp"
#include "CSR_Matrix.hpp"
namespace GENSOL{
   class Intel_ILU0_GMRES
   {
   public:
      typedef double                                     double_type;
      typedef MKL_INT                                    int_type;
      typedef CSR_Matrix< double_type, int_type >        A_type;
      typedef fastl::pod_vector_unbounded< double_type > x_type;
   
   public:
      //.............................  LIFECYCLE  ..........................//
      Intel_ILU0_GMRES( std::size_t _M_MAX, std::size_t _N_MAX, std::size_t _NNZ_MAX ) : 
	 p_tmp_buffer( new double_type [ _N_MAX*(2*_N_MAX+1)+(_N_MAX*(_N_MAX+9))/2 + 1 ] ),
	 p_trvec( new double_type [ _N_MAX ] ),
	 p_bilu0( new double_type [ _NNZ_MAX ] ),
	 p_residual( new double_type [ _N_MAX ] ),
	 TOL( 1.0e-3 ),
	 N( _N_MAX )
      {    };
   
      ~Intel_ILU0_GMRES( )
      {
	 delete [] p_residual;
	 delete [] p_bilu0;
	 delete [] p_trvec;
	 delete [] p_tmp_buffer;	 
      };

      static int_type offset( ) { return 1; }

      //.............................  OPERATORS  .........................//
      int solve( A_type &_A, x_type &_x, x_type & _b )
      {
	 int result = 1;

	 // initialize extra copies
	 iterate = _b;
	 for ( std::size_t i=0; i<_x.size(); ++i ) _x[i]= 0.0;

	 // Initialize solver
	 int_type RCI_request;
	 dfgmres_init (&N, _x.begin(), _b.begin(), &RCI_request, ipar, dpar, p_tmp_buffer);
	 if (RCI_request!=0) 
	 {
	    result = -1;
	 }
	 else 
	 {
	    ipar[7] = 1;         // auto check maximum number of iterations
	    ipar[10] = 1;        // use preconditioner
	    ipar[14] = 5;        // restarts
	    dpar[0]  = 1.0e-3;   // relative error tolerance for convergence
	    dpar[7]  = TOL;

	    dfgmres_check (&N, _x.begin(), _b.begin(), &RCI_request, ipar, dpar,p_tmp_buffer);
	    if (RCI_request==-1100) 
	    {
	       result -1;
	    }
	    else
	    {
	       // Setup IL0
	       int_type ierr = 0;
	       ipar[30] = 1;
	       dpar[30] = 1.0e-10;
	       dpar[31] = 1.0e-8;
	       dcsrilu0 (&N, 
			 _A.value().begin(), _A.rowptr().begin(), _A.colind().begin(), 
			 p_bilu0, ipar, dpar, &ierr);
	       if (ierr != 0) 
	       {
		  result = -1;
	       }
	       else
	       {
		  bool keep_going = true;
		  while (keep_going)
		  {
		     dfgmres (&N, _x.begin(), _b.begin(), &RCI_request, ipar, dpar, p_tmp_buffer);
		     switch (RCI_request)
		     {
			case 0:
			   keep_going = false;
			   break;
			case 1: 
			   multiply_with_temp( _A );
			   break;
			case 2:
			   if ( test_residual( _A, _x, _b ) < TOL ) keep_going=false;
			   break;
			case 3:
			   apply_preconditioner( _A );
			   break;
			case 4:
			   if (dpar[6] < std::sqrt(std::numeric_limits<double>::epsilon( ))) keep_going = false;
			   break;
			default:
			   keep_going = false;
			   result = -1;
			   break;
		     } //switch
		  } // while
	       }// ilu0 ok
	    } // gmres check ok
	 } // gmres init ok
//	 if ( !test_residual( _A, _x, _b ) ) result = -1;
	 std::cout << test_residual( _A, _x, _b ) << std::endl;
	 int_type itercount;
	 ipar[12] = 0;
	 dfgmres_get (&N, _x.begin(), _b.begin(), &RCI_request, ipar, dpar, p_tmp_buffer, &itercount);
	 std::cout << "NITER = " << itercount << std::endl;
	 if ( result != -1 )
	    std::cout << "CONVERGED" << std::endl;
	 else
	    std::cout << "FAILED" << std::endl;

	 return result;
      }
   protected:

      void apply_preconditioner( A_type &_A )
      {		
	 char cvar1 = 'L';
	 char cvar = 'N';
	 char cvar2 = 'U';
	 mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N, p_bilu0, _A.rowptr().begin(), _A.colind().begin(),
		       &p_tmp_buffer[ipar[21] - 1], p_trvec);
	 cvar1 = 'U';
	 cvar2 = 'N';
	 mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N, p_bilu0, 
		       _A.rowptr().begin(), _A.colind().begin(), p_trvec, &p_tmp_buffer[ipar[22] - 1]);
      }
      
      void multiply_with_temp( A_type &_A )
      {
	 char cvar = 'N';
	 mkl_dcsrgemv (&cvar, 
		       &N, 
		       _A.value().begin(), _A.rowptr().begin(), _A.colind().begin(), 
		       &p_tmp_buffer[ipar[21] - 1], &p_tmp_buffer[ipar[22] - 1]);
      }

      double test_residual( A_type &_A, x_type &_x, x_type & _b )
      {
	 // uses: N, ipar, dpar, iterate, p_tmp_buffer, p_residual, TOL

	 // local variables
	 int_type RCI_request;
	 int_type itercount;

	 // put current iterate in b
	 int origpar = ipar[12];
	 ipar[12] = 1;
	 dfgmres_get (&N, _x.begin(), iterate.begin(), &RCI_request, ipar, dpar, p_tmp_buffer, &itercount);
	 ipar[12] = origpar;

	 // p_residual = A . iterate
	 char cvar = 'N';
	 mkl_dcsrgemv (&cvar, 
		       &N, 
		       _A.value().begin(), _A.rowptr().begin(), _A.colind().begin(), 
		       iterate.begin(), 
		       p_residual);

	 // p_residual = A . iterate - b and compute norm ||p_residual||
	 double nrm = 0.0;
	 for ( std::size_t i=0; i<N; ++i ) 
	 {
	    p_residual[i] -= _b[i];
	    nrm += p_residual[i] * p_residual[i];
	 }
	 nrm = std::sqrt(nrm);
	 return nrm;
      }
      
   private:
      double_type * p_tmp_buffer;
      double_type * p_trvec;
      double_type * p_bilu0;
      double_type * p_residual;
      const double TOL;
      int_type     N;
      x_type      iterate;
      int_type    ipar [ 128  ];
      double_type dpar [ 128  ];
   };

};

