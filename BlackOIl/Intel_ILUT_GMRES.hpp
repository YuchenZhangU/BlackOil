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
   class Intel_ILUT_GMRES
   {
   public:
      typedef double                                     double_type;
      typedef MKL_INT                                    int_type;
      typedef CSR_Matrix< double_type, int_type >        A_type;
      typedef fastl::pod_vector_unbounded< double_type > x_type;
   
   public:
      //.............................  LIFECYCLE  ..........................//
      Intel_ILUT_GMRES( std::size_t _M_MAX, std::size_t _N_MAX, std::size_t _NNZ_MAX ) : 
	 N( _N_MAX ),
	 MAXFIL( 30 ),
	 ILUTOL( 1.0e-6 ),
	 p_tmp_buffer( new double_type [ _N_MAX*(2*_N_MAX+1)+(_N_MAX*(_N_MAX+9))/2 + 1 ] ),
	 p_trvec( new double_type [ _N_MAX ] ),
	 bilut( new double_type [ (2*MAXFIL+1)*_NNZ_MAX - MAXFIL*(MAXFIL+1) +1 ] ),
	 ibilut( new MKL_INT [ _NNZ_MAX + 1 ] ),
	 jbilut( new MKL_INT [ (2*MAXFIL+1)*_NNZ_MAX - MAXFIL*(MAXFIL+1) +1 ] ),
	 p_residual( new double_type [ _N_MAX ] )
      {    };
   
      ~Intel_ILUT_GMRES( )
      {
	 delete [] p_residual;
	 delete [] jbilut;
	 delete [] ibilut;
	 delete [] bilut;
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
//	    ipar[14] = 5;        // restarts
	    dpar[0]  = 1.0e-6;   // relative error tolerance for convergence
	    dpar[7]  = 1.0e-6;

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
	       dpar[30] = 10*ILUTOL;
	       dcsrilut (&N, 
			 _A.value().begin(), _A.rowptr().begin(), _A.colind().begin(), 
			 bilut, ibilut, jbilut, &ILUTOL, &MAXFIL, ipar, dpar, &ierr);
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
			   if ( test_residual( _A, _x, _b ) < dpar[7] ) keep_going=false;
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
	 if ( test_residual( _A, _x, _b ) >= dpar[7]  ) result = -1;
	 int_type itercount;
	 ipar[12] = 0;
	 dfgmres_get (&N, _x.begin(), _b.begin(), &RCI_request, ipar, dpar, p_tmp_buffer, &itercount);
	 return result;
      }
   protected:

      void apply_preconditioner( A_type &_A )
      {		
	 char cvar1 = 'L';
	 char cvar = 'N';
	 char cvar2 = 'U';
	 mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N, 
		       bilut,ibilut,jbilut, 
		       &p_tmp_buffer[ipar[21] - 1], p_trvec);
	 cvar1 = 'U';
	 cvar2 = 'N';
	 mkl_dcsrtrsv (&cvar1, &cvar, &cvar2, &N,
		       bilut,ibilut,jbilut, 
		       p_trvec, &p_tmp_buffer[ipar[22] - 1]);
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
      int_type      N;
      MKL_INT       MAXFIL;
      double_type   ILUTOL;
      double_type * p_tmp_buffer;
      double_type * p_trvec;
      double_type * bilut;
      MKL_INT     * ibilut;
      MKL_INT     * jbilut;
      double_type * p_residual;
      x_type      iterate;
      int_type    ipar [ 128  ];
      double_type dpar [ 128  ];
   };

};

