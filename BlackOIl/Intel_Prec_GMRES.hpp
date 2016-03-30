#pragma once

#include "mkl_rci.h"

#include <vector>
#include "CSR_Matrix.hpp"

namespace GENSOL{

  template< typename T_PRECOND >
  class Intel_Prec_GMRES
  {

  public:
    typedef double                                     double_type;
    typedef MKL_INT                                    int_type;
    typedef CSR_Matrix< double_type, int_type >        A_type;
    typedef std::vector< double_type >                 x_type;
   
  public:
    //.............................  LIFECYCLE  ..........................//
    Intel_Prec_GMRES( std::size_t _N_MAX, std::size_t _NNZ_MAX ) : 
      mPreconditioner( _N_MAX, _NNZ_MAX ),
      N( _N_MAX ),
      RCIreq( ),
      CVAR('N')
    { 
      if ( _N_MAX > 100000 ) 
	{
	  do_restart = true;
	  double NTMP = _N_MAX;
	  NRESTART   = 2 * ( sqrt( NTMP*(NTMP+4) + 15055967.56 ) - NTMP ) - 4.5;
	  NRESTART   = (NRESTART < 3 ? 3 : NRESTART );
	}
      else
	{
	  NRESTART = ( N > 150 ? 150 : N );
	}
      tmp_vec.resize( (2*NRESTART+1)*N+NRESTART*(NRESTART+9)/2+1  );
    };
   
    ~Intel_Prec_GMRES( )
    {
    };

    static int_type offset( ) { return 1; }
    
    //.............................  OPERATORS  .........................//
    int solve( A_type &_A, x_type &_x, x_type & _b, const EllipticInfo &ainfo=EllipticInfo() )
    {
      N = _A.N();

      for ( std::size_t i=0; i<_x.size(); ++i ) _x[i]= _b[i];

      int result = initialize_solver( _x, _b );
      if (result >= 0 )
	{
	  result = mPreconditioner.setup_A( _A  );
	  if ( result >= 0 )
	    {
	      bool keep_going = true;
	      while (keep_going)
		{

		  dfgmres (&N, _x.data(), _b.data(), &RCIreq, ipar, dpar, tmp_vec.data() );
		  switch (RCIreq)
		    {
		    case 0:
		      keep_going = false;
		      break;
		    case 1: 
		      mkl_dcsrgemv (&CVAR, 
				    &N, 
				    _A.value().data(), _A.rowptr().data(), _A.colind().data(), 
				    &tmp_vec[ipar[21] - 1 ], &tmp_vec[ipar[22] - 1 ]);
		      break;
		    case 3:
		      mPreconditioner.apply( &tmp_vec[ipar[22]-1], &tmp_vec[ipar[21]-1] );
		      break;
		    default:
		      keep_going = false;
		      result = -1;
		      break;
		    } //switch
		} // while
	    }// ilu0 ok
	} // gmres ok

      int_type itercount;
      ipar[12] = 0;
      dfgmres_get (&N, _x.data(), _b.data(), &RCIreq, ipar, dpar, tmp_vec.data(), &itercount);
      std::cout << "N LN ITER = " << itercount << std::endl;
      return result;
    }

  protected:

    int initialize_solver( x_type &_x, x_type & _b  )
    {
      int result = -1;
      dfgmres_init ( &N, _x.data(), _b.data(), &RCIreq, ipar, dpar, tmp_vec.data() );
      if (RCIreq == 0)
	{
	  setup_options( );
	  dfgmres_check (&N, _x.data(), _b.data(), &RCIreq, ipar, dpar, tmp_vec.data() );
	  if (RCIreq != -1100) result = 0;
	}
      return result;
    }

    void setup_options( )
    {
      ipar[ 5] = 0; // no error output
      ipar[ 6] = 0; // no error output

      ipar[ 7] = 1; // auto test for max iter
      ipar[ 8] = 1; // auto test for residual criterion
      ipar[ 9] = 0; // no manual tests
      ipar[10] = 1; // use a preconditioner
      ipar[11] = 1; // auto test for zero orthogonal vector norm
      ipar[14] = NRESTART;
    }
      
  private:
    T_PRECOND     mPreconditioner;
    int_type      N;
    int_type      RCIreq;
    char          CVAR;
    x_type        tmp_vec;
    int_type      ipar [ 128 ];
    double_type   dpar [ 128 ];
    std::size_t   NRESTART;
    bool          do_restart;
  };

};

