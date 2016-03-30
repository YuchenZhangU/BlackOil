#pragma once

#include "CSR_Matrix.hpp"
#include <vector>

namespace GENSOL{

  template< typename SOL >
  class R_Precond_Solver
  {
  public:
    typedef double                                              double_type;
    typedef int                                                 int_type;

    typedef CSR_Matrix< double_type, int_type >                 A_type;
    typedef std::vector< double_type >                          x_type;
    typedef std::vector< std::size_t >                          i_type;

  public:
    //.............................  LIFECYCLE  ..........................//
    R_Precond_Solver( std::size_t _N_MAX, std::size_t _NNZ_MAX ) : solver( _N_MAX, _NNZ_MAX )
    {
    };
   
    ~R_Precond_Solver( )
    {
    };

    static int_type offset( ) { return 1; }

    //.............................  OPERATORS  .........................//
    int solve( const A_type &_A, x_type &_x, const x_type & _b )
    {
      mA = _A;
      mb = _b;
      N = _A.N();
      calculate_r_norms( );
      apply_scaling ( );
      int ierr = solver.solve( mA, _x, mb );
      return ierr;
    }

  protected:

    void
    apply_scaling ( )
    {
      // apply row and column scaling
      for ( std::size_t r=0; r<N; ++r ) // for each row
	{
	  mb[ r ] /= nrm_row[ r ];

	  const std::size_t inz1 = mA.rowptr()[r]  -mA.offset();
	  const std::size_t inz2 = mA.rowptr()[r+1]-mA.offset();
	  for ( std::size_t inz=inz1; inz<inz2; ++inz ) // for each nonzero in this row
	    {
	      const std::size_t c = mA.colind()[inz] - mA.offset();
	      mA.value()[inz] /= nrm_row[r];
	    }
	}
    }

    void
    calculate_r_norms( )
    {
      // initialize arrays
      nrm_row.resize( N );
      for ( std::size_t r=0; r<N; ++r ) 
	nrm_row[r] = 0.0;

      // calculate row and column infinity norms
      for ( std::size_t r=0; r<N; ++r ) // for each row
	{
	  double max_val = 0.0;
	  const std::size_t inz1 = mA.rowptr()[r]  -mA.offset();
	  const std::size_t inz2 = mA.rowptr()[r+1]-mA.offset();
	  for ( std::size_t inz=inz1; inz<inz2; ++inz ) // for each nonzero in this row
	    {
	      const double      v = mA.value()[inz];
	      if ( std::fabs(v) > max_val )    max_val = std::fabs(v); // build row inf norm
	    }
	  nrm_row[r] = max_val;
	}
    }

  private:
    SOL                 solver;
    A_type              mA;
    x_type              mb;
    x_type              nrm_row;
    std::size_t         N;
  };

};

