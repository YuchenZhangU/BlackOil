#pragma once

#include "mkl.h"
#include "mkl_rci.h"
#include "CSR_Matrix.hpp"
#include <vector>

namespace GENSOL{

  class Intel_ILU0
  {
  public:
    typedef double                                     double_type;
    typedef MKL_INT                                    int_type;

    typedef CSR_Matrix< double_type, int_type >        A_type;
    typedef std::vector< double_type >                 x_type;
   
  public:
    //.............................  LIFECYCLE  ..........................//
    Intel_ILU0( std::size_t _N_MAX, std::size_t _NNZ_MAX ) :
      bilu0( _NNZ_MAX ),
      bilu0_rowptr( _N_MAX + 1 ),
      bilu0_colind( _NNZ_MAX ),
      tmp_vec( _N_MAX ),
      mN( _N_MAX )
    {
      for ( std::size_t i=0; i<128; ++i )
	{
	  ipar[i] = 0;
	  dpar[i] = 0.0;
	}
      ipar[ 1] = 6;       // if error messages allowed display to screen
      ipar[ 5] = 0;       // no error messages
      ipar[30] = 1;       // check for zero pivot
      dpar[30] = 1.0e-16; // this is condered zero
      dpar[31] = 1.0e-10; // replace zero with this number
    };
   
    ~Intel_ILU0( )
    {
    };

    static int_type offset( ) { return 1; }

    //.............................  OPERATORS  .........................//
    int setup_A( A_type &_A )
    {
      int ierr = 0;

      mN       = _A.N();

      tmp_vec.resize( mN );
      bilu0.resize( _A.value().size()  );
      bilu0_rowptr = _A.rowptr();
      bilu0_colind = _A.colind();

      dcsrilu0 (&mN, 
		_A.value().data(), bilu0_rowptr.data(), bilu0_colind.data(), 
		bilu0.data(), ipar, dpar, &ierr);

      if (ierr != 0) 
	{
	  std::cout << "dcsrilu0 gave error= " << ierr <<std::endl;
	  ierr = -1;
	}

      return ierr;
    }

    int apply( double_type *_x, double_type *_b )
    {
      int  ierr = 0;
      char cvar = 'N'; // no transpose

      char cvar1 = 'L'; // invert lower triangular part
      char cvar2 = 'U'; // assume unit lower traingular
      mkl_dcsrtrsv ( &cvar1, &cvar, &cvar2, 
		     &mN, bilu0.data(), bilu0_rowptr.data(), bilu0_colind.data(),
		     _b, tmp_vec.data());

      cvar1 = 'U'; // take upper traingular part
      cvar2 = 'N'; // it is not unit triangular
      mkl_dcsrtrsv ( &cvar1, &cvar, &cvar2, 
		     &mN, bilu0.data(), bilu0_rowptr.data(), bilu0_colind.data(),
		     tmp_vec.data(), _x );
      return ierr;
    }

    int solve( A_type &_A, x_type &_x, x_type & _b )
    {
      int ierr = setup_A( _A );
      apply( _x.data( ), _b.data( ) );
      return ierr;
    }

  private:
    std::vector<double_type> bilu0;
    std::vector<int_type>    bilu0_rowptr;
    std::vector<int_type>    bilu0_colind;
    std::vector<double_type> tmp_vec;
    int_type      mN;
    int_type      ipar [ 128  ];
    double_type   dpar [ 128  ];
  };

};

