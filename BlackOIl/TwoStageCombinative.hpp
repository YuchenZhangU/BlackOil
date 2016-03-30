#pragma once

#include "Intel_ILU0.hpp"
#include "SAMGwrapper.hpp"
#include "CSR_Matrix.hpp"
#include <vector>

namespace GENSOL{

  template< typename M1=GENSOL::SAMGwrapper, typename M2=GENSOL::Intel_ILU0 >
  class TwoStageCombinative
  {
  public:
    typedef double                                     double_type;
    typedef int                                        int_type;

    typedef CSR_Matrix< double_type, int_type >        A_type;
    typedef std::vector< double_type >                 x_type;
   
  public:
    //.............................  LIFECYCLE  ..........................//
    TwoStageCombinative( std::size_t _N_MAX, std::size_t _NNZ_MAX ) : 
      m_prec1( _N_MAX, _NNZ_MAX ), m_prec2( _N_MAX, _NNZ_MAX )
    {
    };
   
    static int_type offset( ) { return 1; }

    //.............................  OPERATORS  .........................//
    int setup_A( A_type &_A  )
    {
      m_A     =_A;
      int err = 0;
      err     = restrict_operator( );
      err     = m_prec1.setup_A( m_Ap );
      err     = m_prec2.setup_A( _A );
      return err;
    }

    int apply( double_type *_x, double_type *_b )
    {
      int  ierr = 0;

      restrict_rhs( _b );
      m_xp.resize( m_A.ainfo.n_blocks);
      for ( std::size_t blk = 0 ; blk < m_A.ainfo.n_blocks; ++blk )
	  m_xp[blk] = _x[ blk * m_A.ainfo.n_vars_per_block + m_A.ainfo.elliptic_var_id ];
      m_prec1.apply( m_xp.data(), m_bp.data() );

      project_xp( );
      calculate_residual( _b );

      std::size_t N =  m_A.N();
      m_xs.resize( N );
      m_prec2.apply( m_xs.data(), m_bs.data() );

      for ( std::size_t i=0; i < N; ++i ) _x[i] = m_xs[i] + m_xp_long[i];

      return ierr;
    }

  protected:

    int
    restrict_operator( )
    {
      m_Ap.value().clear();
      m_Ap.colind().clear();
      m_Ap.rowptr().clear();
      m_Ap.rowptr().push_back( m_A.offset() );
      for ( std::size_t blk = 0 ; blk < m_A.ainfo.n_blocks; ++blk )
	{
	  std::size_t r = blk * m_A.ainfo.n_vars_per_block + m_A.ainfo.elliptic_var_id;
	  std::size_t nnz_in_row = 0;
	  for( std::size_t inz = m_A.rowptr()[r]-m_A.offset(); 
	       inz < m_A.rowptr()[r+1]-m_A.offset();
	       ++inz )
	    {
	      std::size_t var_id = (m_A.colind()[inz]-m_A.offset()) % m_A.ainfo.n_vars_per_block;
	      if ( (var_id == m_A.ainfo.elliptic_var_id) &&
		   ((m_A.colind()[inz]-m_A.offset())< (m_A.ainfo.n_blocks*m_A.ainfo.n_vars_per_block) ) )
		{
		  // extract entry
		  ++nnz_in_row;
		  std::size_t reduced_i = std::floor(m_A.colind()[inz]/m_A.ainfo.n_vars_per_block ) + m_A.offset();
		  m_Ap.value().push_back( m_A.value()[inz] );
		  m_Ap.colind().push_back( reduced_i );
		}
	    }
	  m_Ap.rowptr().push_back( m_Ap.rowptr().back() + nnz_in_row );
	}
      m_Ap.check_size( );
      m_Ap.ainfo.elliptic_var_id  = 0;
      m_Ap.ainfo.n_vars_per_block = 1;
      m_Ap.ainfo.n_blocks         = m_Ap.N();
      m_Ap.ainfo.n_unblocked_vars = 0;

      return 0;
    }

    int 
    restrict_rhs( double_type *_b )
    {
      m_bp.resize( m_A.ainfo.n_blocks );
      for ( std::size_t blk = 0 ; blk < m_A.ainfo.n_blocks; ++blk )
	{
	  m_bp[blk] = _b[ blk * m_A.ainfo.n_vars_per_block + m_A.ainfo.elliptic_var_id ];
	}
      return 0;
    }

    void
    project_xp( )
    {
      std::size_t N =  m_A.N();
      m_xp_long.resize( N );
      for ( std::size_t i=0; i<N; ++i) m_xp_long[i] = 0.0;
      for ( std::size_t i=0; i<m_A.ainfo.n_blocks; ++i) 
	m_xp_long[i*m_A.ainfo.n_vars_per_block+m_A.ainfo.elliptic_var_id] = m_xp[i];
    }

    void
    calculate_residual( double_type *_b  )
    {
      char TRANS   = 'N';      // no transpose
      int  NA      = m_A.N();  // dimension of square matrices
      m_bs.resize( NA );
      mkl_dcsrgemv (&TRANS, 
		    &NA, 
		    m_A.value().data(), m_A.rowptr().data(), m_A.colind().data(), 
		    m_xp_long.data(), 
		    m_bs.data() );
      for ( std::size_t i=0; i<NA; ++i) 
	{
	  m_bs[i] -= _b[i];
	  m_bs[i] *= -1.0;
	}
    }

  private:
    M1                   m_prec1;
    M2                   m_prec2;
    A_type               m_Ap;
    A_type               m_A;
    x_type               m_bp;
    x_type               m_xp;
    x_type               m_xp_long;
    x_type               m_bs;
    x_type               m_xs;
  };

};

