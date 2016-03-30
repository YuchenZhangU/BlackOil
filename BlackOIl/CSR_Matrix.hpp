// -------------------------------- *- C++ -* ------------------------------ //
// -------------------------------* TEMPLATE *------------------------------ //
// ------------------------------------------------------------------------- //
//! \file  CSR_Matrix.hpp
//! \brief Compressed Sparse Row Storage Matrix datatype
//
// ------------------------------------------------------------------------- //
/** *************************************************************************
 *  \author     : Rami M. Younis
 *  \created    : 07/22/2013
 *  \revised    : 
 *  \warning    : 
 *  Target MCU  : Generic
 *  Editor TABs : 3
 *  Editor      : emacs
 *  Auto Style  : ellemtel
 ** ************************************************************************ */
/** **************************************************************************
 * Copyright (c) 2013, all rights reserved
 * FUture Reservoir Simulation Systems & Technology
 * McDougall School of Petroleum Engineering
 * The University of Tulsa, Tulsa, Oklahoma 74104, USA
 ** ************************************************************************ */
#ifndef __CSRMATRIX_HPP_INCLUDED_
#define __CSRMATRIX_HPP_INCLUDED_

#include <vector>
#include <ostream>

namespace GENSOL { // -------------------------------------------- BEGIN NAMESPACE 

  class EllipticInfo
  {
  public: 
    std::size_t elliptic_var_id;
    std::size_t n_vars_per_block;
    std::size_t n_blocks;
    std::size_t n_unblocked_vars;
  };

  // ----------------------------  CSR_Matrix  ------------------------------ //
  /** \class CSR_Matrix
   *  
   *  Compressed Sparse Row Storage Matrix datatype
   **/
  // ------------------------------------------------------------------------- //
  template< typename __T, typename __I >
  class CSR_Matrix
  {
  public:
    typedef __T *                             T_pointer;
    typedef __I *                             I_pointer;

  private:
    typedef std::vector<__T>                  T_collection;
    typedef std::vector<__I>                  I_collection;

  public:
    //.............................  LIFECYCLE  ..........................//
    CSR_Matrix(  ) : m_val( ), m_colind( ), m_rowptr( ) {  }

    CSR_Matrix(  __I _N, __I _NNZ ) { resize(_N,_NNZ); }

      
    //...........................  ASSIGNMENT  ...........................//      
    
    //.............................  ACCESS  .............................//
    __I N() const { return mN; }
    __I NNZ() const { return mNNZ; }
    __I offset() const { return m_rowptr[0]; }

    const T_collection & value()  const { return m_val; }
    const I_collection & colind() const { return m_colind; }
    const I_collection & rowptr() const { return m_rowptr; }
    __T                  operator()( std::size_t i, std::size_t j ) const
    {
      const __I offset = m_rowptr[0];
      double    value  = 0;
      if ( (i<mN) && (j<mN) ) 
	{
	  const __I r1 = m_rowptr[i] - offset;
	  const __I r2 = m_rowptr[i+1] - offset - 1;
	  if ( ( j >= (m_colind[r1]-offset) ) && ( j <= ( m_colind[r2]-offset ) ) )
	    {
	      __I r = r1;
	      while ( j > (m_colind[r]-offset) ) ++r;
	      if ( j == ( m_colind[r]-offset) ) value = m_val[r];
	    }
	}
      return value;
    };

    T_collection & value()  { return m_val; }
    I_collection & colind() { return m_colind; }
    I_collection & rowptr() { return m_rowptr; }   

    //.............................  OPERATORS  .........................//
    void resize( __I _N, __I _NNZ )
    {
      mN = _N; mNNZ = _NNZ;
      m_val.resize( mNNZ );
      m_colind.resize( mNNZ );
      m_rowptr.resize( mN + 1 );
    }
   
    void check_size( )
    {
      mN   = m_rowptr.size()-1;
      mNNZ = m_val.size();
    }

  public:
    EllipticInfo ainfo;
  private:
    __I mN;
    __I mNNZ;
    T_collection m_val;
    I_collection m_colind;
    I_collection m_rowptr;
  }; // CSR_Matrix

  template< typename A, typename B >
  std::ostream & spy ( std::ostream & ostr, 
		       const GENSOL::CSR_Matrix<A,B> & _out )
  {
    for ( std::size_t r = 0; r < _out.N( ); ++r )
      {
	ostr << "[";
	const std::size_t I1  = _out.rowptr()[r] - 1;
	const std::size_t I2  = _out.rowptr()[r+1] - 1;
	std::size_t cc  = 0;
	for ( std::size_t nz = I1; nz < I2; ++nz )
	  {
	    const std::size_t c = _out.colind()[nz];
	    for ( ; cc < c; ++cc ) ostr << " ";
	    if (_out.value()[nz]!=0.0) ostr << "1";
	    ++cc;
	  }
	ostr << "]" << std::endl;
      }
    return ostr;
  };

}; // ---------------------------------------------------------- END NAMESPACE

template< typename A, typename B >
std::ostream & operator << ( std::ostream & ostr, 
			     const GENSOL::CSR_Matrix<A,B> & _out )
{
  ostr << _out.N() << " X " << _out.N() << "\t NNZ = " << _out.NNZ() << std::endl;
  ostr << "VAL = [";
  for (std::size_t i = 0; i<_out.NNZ(); ++i )
    ostr << _out.value()[i] << ", ";
  ostr << " ]" << std::endl;
  ostr << "COL = [";
  for (std::size_t i = 0; i<_out.NNZ(); ++i )
    ostr << _out.colind()[i] << ", ";
  ostr << " ]" << std::endl;
  ostr << "ROW = [";
  for (std::size_t i = 0; i<_out.N()+1; ++i )
    ostr << _out.rowptr()[i] << ", ";
  ostr << " ]";
  return ostr;
}



#endif // __CSR_Matrix_HPP_INCLUDED_

//---------------------------------------------------------------------------//
//                           EOF CSR_Matrix.hpp
//---------------------------------------------------------------------------//

