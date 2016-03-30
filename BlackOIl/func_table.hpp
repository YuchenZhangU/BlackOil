// -------------------------------- *- C++ -* ------------------------------ //
// -------------------------------* TEMPLATE *------------------------------ //
// ------------------------------------------------------------------------- //
//! \file  func_table.hpp
//! \brief Lookup table for function approximations. 1 key, 1 value column
//
// ------------------------------------------------------------------------- //
/** *************************************************************************
 *  \author     : Rami M. Younis
 *  \created    : 03/18/2007
 *  \revised    : 08/23/2007
 *  \version    : 0.0.0
 *  \warning    : 
 *  Target MCU  : Generic
 *  Editor TABs : 3
 *  Editor      : emacs
 *  Auto Style  : stroustrup
 ** ************************************************************************ */
/** **************************************************************************
 * Copyright (c) 2006, 2007 by Rami Younis,
 * Stanford University Petroleum Research Institute - B,
 * Stanford University, California 94305, USA
 * ryounis_AT_stanford dot edu
 *
 * Permission to use, copy, modify, and/or distribute this software and/or its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice and this permission notice
 * appear in all copies, modifications, and distributions. Otherwise, all
 * rights are reserved.
 *
 * This software is distributed in the hope that it will be useful, but no
 * representations are made about the suitability of this software for any
 * purpose. It is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  It is
 * provided "as is" without express or implied warranty.
 ** ************************************************************************ */
#ifndef __FUNC_TABLE_HPP_INCLUDED_
#define __FUNC_TABLE_HPP_INCLUDED_

#include <cstddef> // provides size_t difference_t
#include <algorithm> // provides is_sorted
#include "interpolators.hpp"
#include <vector>
#include "fastl/containers/pod_array.hpp"

#define __restrict__

namespace fastl { // ----------------------------------------- BEGIN NAMESPACE 

   // ---------------------------  sorted_key  ----------------------------- //
   /** \class sorted_key
    * 
    *  Encapsulates the key data (domain) of a univariate tabular function.
    *  The key data is provided without an equi-spacing requirement. 
    *  Subsequently it is initialized in sorted order (ascending), and find
    *  functionality uses binary search.
    *  __T must be built-in numeric
    *  __C must have
    *      ctor( size_t, __T * )
    *      __T * begin()
    *      __T * end()
    *      __T  opertor [] ()
    *      size_t size()
    **/
   // ---------------------------------------------------------------------- //
   template< typename __T = float > 

   class sorted_key
   {
   public:
      sorted_key ( ) : mc_data ( ) { };

      sorted_key ( const std::size_t _dim, 
		   const __T * __restrict__ _p_data )
	 : mc_data( _p_data, _p_data + _dim )
         { 
	 }

      void
      assign ( const std::size_t _dim, 
	       const __T * __restrict__ _p_data )
      { 
	 mc_data.clear();
	 mc_data.assign( _p_data, _p_data + _dim );
      }

      std::size_t
      size ( ) const
	 { return mc_data.size( ); }

      __T
      operator[] ( const std::size_t _i ) const
	 {
	    return mc_data[_i];
	 }

      __T &
      operator[] ( const std::size_t _i )
	 {
	    return mc_data[_i];
	 }

      void push_back( __T v )
      {
	 mc_data.push_back( v );
      }

      template< typename __T_ans >
      std::size_t
      operator( ) ( const __T_ans & __restrict__ _x ) const
	 {
	    const __T x( _x ); // note use of conversion operator
	    std::size_t r_ind;
	    if ( mc_data.front( ) >= x ) {
	       r_ind = 0;
	    } else if ( mc_data.back( ) <= x ) {
	       r_ind = mc_data.size( );
	    } else {
	       r_ind = interval_find( mc_data.begin( ), (mc_data.end( )-1), x );
	    }
	    return r_ind;
	 }

   protected:
      
      std::size_t
      interval_find ( typename std::vector<__T>::const_iterator _left,
		      typename std::vector<__T>::const_iterator _right,
		      const __T & __restrict__ _x ) const
	 {
	    std::ptrdiff_t len = _right - _left;
	    
	    while ( len>1 )
	    {
	       typename std::vector<__T>::const_iterator mid = _left + (len >> 1);
	       const unsigned short jmp = ((*mid<_x)<<1) + (*mid>_x);
	       switch ( jmp )
	       {
	       case 0:
		  _right = _left = mid;
		  break;
	       case 1:
		  _right = mid;
		  break;
	       default:
		  _left = mid;
		  break;
	       }
	       len = _right - _left;
	    }
	    return static_cast<std::size_t>( _right - mc_data.begin( ) );
	 }
      
   private:
      std::vector<__T> mc_data;
   };
      

   // ------------------------  equidistant_key  --------------------------- //
   /** \class equidistant_key
    * 
    *  Encapsulates the key data (domain) of a univariate tabular function.
    *  The key data is assumed to be equdistant between a MIN and MAX val with
    *  constant interval size DX. Subsequently look-ups are FAST.
    *  __T must be built-in numeric
    **/
   // ---------------------------------------------------------------------- //
   template< typename __T = float > 
   class equidistant_key
   {
   public:
      equidistant_key ( const __T &_min, 
			const __T &_max, 
			const std::size_t _dim )
	 : mc_min(_min), mc_max(_max), 
	   mc_dx( (_max-_min)/(_dim-1) ), mc_dim(_dim)
	 { 
	 }

      std::size_t
      size ( ) const
	 { return mc_dim; }

      const __T
      operator[] ( const std::size_t _i ) const
	 {
	    return mc_min + mc_dx*_i;
	 }

      template< typename __T_ans >
      std::size_t
      operator( ) ( const __T_ans & __restrict__ _x ) const
	 {
	    const __T tmp = __T(_x) - mc_min;
	    std::size_t i = 0;
	    if ( tmp>0 )
	    {
	       i = static_cast<std::size_t>( ceil(tmp/mc_dx)  );
	       if ( i>mc_dim )  i=mc_dim;
	    }
	    return i;
	 }
     
   private:
      const __T mc_min;
      const __T mc_max;
      const __T mc_dx;
      const std::size_t mc_dim;
   };


   // ------------------------  interpolatable_range  ---------------------- //
   /** \class tabular_function
    *  PUT IN ALLOC, MANIP
    *  let it manage raw buffer
    *  lookup must pass in key data pointer
    *  client manages key and value cols independently
    *  spartan but fast
    **/
   // ---------------------------------------------------------------------- //
   template< typename __T       = float ,
	     typename __F_l     = forward_linear, 
	     typename __F_c     = centered_linear, 
	     typename __F_r     = backward_linear >
   class interpolatable_range
      : public interpolator_1D< __F_l, __F_c, __F_r >
   {
   private:
      typedef interpolator_1D<__F_l, __F_c, __F_r > interp_func;
      
   public:
      interpolatable_range ( ) : mc_data( ) { };

      interpolatable_range ( const std::size_t _dim, 
			     const __T * __restrict__ _p_data )
	 : mc_data( _p_data, _p_data + _dim )
	 { 
	 }

      void
      assign ( const std::size_t _dim, 
	       const __T * __restrict__ _p_data )
      { 
	 mc_data.clear();
	 mc_data.assign( _p_data, _p_data + _dim );
      }

      void push_back( __T v )
      {
	 mc_data.push_back( v );
      }

      std::size_t
      size ( ) const
	 { return mc_data.size(); }

      const __T &
      operator[] ( const std::size_t _i ) const
	 {
	    return mc_data[_i];
	 }

      __T &
      operator[] ( const std::size_t _i )
	 {
	    return mc_data[_i];
	 }

      template< typename __T_ans, template<typename> class __C_key >
      inline void
      operator ( ) ( const std::size_t _i, const __C_key<__T> &_X,
		     const __T_ans & _x, __T_ans &  y_ ) const
	 {
	    interp_func::apply( 
	       mc_data.size(), _X, mc_data, _x, y_, _i 
	       );
	 }

   private:
      std::vector<__T> mc_data;
   }; // class interpolatable_range

   // ---------------------------  tabular_function  ----------------------- //
   /** \class tabular_function
    **/
   // ---------------------------------------------------------------------- //
   template< typename __T       = float ,
	     typename __C_key   = sorted_key<__T>,
	     typename __F_l     = forward_linear, 
	     typename __F_c     = centered_linear, 
	     typename __F_r     = backward_linear >
   class tabular_function
   {
   private:
      typedef interpolatable_range<__T,__F_l, __F_c, __F_r > __C_val;
      
   public:
      
      tabular_function ( const std::size_t _dim,
			 const __C_key &_X,
			 const __C_val &_Y )
	 : mc_X(_X), mc_Y(_Y)
	 { /* done */ }

      template< typename __T_ans >
      std::size_t
      operator( ) ( const __T_ans & __restrict__ _x, __T_ans & __restrict__ y_ ) const
	 {
	    const std::size_t i = mc_X( _x );
	    mc_Y( i, mc_X, _x, y_);
	    return i;
	 }

   private:
      const __C_key mc_X;
      const __C_val mc_Y;

   }; // class tabular_function

} // ---------------------------------------------------------- END NAMESPACE

#endif // __FUNC_TABLE_HPP_INCLUDED_

//---------------------------------------------------------------------------//
//                           EOF func_table.hpp
//---------------------------------------------------------------------------//
