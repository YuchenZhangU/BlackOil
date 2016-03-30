// -------------------------------- *- C++ -* ------------------------------ //
// -------------------------------* TEMPLATE *------------------------------ //
// ------------------------------------------------------------------------- //
//! \file  interpolators.hpp
//! \brief Functor templates for interpolation and extrapolation.
//
//! Several functors are defined to perform 1D interpolation/extrapolation.
//! All functors have a common apply(...) function interface. They differ
//! in the type of approximation used to compute values.
//!
//! TEMPLATE PARAMETERS::
//!
//!   __T_tbl_col 
//!       Data-type of a table column. Cannot be a pointer or reference. Must
//!       a model of some sort of Indexable concept, i.e.
//!              +  is a model of a collection container concept with
//!                 element type __T_tbl_col::element_type::value_type
//!              +  elem type is a field that is orderable & comparable
//!              +  elem type is closed under arithmetic ops
//!              +  collcetion is indexable, i.e. has 
//!                      elem_type  operator[]( std::size_t ) and
//!                      size_t     size      ( void)
//!
//!   __T_ans
//!       Querry key and return value data-type. Can be any object/built-in
//!       type that supports heterogeneous arithmetic and assignment with 
//!       typename __T_tbl_col::element_type::value_type
//!
//! COMMON APPLYER PROTOTYPE::
//!   void apply( std::size_t _i,
//!               const __T_ans &_x,     __T_ans &y_,
//!               const __T_tbl_col &_X, __T_tbl_col &_Y );
//!
//!   @param _i input   Index such that _X[_i-1] <= _x <= _X[_i]
//!   @param _x input   Target key, at which to evaluate table estimate
//!   @param y_ output  Resulting look-up value stored in this
//!   @param _X input   SORTED collection of table keys in INCREASING order
//!   @param _Y input   Corresponding collection of table values 
//!
//! USAGE NOTES::
//!   Not all functors support both interpolation and extrapolation.
//!
//! DEV NOTES::
//!   A few todos... look into second order approximations. Make container
//!   concepts "concept-checkable".
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
#ifndef __INTERPOLATORS_HPP_INCLUDED_
#define __INTERPOLATORS_HPP_INCLUDED_

#include <cstddef> // provides size_t

namespace fastl { // ----------------------------------------- BEGIN NAMESPACE 

   // --------------------------  interpolator_1D  ------------------------- //
   /** \class interpolator_1D
    *  
    *  This is a dispatch functor whose sole purpose is to combine 3 
    *  interpolation functors. The first is used to extrapolate data on 
    *  the lower bound of a sorted collection. The second interpolates in the 
    *  interior, The third extrapolates on the upper boundary.
    **/
   // ---------------------------------------------------------------------- //
   template< typename __F_Extrap_l,
	     typename __F_Interp,
	     typename __F_Extrap_u >
   struct interpolator_1D
      : public __F_Extrap_l, public __F_Interp, public __F_Extrap_u
   { 

      template< typename __X, typename __Y, typename __T_ans  >
      static inline void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
	 {
	    // Dev notes:
	    // Jump function with a little instruction scheduling in case
	    //..todo.. Watch for register spillage on AMD's
	    const unsigned short jmp_1 = static_cast<unsigned short>( _i==_n );
	    unsigned short jmp_2 = static_cast<unsigned short>( _i==0 );
	    jmp_2 <<= 1;
	    const unsigned short jmp = jmp_1 + jmp_2;
	    
	    switch ( jmp )
	    {
	    case 0:
	       __F_Interp::apply ( _n, _X, _Y, _x, y_, _i );
	       break;

	    case 1:	       
	       __F_Extrap_u::apply ( _n, _X, _Y, _x, y_, _i );
	       break;

	    default:
	       __F_Extrap_l::apply ( _n, _X, _Y, _x, y_, _i );
	       break;
	    } // compiler makes a jmp table with better locality (hopefully)

	 } // apply

   protected:
      ~interpolator_1D ( ) { /* no polymorphic deletion */ }

   }; // struct interpolator_1D


   // ---------------------------  forward_constant  ----------------------- //
   /** \class forward_constant
    *  
    *  Piecewise constant function with the value from the right, i.e.
    *  Given, i : X_{i-1} <= x <= X_i , get,
    *         y(x) = Y_i , 
    * 
    *  PARAMS:  See file header comment for important assumptions on inputs.
    *  WARNING: Not a valid upper extrapolator.
    *  ASSERTIONS: Requires i < length(X)
    **/
   // ---------------------------------------------------------------------- //
   struct forward_constant
   {
      template< typename __X, typename __Y, typename __T_ans  >
      static inline void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
	 {
	    y_ = __T_ans( _Y [ _i ] );
	 }

   protected:
      ~forward_constant ( ) { /* no polymorphic deletion */ }
   };

   // -----------------------------  forward_linear  ----------------------- //
   /** \class forward_linear
    *  
    *  Piecewise linear function with slope taken from the right, i.e.
    *  Given,  i : X_{i-1} <= x <= X_i , get,
    *         \sigma = (Y_{i+1} - Y_{i})/(X_{i+1} - X_{i})
    *          y(x)  = Y_{i} + \sigma (x - X_{i})
    * 
    *  PARAMS:  See file header comment for important assumptions on inputs.
    *  WARNING: Not a valid upper extrapolator. Not valid interpolator.
    *  ASSERTIONS: Requires i < length(X)-1
    **/
   // ---------------------------------------------------------------------- //
   struct forward_linear
   {
      template< typename __X, typename __Y, typename __T_ans  >
      static void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
      {
	 const std::size_t i_c = _i;
	 const std::size_t i_r = _i + 1;
	 
	 y_ = ::operator+( _Y[i_c] ,
			   ::operator*( ( _Y[i_r] - _Y[i_c] ) / ( _X[i_r] - _X[i_c] ) , 
					::operator-( _x , _X[i_c]) ) );
      }
      
   protected:
      ~forward_linear ( ) { /* no polymorphic deletion */ }
   };


   // ----------------------------  centered_linear  ----------------------- //
   /** \class centered_linear
    *  
    *  Piecewise linear function with centered slope, i.e.
    *  Given , i : X_{i-1} <= x <= X_i , get , 
    *          \sigma = (Y_{i} - Y_{i-1})/(X_{i} - X_{i-1})
    *           y(x)  = Y_{i-1} + \sigma (x - X_{i-1}) 
    * 
    *  PARAMS:  See file header comment for important assumptions on inputs.
    *  WARNING: Not a valid extrapolator (neither upper nor lower).
    *  ASSERTIONS: Requires 0 < i < length(X)
    **/
   // ---------------------------------------------------------------------- //
   struct centered_linear
   {
      template< typename __X, typename __Y, typename __T_ans  >
      static inline void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
	 {
	    const std::size_t i_r = _i;
	    const std::size_t i_l = _i - 1;
	    y_ = ::operator+( _Y[i_l] , 
			      ::operator*( (( _Y[i_r] - _Y[i_l] ) / ( _X[i_r] - _X[i_l] )),
					   ::operator-( _x, _X[i_l] ) ) );
	 }

   protected:
      ~centered_linear ( ) { /* no polymorphic deletion */ }
   };

   // ----------------------------  backward_linear  ----------------------- //
   /** \class backward_linear
    *  
    *  Piecewise linear function with slope taken from the left, i.e.
    *  Given,  i : X_{i-1} <= x <= X_i , get,
    *         \sigma = (Y_{i-1} - Y_{i-2})/(X_{i-1} - X_{i-2})
    *          y(x)  = Y_{i-1} + \sigma (x - X_{i-1})
    * 
    *  PARAMS:  See file header comment for important assumptions on inputs.
    *  WARNING: Not a valid lower extrapolator. Not valid interpolator
    *  ASSERTIONS: Requires 1 < i <= length(X)
    **/
   // ---------------------------------------------------------------------- //
   struct backward_linear
   {

      template< typename __X, typename __Y, typename __T_ans  >
      static inline void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
	 {
	    const std::size_t i_l  = _i - 1;
	    const std::size_t i_ll = i_l - 1;

	    y_ = ::operator+( _Y[i_l] ,
			      ::operator*( ( _Y[i_l] - _Y[i_ll] ) / ( _X[i_l] - _X[i_ll] ) ,
					   ::operator-( _x , _X[i_l]) ) );
	 }

   protected:
      ~backward_linear ( ) { /* no polymorphic deletion */ }
   };

   // ---------------------------  backward_constant  ---------------------- //
   /** \class backward_constant
    *  
    *  Piecewise constant function with the value from the left, i.e.
    *  Given, i : X_{i-1} <= x <= X_i , get,
    *         y(x) = Y_{i-1} , 
    * 
    *  PARAMS:  See file header comment for important assumptions on inputs.
    *  WARNING: Not a valid lower extrapolator.
    *  ASSERTIONS: Requires 0 < i <= length(X)
    **/
   // ---------------------------------------------------------------------- //
   struct backward_constant
   {
      template< typename __X, typename __Y, typename __T_ans  >
      static inline void
      apply( const std::size_t _n, const __X &_X, const __Y &_Y,
	     const __T_ans &_x, __T_ans &y_, const std::size_t _i )
	 {
	    y_ = __T_ans( _Y [ _i-1 ] );
	 }

   protected:
      ~backward_constant ( ) { /* no polymorphic deletion */ }
   };

} // ---------------------------------------------------------- END NAMESPACE

#endif // __INTERPOLATORS_HPP_INCLUDED_

//---------------------------------------------------------------------------//
//                           EOF interpolators.hpp
//---------------------------------------------------------------------------//
