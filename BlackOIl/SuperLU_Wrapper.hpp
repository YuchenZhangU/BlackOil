// -------------------------------- *- C++ -* ------------------------------ //
// -------------------------------* TEMPLATE *------------------------------ //
// ------------------------------------------------------------------------- //
//! \file  SuperLU_Wrapper.hpp
//! \brief A wrapper class around the serial superlu solver dgssv
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
#ifndef __SUPERLU_WRAPPER_HPP_INCLUDED_
#define __SUPERLU_WRAPPER_HPP_INCLUDED_

#include "slu_ddefs.h"
#include "fastl/containers/pod_vector_unbounded.hpp"
#include "CSR_Matrix.hpp"

namespace GENSOL { // -------------------------------------------- BEGIN NAMESPACE 

   // ----------------------------  SuperLU_Wrapper  ----------------------- //
   /** \class SuperLU_Wrapper
    *  
    *  Wrapper around serial version of SuperLU solver
    **/
   // ------------------------------------------------------------------------- //
   class SuperLU_Wrapper
   {
   public:
      typedef double                                     double_type;
      typedef int                                        int_type;
      typedef CSR_Matrix< double_type, int_type >        A_type;
      typedef fastl::pod_vector_unbounded< double_type > x_type;

   public:
      //.............................  LIFECYCLE  ..........................//
      SuperLU_Wrapper( int _M, int _N, int _NNZ) : 
	 m_PERMC_SPEC(1), mM(_M), mN(_N), mNNZ(_NNZ), is_first( true )
      { 
	 set_default_options( &moptions );
	 moptions.ColPerm = NATURAL;
	 mlwork  = 1000*mM*mN;
	 mwork   = SUPERLU_MALLOC( mlwork );
	 mR      = (double *) SUPERLU_MALLOC( mM * sizeof( double ) );
	 mC      = (double *) SUPERLU_MALLOC( mN * sizeof( double ) );
	 mferr   = (double *) SUPERLU_MALLOC( sizeof( double ) );
	 mberr   = (double *) SUPERLU_MALLOC( sizeof( double ) );
	 metree  = intMalloc(mN);
	 mperm_r = intMalloc(mM);
	 mperm_c = intMalloc(mN);
	 mequed[0] = 'N';
      }

      ~SuperLU_Wrapper( )
      {
	 SUPERLU_FREE (mperm_c);
	 SUPERLU_FREE (mperm_r);
	 SUPERLU_FREE (metree);
	 SUPERLU_FREE (mberr);
	 SUPERLU_FREE (mferr);
	 SUPERLU_FREE (mC);
	 SUPERLU_FREE (mR);
	 SUPERLU_FREE ( mwork );
      } 
      
      //...........................  ASSIGNMENT  ...........................//      

      //.............................  ACCESS  .............................//

      //.............................  OPERATORS  .........................//
      static int_type offset( ) { return 0; }

      int solve( A_type &_A, x_type &_x, x_type & _b )
      {
	 mM = _A.M();
	 mN = _A.N();
	 mNNZ = _A.NNZ();

//	 if (! is_first) moptions.Fact = SamePattern;

	 StatInit(&mstat);	 

	 // create RHS
	 dCreate_Dense_Matrix(&mB, mM, 1, _b.begin(), mM, SLU_DN, SLU_D, SLU_GE);
	 dCreate_Dense_Matrix(&mX, mM, 1, _x.begin(), mM, SLU_DN, SLU_D, SLU_GE );

	 set_A(_A);

	 // solve
	 int info;
	 mem_usage_t mem_usage;
	 double rcond, rpg;
	 dgssvx(&moptions, &mA, mperm_c, mperm_r, metree, mequed, mR, mC,
		&mL, &mU, mwork, mlwork, &mB, &mX, &rpg, &rcond, mferr, mberr,
		&mem_usage, &mstat, &info);

	 StatFree(&mstat);	 

//	 Destroy_SuperNode_Matrix(&mL);
//	 Destroy_CompCol_Matrix(&mU);

	 return -info;
      }

   protected:
      void set_A(  A_type &_A )
      {
	 // Create Operator
	 dCreate_CompRow_Matrix(&mA, 
				mM, mN, mNNZ, 
				_A.value().begin(), _A.colind().begin(), _A.rowptr().begin(), 
				SLU_NR, SLU_D, SLU_GE);
      }

   private:
      const int   m_PERMC_SPEC; // default 1: minimum degree ordering on structure of A'*A
      int   mM;
      int   mN;
      int   mNNZ;
      bool  is_first;
      int        *mperm_r; // row permutations from partial pivoting
      int        *mperm_c; // column permutation vector 
      superlu_options_t moptions;
      SuperLUStat_t mstat;
      SuperMatrix mA;
      SuperMatrix mB;
      SuperMatrix mL;
      SuperMatrix mU;
      SuperMatrix mX;
      char        mequed[1];
      double     *mferr;
      double     *mberr;
      double     *mR;
      double     *mC;
      int        *metree;
      int         mlwork;
      void       *mwork;
   }; // SuperLU_Wrapper

}; // ---------------------------------------------------------- END NAMESPACE

#endif // __SuperLU_Wrapper_HPP_INCLUDED_

//---------------------------------------------------------------------------//
//                           EOF SuperLU_Wrapper.hpp
//---------------------------------------------------------------------------//

