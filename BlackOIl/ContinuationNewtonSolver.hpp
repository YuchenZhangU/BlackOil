#pragma once
#include <iostream>
#include "NewtonSolver.hpp"

template< typename Model, typename LNsolver >
class ContinuationNewtonSolver
{
public:
   typedef struct
   {
      bool                           is_converged;
      std::size_t                    niter;
      std::size_t                    npcstages;
      std::size_t                    ncorrector;
      std::size_t                    npredictor;
      double                         dt;
   } report_t;

public:

   ContinuationNewtonSolver( Model &model, std::size_t maxcorr=10, int verbosity=0 ) :
      NTARGCORR( 9 ),
      J( model.max_num_eqns(), model.max_num_eqns(), model.max_num_nnz() ),
      drdt( model.max_num_eqns() ), 
      tangent( model.max_num_eqns() ), 
      tangent_trial( model.max_num_eqns() ), 
      verbose( verbosity ),
      corrector( model, maxcorr, verbosity )
   { 
      MAXERR[0] = 25.0;
      MAXERR[1] = 0.2;
      MAXERR[2] = 0.2;
      drdt.resize( model.max_num_eqns() );
      tangent.resize( model.max_num_eqns() );
      tangent_trial.resize( model.max_num_eqns() );
   }

   report_t solve_timestep( typename Model::StateVector &newState, 
			    const typename Model::StateVector &oldState,
			    double DT_TARG,
			    Model &model,
			    LNsolver &lnsolver )
   {
//      MAXERR[0] = 25.0;
//      MAXERR[1] = 0.2;
//      MAXERR[2] = 0.2;

      if (verbose > 0)
	 std::cout << "Predctor-Corrector Iterations" << std::endl;
	 
      double alpha = 0.0;
      double dt    = 0.0;
      std::size_t ncorr = 0;
      newState     = oldState;
      model.bind_to_old_state( oldState );
      initialize_report( );
      while( dt < DT_TARG )
      {
	 compute_tangent( newState, dt, tangent, model, lnsolver );
	 compute_tangent_norms( newState );
	 alpha = select_steplength( dt, DT_TARG );
	 status.niter      += 1;
	 status.npcstages  += 1;
	 status.npredictor += 1;

	 compute_tangent( newState, dt+alpha, tangent_trial, model, lnsolver );
	 for (std::size_t i=0; i<tangent.size(); ++i)
	    tangent[i] = 0.5 * ( tangent[i] + tangent_trial[i] );
	 compute_tangent_norms( newState );
	 alpha = select_steplength( dt, DT_TARG );
	 status.niter      += 1;
	 status.npredictor += 1;

	 int  ncuts = 0;
	 bool keep_trying = true;
	 while (keep_trying)
	 {
	    double dt_trial = dt + alpha;
	    newState_trial  = newState;

	    // I/O
	    if (verbose > 0)
	    {
	       std::cout << dt << " ----> " << dt_trial << std::endl;
	       std::cout << "\t ALPHA = " << alpha << std::endl;
	       std::cout << "\t TANGENT NORMS: " 
			 << tnorm[0] << "\t" 
			 << tnorm[1] << "\t" 
			 << tnorm[2] << std::endl;
	       std::cout << "\t PREDICTED ERRORS: " 
			 << MAXERR[0] << "\t" 
			 << MAXERR[1] << "\t" 
			 << MAXERR[2] << std::endl;
	    }
	    
	    // update state
	    for (std::size_t i=0; i<tangent_trial.size(); ++i ) tangent_trial[i] = tangent[i]*alpha;
	    model.update_state( newState_trial, tangent_trial, false );
	    
	    // corrector
	    corrstatus = corrector.solve_timestep( newState_trial, oldState, dt_trial, model, lnsolver );
	    
	    // update report
	    status.niter      += corrstatus.niter;
	    status.ncorrector += corrstatus.niter;

	    if (corrstatus.is_converged)
	    {
	       status.dt = dt_trial;
	       keep_trying = false;
	       newState    = newState_trial;
	       dt          = dt_trial;
	       ncorr       = corrstatus.niter;
	       adjust_constraints( ncorr );
	       if ( verbose > 0) std::cout << "PASSED" << std::endl;
	    }
	    else
	    { 
	       ++ncuts;
	       if ( ncuts > 4 ) // terminate
	       {
		  keep_trying = false;
		  dt = DT_TARG;
	       }
	       tighten_constraints( );
	       alpha = select_steplength( dt, DT_TARG );
	       if ( verbose > 0) std::cout << "FAILED... cutting alpha" << std::endl;
	    }
	 }
	 status.is_converged  = corrstatus.is_converged;
      }
      
      if (verbose > 0)
	 std::cout << "==========================================" << std::endl;
      return status;
   }

protected:
   void compute_tangent( typename Model::StateVector &newState,
			 double dt,
			 typename LNsolver::x_type &tang,
			 Model &model,
			 LNsolver &lnsolver  )
   {
      model.discretize( newState, dt ); // assume this always works
      model.extract_R_J( drdt, J, LNsolver::offset() );
      model.extract_dR_dDT( drdt );
      for (std::size_t i=0; i<drdt.size(); ++i ) drdt[i] *= -1.0;
      lnsolver.solve( J, tang, drdt ); // assume this always works
   }

   void compute_tangent_norms( typename Model::StateVector &newState )
   {
      tnorm[0] = 0.0, tnorm[1]=0.0, tnorm[2]=0.0;
      const std::size_t N = tangent.size()/3;
      for (std::size_t i=0; i<N; ++i )
      {
	 tnorm[0] += std::pow( tangent[3*i], 2);
	 tnorm[1] = std::max( tnorm[1], std::abs(tangent[3*i+1]) );
	 tnorm[2] = std::max( tnorm[2], std::abs(tangent[3*i+2]) );
      }
      tnorm[0]  = std::sqrt(tnorm[0]);
   }

   double select_steplength( double dt, double DT_TARG )
   {
      is_constrained[0] = is_constrained[1] = is_constrained[2] = false;

      double alpha = MAXERR[0]/tnorm[0];
      is_constrained[0] = true;

      if ( MAXERR[1]/tnorm[1] < alpha )
      {
	 alpha = MAXERR[1]/tnorm[1];
	 is_constrained[1] = true;
	 is_constrained[0] = false;
      }

      if ( MAXERR[2]/tnorm[2] < alpha )
      {
	 alpha = MAXERR[2]/tnorm[2];
	 is_constrained[2] = true;
	 is_constrained[1] = false;
	 is_constrained[0] = false;
      }

      if ( (DT_TARG - dt) < alpha )
      {
	 alpha = DT_TARG - dt;
	 is_constrained[2] = true;
	 is_constrained[1] = true;
	 is_constrained[0] = true;
	 MAXERR[0] = alpha * tnorm[0];
	 MAXERR[1] = alpha * tnorm[1];
	 MAXERR[2] = alpha * tnorm[2];
      }
	 
      return alpha;
   }

   void adjust_constraints( std::size_t niter )
   {
      if (niter==0) ++niter;

      if (is_constrained[0]) MAXERR[0] *= NTARGCORR/static_cast<double>(niter); 
      if (is_constrained[1]) MAXERR[1] *= NTARGCORR/static_cast<double>(niter); 
      if (is_constrained[2]) MAXERR[2] *= NTARGCORR/static_cast<double>(niter); 
   }

   void tighten_constraints( )
   {
      const double CUTFACTOR = 0.5;
      if (is_constrained[0]) MAXERR[0] *= CUTFACTOR;
      if (is_constrained[1]) MAXERR[1] *= CUTFACTOR; 
      if (is_constrained[2]) MAXERR[2] *= CUTFACTOR; 
   }

   void initialize_report( )
   {
      status.is_converged = true;
      status.niter        = 0;
      status.npcstages    = 0;
      status.ncorrector   = 0;
      status.npredictor   = 0;
      status.dt           = 0.0;
   }

private:
   const std::size_t                               NTARGCORR;
   typename LNsolver::A_type                       J;
   typename LNsolver::x_type                       drdt;
   typename LNsolver::x_type                       tangent;
   typename LNsolver::x_type                       tangent_trial;
   const int                                       verbose;
   NewtonSolver< Model, LNsolver >                 corrector;
   double                                          MAXERR[3];
   report_t                                        status;
   typename NewtonSolver<Model,LNsolver>::report_t corrstatus;
   double                                          tnorm[3];
   typename Model::StateVector                     newState_trial;
   bool                                            is_constrained[3];
};
