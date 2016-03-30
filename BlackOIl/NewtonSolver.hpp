#pragma once
#include <iostream>
#include<string>

namespace GENSOL{
  
  template< typename Model, typename LNsolver >
  class NewtonSolver
  {
  public:
    typedef struct
    {
      std::size_t                    niter;       
      int                            ln_error_code;
      bool                           is_badresidual;
      bool                           is_converged;         
      bool                           is_toomanyiter;
      bool                           is_badlinearsolve;
      bool                           is_failed;
    }                                report_t;

  public:

    NewtonSolver( Model &model, std::size_t maxiter=1, int verbosity=0 ) :
      MAXITER(maxiter), J( model.max_num_eqns()   , model.max_num_nnz() ),
      r( model.max_num_eqns()   ), update( model.max_num_eqns() ), verbose( verbosity )
    {     }

    report_t solve_timestep( typename Model::StateVector &newState, 
			     const typename Model::StateVector &oldState,
			     double DT,
			     Model &model,
			     LNsolver &lnsolver )
    {

      status.niter = 0;
      status.ln_error_code = 0;

      model.bind_to_old_state( oldState );

      status.is_badresidual    = model.discretize( newState, DT );
      status.is_converged      = model.is_converged( nrm );
      status.is_toomanyiter    = ( status.niter > MAXITER );
      status.is_badlinearsolve = ( status.ln_error_code < 0);
      status.is_failed         = ( (status.is_toomanyiter || status.is_badlinearsolve) || 
				   status.is_badresidual );
      //if (verbose > 0) 
	  //std::cout << "\t" << status.niter << "\t" << nrm << std::endl;

      while ( (!status.is_failed) && (!status.is_converged) )
	{
	  model.extract_R_J( r, J, LNsolver::offset() );
	  for (std::size_t i=0; i<r.size(); ++i ) r[i] *= -1.0;
	  update.resize(r.size(), 0);
	  status.ln_error_code = lnsolver.solve( J, update, r );
	  model.update_state( newState, update, true );

	  ++status.niter;

	  status.is_badresidual    = model.discretize( newState, DT );
	  status.is_converged      = model.is_converged( nrm );
	  status.is_toomanyiter    = ( status.niter > MAXITER );
	  status.is_badlinearsolve = ( status.ln_error_code < 0);
	  status.is_failed         = ( (status.is_toomanyiter || status.is_badlinearsolve) || 
				       status.is_badresidual );
	  //std::cout << "failed" << status.is_badlinearsolve << std::endl;
	 // if (verbose > 0) 
	    //std::cout << "\t" << status.niter << "\t" << nrm << std::endl;

	}
      return status;
    }

  private:
    const std::size_t               MAXITER;
    typename LNsolver::A_type       J;
    typename LNsolver::x_type       r;
    typename LNsolver::x_type       update;
    const int                       verbose;
    typename Model::ConvergenceInfo nrm;
    report_t                        status;

  };
};
