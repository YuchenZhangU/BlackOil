#pragma once

#include <vector>
#include "adetl/systems/ADvector.hpp"
#include "PropertyCalculator.hpp"
#include "DiscreteDomain.hpp"

class DiscreteProblem
{
private:
	typedef PropertyCalculator::Phase     PhaseID;
	typedef PropertyCalculator::StatusID  StatusID;
	typedef struct
	{
		adetl::ADscalar<> Po;
		adetl::ADscalar<> Sw;
		adetl::ADscalar<> Sg;  // active for OWG
		adetl::ADscalar<> Rso; // active for OW
		StatusID   status;
	}                                            StateElement;

	typedef struct
	{
		bool        is_producer;
		std::size_t loc;
		double      WI;
		double      Pbh;
		double      QINJ[3];
		double      Rso;
	}                                            Well;
public:
	typedef std::vector< StateElement >          StateVector;
	typedef struct
	{
		double MatBal[3];
		double NormSat[3];
	}                                            ConvergenceInfo;
private:
	typedef PropertyCalculator::CellProps        CellProps;
	typedef std::vector< CellProps >             CellVector;

	typedef struct
	{
		double            T;
		adetl::ADscalar<> L[3];
		adetl::ADscalar<> Pot[3];
		adetl::ADscalar<> Rso;
	}                                           FaceProps;
	typedef std::vector< FaceProps >            FaceVector;

	typedef DiscreteDomain                      MeshType;

public:
	DiscreteProblem(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX,
		const std::vector<double> &vKY,
		const std::vector<double> &vKZ,
		const std::vector<double> &Phi_ref);

	std::size_t  max_num_eqns() const { return 3 * mMesh.size_cells(); }
	std::size_t  max_num_nnz() const { return 63 * mMesh.size_cells(); }

	void initialize_state(StateVector & state);

	void bind_to_old_state(const StateVector &_old_state);

	bool discretize(const StateVector & state, double DT);

	bool is_converged(ConvergenceInfo & nrm);

	template< typename V, typename M >
	void extract_R_J(V &residual, M &jacobian, std::size_t i);

	template< typename V >
	void extract_dR_dDT(V &dR_dDT);

	template< typename R >
	void update_state(StateVector &_state, const R & update, bool do_safeguard);

	void dump_residual(std::ofstream &out)
	{
		fastl::spy(out, mResidual);
	}
protected:
	void setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ,
		const std::vector<double> &vKX,
		const std::vector<double> &vKY,
		const std::vector<double> &vKZ);

	std::size_t eqnID(std::size_t c, PhaseID ph) const
	{
		std::size_t incr;
		switch (ph)
		{
		case PhaseID::O:
			incr = 0;
			break;
		case PhaseID::W:
			incr = 1;
			break;
		case PhaseID::G:
			incr = 2;
			break;
		}
		return 3 * c + incr;
	}

	std::size_t eqnID(std::size_t c, std::size_t ph) const
	{
		PhaseID incr;
		switch (ph)
		{
		case 0:
			incr = PhaseID::O;
			break;
		case 1:
			incr = PhaseID::W;
			break;
		case 2:
			incr = PhaseID::G;
			break;
		}
		return eqnID(c, incr);
	}

	void initialize_transmissibility(const std::vector<double> & KX,
		const std::vector<double> & KY,
		const std::vector<double> & KZ);
	double safeguard_MAC(double upd);
	void compute_accumulation();
	void compute_flow();
	void compute_wells();
	void compute_face_properties();
	void compute_cell_properties(const StateVector & state);

private:
	MeshType            mMesh;
	PropertyCalculator  mPropCalc;
	CellVector          mCells;
	FaceVector          mFaces;
	std::vector<Well>   mWells;
	std::vector<double> mPhi_ref;
	std::vector<double> mAccum_old;
	adetl::ADvector     mAccum;
	adetl::ADvector     mFlow;
	adetl::ADvector     mResidual;
	double              mDT;
	adetl::ADscalar<>   mTmpVars[3];
};

std::ostream &
operator << (std::ostream & ostr,
const DiscreteProblem::ConvergenceInfo & _out);

#include "DiscreteProblem_IMPL.hpp"
