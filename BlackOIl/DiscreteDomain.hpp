#include <cstddef>
#include "fastl/containers/pod_vector_unbounded.hpp"
#include "Point.hpp"

class DiscreteDomain
{
public:
	typedef std::size_t  Uint;

	typedef struct
	{
		Point  coord;
		Point  unit_normal; // pointing from cell1 to cell2
		double measure;
		Uint   cell1;
		Uint   cell2;
	} Face;

	typedef struct
	{
		Point  coord;
		double measure;
		std::vector<Uint> faces;
	} Cell;

public:
	template< typename B >
	DiscreteDomain(B builder) : m_gravity(), m_cells(), m_faces()
	{
		builder.build(m_gravity, m_cells, m_faces);
	}

	Uint size_cells() const { return m_cells.size(); }
	Uint size_faces() const { return m_faces.size(); }

	template< typename V >
	void   cell_faces(Uint _c, V & v_) const
	{
		const Uint N = m_cells[_c].faces.size();
		v_.resize(N);
		for (Uint i = 0; i < N; ++i) v_[i] = m_cells[_c].faces[i];
	}

	template< typename V >
	void   cell_neighbors(Uint _c, V & v_) const
	{
		const Uint Nf = m_cells[_c].faces.size();
		v_.resize(Nf);
		for (Uint f = 0; f < Nf; ++f)
		{
			const Uint c1 = m_faces[m_cells[_c].faces[f]].cell1;
			const Uint c2 = m_faces[m_cells[_c].faces[f]].cell2;
			v_[f] = c1 + c2 - _c;
		}
	}

	template< typename V >
	void   cell_extended_neighbors(Uint _c, V  & v_) const
	{
		// TODO
	}

	Point  cell_coord(Uint _c) const { return m_cells[_c].coord; }
	double cell_measure(Uint _c) const { return m_cells[_c].measure; }

	Uint   face_adj_cell_1(Uint _f) const { return m_faces[_f].cell1; }
	Uint   face_adj_cell_2(Uint _f) const { return m_faces[_f].cell2; }

	template< typename V >
	void   face_extended_adj_cells(Uint _f, V & v_) const
	{
		// TODO
	}

	Point  face_coord(Uint _f)   const { return m_faces[_f].coord; }
	double face_measure(Uint _f) const { return m_faces[_f].measure; }
	Point  gravity_vector()       const { return m_gravity; }
	Point  unit_normal(Uint _f)  const { return m_faces[_f].unit_normal; }

	Point  outward_unit_normal(Uint _f, Uint _c) const
	{
		Point nrml = m_faces[_f].unit_normal;
		if (_c != face_adj_cell_1(_f))
		{
			nrml.p[0] *= -1;
			nrml.p[1] *= -1;
			nrml.p[2] *= -1;
		}
		return nrml;
	};

private:
	Point m_gravity;
	std::vector< Cell > m_cells;
	std::vector< Face > m_faces;
};
