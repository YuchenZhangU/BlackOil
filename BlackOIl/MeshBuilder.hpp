#pragma once
//---------------- Build Mesh and Connection list ------------------
// provide method called 'build' to updata elements cell and face
#include "Point.hpp"
#include <cstddef>

template< typename topology >
class MeshBuilder
{
};

class UniformCartesian;
template< >
class MeshBuilder< UniformCartesian >
{
public:
	MeshBuilder(std::size_t NX, std::size_t NY, std::size_t NZ,
		double LX, double LY, double LZ)
	{
		N[0] = NX;
		N[1] = NY;
		N[2] = NZ;
		L[0] = LX;
		L[1] = LY;
		L[2] = LZ;
		DL[0] = LX / static_cast<double>(NX);
		DL[1] = LY / static_cast<double>(NY);
		DL[2] = LZ / static_cast<double>(NZ);
	}

	template< typename T1, typename T2 >
	void build(Point &gravity, T1 & cells, T2 &faces)
	{
		gravity.p[0] = 0;
		gravity.p[1] = 0;
		gravity.p[2] = 1;

		cells.resize(N[0] * N[1] * N[2]);
		for (std::size_t c = 0; c < cells.size(); ++c)
		{
			std::size_t i, j, k;
			loc_to_coord(c, i, j, k);

			cells[c].coord.p[0] = DL[0] * (0.5 + static_cast<double>(i));
			cells[c].coord.p[1] = DL[1] * (0.5 + static_cast<double>(j));
			cells[c].coord.p[2] = DL[2] * (0.5 + static_cast<double>(k));

			cells[c].measure = DL[0] * DL[1] * DL[2];

			// collection of face ids // BUG !!
			cells[c].faces.clear();
			std::size_t offset = (N[0] - 1)*N[1] * k + (N[0] - 1)*j + i;
			if (i < (N[0] - 1)) cells[c].faces.push_back(offset);
			if (i > 0) cells[c].faces.push_back(offset - 1);

			offset = (N[0] - 1)*N[1] * N[2] + N[0] * (N[1] - 1)*k + N[0] * j + i;
			if (j < (N[1] - 1)) cells[c].faces.push_back(offset);
			if (j > 0) cells[c].faces.push_back(offset - N[0]);

			offset = (N[0] - 1)*N[1] * N[2] + N[0] * (N[1] - 1)*N[2];
			if (k < (N[2] - 1)) cells[c].faces.push_back(offset + c);
			if (k > 0) cells[c].faces.push_back(offset + c - N[0] * N[1]);
		}

		faces.resize((N[0] - 1)*N[1] * N[2] + N[0] * (N[1] - 1)*N[2] + N[0] * N[1] * (N[2] - 1));
		std::size_t f = 0;
		for (std::size_t k = 0; k < N[2]; ++k)
			for (std::size_t j = 0; j < N[1]; ++j)
				for (std::size_t i = 0; i < (N[0] - 1); ++i)
				{
					faces[f].coord.p[0] = DL[0] * (1 + i);
					faces[f].coord.p[1] = DL[1] * (0.5 + j);
					faces[f].coord.p[2] = DL[2] * (0.5 + k);
					faces[f].unit_normal.p[0] = 1;
					faces[f].unit_normal.p[1] = 0;
					faces[f].unit_normal.p[2] = 0;
					faces[f].measure = DL[1] * DL[2];
					faces[f].cell1 = coord_to_loc(i, j, k);
					faces[f].cell2 = coord_to_loc(i + 1, j, k);
					++f;
				}
		for (std::size_t k = 0; k < N[2]; ++k)
			for (std::size_t j = 0; j < (N[1] - 1); ++j)
				for (std::size_t i = 0; i < N[0]; ++i)
				{
					faces[f].coord.p[0] = DL[0] * (0.5 + i);
					faces[f].coord.p[1] = DL[1] * (1 + j);
					faces[f].coord.p[2] = DL[2] * (0.5 + k);
					faces[f].unit_normal.p[0] = 0;
					faces[f].unit_normal.p[1] = 1;
					faces[f].unit_normal.p[2] = 0;
					faces[f].measure = DL[0] * DL[2];
					faces[f].cell1 = coord_to_loc(i, j, k);
					faces[f].cell2 = coord_to_loc(i, j + 1, k);
					++f;
				}
		for (std::size_t k = 0; k < (N[2] - 1); ++k)
			for (std::size_t j = 0; j < N[1]; ++j)
				for (std::size_t i = 0; i < N[0]; ++i)
				{
					faces[f].coord.p[0] = DL[0] * (0.5 + i);
					faces[f].coord.p[1] = DL[1] * (0.5 + j);
					faces[f].coord.p[2] = DL[2] * (1.0 + k);
					faces[f].unit_normal.p[0] = 0;
					faces[f].unit_normal.p[1] = 0;
					faces[f].unit_normal.p[2] = 1;
					faces[f].measure = DL[0] * DL[1];
					faces[f].cell1 = coord_to_loc(i, j, k);
					faces[f].cell2 = coord_to_loc(i, j, k + 1);
					++f;
				}
	}

protected:
	void loc_to_coord(std::size_t c, std::size_t &i, std::size_t &j, std::size_t &k)
	{
		const std::size_t Nplane = N[0] * N[1];
		k = c / Nplane;
		j = (c - k*Nplane) / N[0];
		i = (c - k*Nplane) % N[0];
	}

	std::size_t coord_to_loc(std::size_t i, std::size_t j, std::size_t k)
	{
		return i + j*N[0] + k*(N[0] * N[1]);
	}

private:
	std::size_t N[3];
	double      L[3];
	double      DL[3];
};
