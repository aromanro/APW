#pragma once

#include <vector>
#include <atomic>

#include "Vector3D.h"

#include "SymmetryPoints.h"
#include "Options.h"

namespace APW
{

	class BandStructure
	{
	public:
		// pass rmax as you want it or negative to have it computed as for touching spheres
		BandStructure(double a = 6.8219117, double rmax = /*2.41191*/-1.);

		SymmetryPoints symmetryPoints;

		std::vector<std::vector<double>> results;
		std::vector<unsigned int> symmetryPointsPositions;

		void Initialize(std::vector<std::string> path, unsigned int nrPoints = 600,  unsigned int nearestNeighborsNumber = 10);
		std::vector<std::vector<double>> Compute(const std::atomic_bool& terminate, const Options& options);

		unsigned int GetPointsNumber() const { return static_cast<unsigned int>(kpoints.size()); }

		const std::vector<std::string>& GetPath() const { return m_path; }
	private:
		std::vector<std::string> m_path;

		std::vector<Vector3D<double>> basisVectors;

		std::vector<Vector3D<double>> kpoints;

		unsigned int nearestNeighbors;
		std::vector<unsigned int> G2;

		double m_a;
		double m_Rmax;

		bool GenerateBasisVectors(unsigned int nearestNeighborsNumber);
	};

}

