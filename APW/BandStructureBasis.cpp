#include "BandStructureBasis.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace APW
{

	BandStructureBasis::BandStructureBasis(double a, double rmax)
		: m_a(a), m_Rmax(rmax)
	{
		// if the passed value was zero or negative, make them touching spheres
		if (m_Rmax <= 0)
			m_Rmax = sqrt(2.) * m_a / 4.;

		basisVectors.reserve(137);
	}


	bool BandStructureBasis::GenerateBasisVectors(unsigned int nearestNeighborsNumber)
	{
		if (nearestNeighborsNumber < 2 || nearestNeighborsNumber > 10) return false;

		static const std::vector<unsigned int> G2{ 0, 3, 4, 8, 11, 12, 16, 19, 20, 24 };
		const unsigned int nearestNeighbors = nearestNeighborsNumber - 1;
		basisVectors.clear();

		const int size = static_cast<int>(ceil(sqrt(static_cast<double>(G2[nearestNeighbors]))));

		// the Bravais lattice is a fcc lattice
		// the reciprocal lattice is a bcc lattice 

		// the basis vectors for the reciprocal space (without the 2 pi / a)

		const Vector3D<double> b1(-1, 1, 1), b2(1, -1, 1), b3(1, 1, -1);

		// you can also get them from the Bravais lattice vectors
		// like this:

		// Bravais lattice vectors for fcc (they should be multiplied by the lattice constant):
		//const Vector3D<double> a1(0, 0.5, 0.5), a2(0.5, 0, 0.5), a3(0.5, 0.5, 0);

		// the volume of the cell is a1 * (a2 % a3) which gives 1/4 (multiplied with a^3, of course)

		// reciprocal lattice (they should be multiplied by 2 pi / lattice constant):
		//const Vector3D<double> b1 = a2 % a3 / (a1 * (a2 % a3)); 		
		//const Vector3D<double> b2 = a3 % a1 / (a2 * (a3 % a1));
		//const Vector3D<double> b3 = a1 % a2 / (a3 * (a1 % a2));
		// the denominator is the volume, mentioned above

		for (int i = -size; i <= size; ++i)
			for (int j = -size; j <= size; ++j)
				for (int k = -size; k <= size; ++k)
				{
					const Vector3D<double> vect = b1 * i + b2 * j + b3 * k; // reciprocal lattice vector

					const double vectSquared = vect * vect;

					if (vectSquared <= G2[nearestNeighbors]) // if it's under the cutoff length, add it
						basisVectors.push_back(vect);
				}

		return true;
	}


	void BandStructureBasis::Initialize(std::vector<std::string> path, unsigned int nrPoints, unsigned int nearestNeighborsNumber)
	{
		kpoints.clear();
		kpoints.reserve(nrPoints);

		m_path.swap(path);

		GenerateBasisVectors(nearestNeighborsNumber);

		const double recVectPre = 2. * M_PI / m_a;

		// adjust 'basis' vectors
		size_t size = basisVectors.size();
		for (unsigned int i = 0; i < size; ++i)
			basisVectors[i] *= recVectPre;

		kpoints = symmetryPoints.GeneratePoints(m_path, nrPoints, symmetryPointsPositions);

		// adjust kpoints
		for (unsigned int i = 0; i < kpoints.size(); ++i)
			kpoints[i] *= recVectPre;
	}

}
