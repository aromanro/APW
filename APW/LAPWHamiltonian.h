#pragma once

#include <iostream>

#include <Eigen/eigen>

#include <vector>

#include "Vector3D.h"

#include "SpecialFunctions.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace LAPW
{
	// Still have to compute them
	// there will be one for each l
	class Values
	{
	public:		
		double Wavefunction = 0;
		double EnergyDerivative = 0; // with a dot in the book
		double RadialDerivative = 0; // with a prime in the book
		double BothDerivative = 0; // derivative of both energy and position
		double Nl = 0; // the norm of the energy derivative of the wavefunction (6.49)


	};

	class Hamiltonian
	{
	public:
		Hamiltonian(const std::vector<Vector3D<double>>& basisVectors, double R, double a, double cellVolume, unsigned int lmax)
			: m_basisVectors(basisVectors), m_R(R), prefactor(4. * M_PI * R * R / cellVolume), m_lMax(lmax)
		{
			// compute U, it's the same as A in APW
			size_t size = basisVectors.size();
			U.resize(size, size);
			const double mtwopref = -prefactor; // -4 * M_PI * R^2 / cellVolume

			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j < i; ++j)
				{
					const double dist = (m_basisVectors[i] - m_basisVectors[j]).Length();

					U(j, i) = U(i, j) = mtwopref * SpecialFunctions::Bessel::j(1, dist * R) / dist;
				}

				U(i, i) = mtwopref * R / 3. + 1; // 1 is from delta
			}


			S.resize(size, size);
			H.resize(size, size);
		}


	protected:
		// technically the basis vectors are Ki + k, those here are only Ki
		const std::vector<Vector3D<double>>& m_basisVectors;
		const double m_R;
		const double prefactor;
		const unsigned int m_lMax;

		// see 6.43c for calculation details
		Eigen::MatrixXd U;

		// see 6.43a
		Eigen::MatrixXd S; 

		// see 6.44a
		Eigen::MatrixXd H;
	};

}

