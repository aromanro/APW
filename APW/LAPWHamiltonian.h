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
	// there will be one for each l up to max l
	class Values
	{
	public:		
		double El = 0; // in general, the energy does not need to be the same for all ls. Big l values matter more for bigger energy

		double Wavefunction = 0;
		double EnergyDerivative = 0; // with a dot in the book
		double RadialDerivative = 0; // with a prime in the book
		double BothDerivative = 0; // derivative of both energy and position
		double Nl = 0; // the norm of the energy derivative of the wavefunction (6.49)

		std::pair<double, double> ComputeSlGammal(unsigned int l, double R, const Vector3D<double>& q1, const Vector3D<double>& q2) const
		{
			const double q1_length = q1.Length();
			const double q2_length = q2.Length();

			const double jlp1 = SpecialFunctions::Bessel::jderiv(l, q1_length * R) * q1_length;
			const double jl1 = SpecialFunctions::Bessel::j(l, q1_length * R);
			
			const double jlp2 = SpecialFunctions::Bessel::jderiv(l, q2_length * R) * q2_length;
			const double jl2 = SpecialFunctions::Bessel::j(l, q2_length * R);

			// 6.42b
			const double al1 = jlp1 * EnergyDerivative - jl1 * BothDerivative;
			const double al2 = jlp2 * EnergyDerivative - jl2 * BothDerivative;

			// 6.42d
			const double bl1 = jl1 * RadialDerivative - jlp1 * Wavefunction;
			const double bl2 = jl2 * RadialDerivative - jlp2 * Wavefunction;

			// 6.43b
			const double sl = al1 * al2 + bl1 * bl2 * Nl; 

			// 6.44b
			const double gammal = EnergyDerivative * RadialDerivative * (jlp1 * jl2 + jl1 * jlp2) -
				(RadialDerivative * BothDerivative * jl1 * jl2 + Wavefunction * EnergyDerivative * jlp1 * jlp2);

			return std::make_pair(sl, gammal);
		}
	};

	class Hamiltonian
	{
	public:
		Hamiltonian(const std::vector<Vector3D<double>>& basisVectors, double R, double a, double cellVolume)
			: m_basisVectors(basisVectors), m_R(R), prefactor(4. * M_PI * R * R / cellVolume)
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


		void Compute(const Vector3D<double>& k, const std::vector<Values>& vals)
		{
			const size_t m_lMax = vals.size();
			const size_t size = m_basisVectors.size();

			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j <= i; ++j)
				{
					const Vector3D<double> qi = m_basisVectors[i] + k;
					const Vector3D<double> qj = m_basisVectors[j] + k;
					const double qiqjscalar = qi * qj;

					const double qilength = qi.Length();
					const double qjlength = qj.Length();

					const double qiqj = qilength * qjlength;
					double cosTheta = qiqjscalar / qiqj;
					// some numerical issues prevented using std::legendre (the other custom implementation in Legendre::p worked fine), this solves it:
					if (isnan(cosTheta) || isinf(cosTheta) || cosTheta > 1. || cosTheta < -1) cosTheta = (cosTheta < 0 ? -1 : 1);

					// 6.43a and 6.44a
					S(i, j) = 0;
					H(i, j) = 0;
					
					for (unsigned int l = 0; l < m_lMax; ++l)
					{
						const double p = (2. * l + 1.) * SpecialFunctions::Legendre::p(l, cosTheta);
						double sl;
						double gammal;
						std::tie(sl, gammal) = vals[l].ComputeSlGammal(l, m_R, qi, qj);
						
						S(i, j) += p * sl;
						H(i, j) += p * (vals[l].El * sl + gammal);
					}
					
					S(i, j) = U(i, j) + prefactor * m_R * m_R * S(i, j); 
					H(i, j) = qiqj * U(i, j) + prefactor * H(i, j);

					if (i != j)
					{
						S(j, i) = S(i, j);
						H(j, i) = H(i, j);
					}
				}
			}
		}


		Eigen::VectorXd GetEnergies() const
		{
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(H, S);
			return es.eigenvalues();
		}


	protected:
		// technically the basis vectors are Ki + k, those here are only Ki
		const std::vector<Vector3D<double>>& m_basisVectors;
		const double m_R;
		const double prefactor;

		// see 6.43c for calculation details
		Eigen::MatrixXd U;

		// see 6.43a
		Eigen::MatrixXd S; 

		// see 6.44a
		Eigen::MatrixXd H;
	};

}

