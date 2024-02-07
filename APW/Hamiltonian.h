#pragma once

#include <iostream>

#include <Eigen/eigen>

#include <vector>

#include "Vector3D.h"

#include "SpecialFunctions.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace APW
{

	class Secular
	{
	public:
		Secular(const std::vector<Vector3D<double>>& basisVectors, double R, double cellVolume, unsigned int lmax)
		: m_basisVectors(basisVectors), m_R(R), prefactor(2. * M_PI * R / cellVolume), m_lMax(lmax)
		{		
			// compute A
			const size_t size = basisVectors.size();
			OverlapInterstitial.resize(size, size);
			const double mtwopref = -2. * prefactor * R; // -4 * M_PI * R^2 / cellVolume

			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j < i; ++j)
				{
					const double dist = (m_basisVectors[i] - m_basisVectors[j]).Length();

					OverlapInterstitial(j, i) = OverlapInterstitial(i, j) = mtwopref * SpecialFunctions::Bessel::j(1, dist * R) / dist;
				}
				
				OverlapInterstitial(i, i) = mtwopref * R / 3. + 1; // 1 is from delta
			}

			// used for computing B
			OverlapInterstitialHalf = 0.5 * OverlapInterstitial; // 0.5 is from hbar / (2 * me), in atomic units becomes 1/2

			B.resize(size, size);
			C.resize(lmax + 1ULL); 
			for (size_t l = 0; l <= lmax; ++l)
				C[l].resize(size, size);
		}

		void ComputeBC(const Vector3D<double>& k)
		{
			const size_t size = m_basisVectors.size();

			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j <= i; ++j)
				{
					// compute B
					const Vector3D qi(m_basisVectors[i] + k);
					const Vector3D qj(m_basisVectors[j] + k);
					const double qiqjscalar = qi * qj;

					B(i, j) = B(j, i) = OverlapInterstitialHalf(i, j) * qiqjscalar;

					// compute C
					const double qilength = qi.Length();
					const double qjlength = qj.Length();

					const double qiqj = qilength * qjlength;
					double cosTheta = qiqjscalar / qiqj;
					// some numerical issues prevented using std::legendre (the other custom implementation in Legendre::p worked fine), this solves it:
					if (isnan(cosTheta) || isinf(cosTheta) || cosTheta > 1. || cosTheta < -1) cosTheta = (cosTheta < 0 ? -1 : 1);

					const double qilengthm_R = qilength * m_R;
					const double qjlengthm_R = qjlength * m_R;

					for (unsigned int l = 0; l <= m_lMax; ++l)
						C[l](i, j) = C[l](j, i) = prefactor * (2. * l + 1.) * SpecialFunctions::Legendre::p(l, cosTheta) * SpecialFunctions::Bessel::j(l, qilengthm_R) * SpecialFunctions::Bessel::j(l, qjlengthm_R);
				}
			}
		}

		void ComputeHamiltonian(double E, const std::vector<double>& ratios)
		{
			assert(ratios.size() == m_lMax + 1);
			assert(C.size() == m_lMax + 1);

			H = B - E * OverlapInterstitial;

			// the -1 below comes from u = r R
			// in 'ratios' there is u'/u
			// but R'/R is needed. If you substitute R = u/r the -1 comes out nicely.
			
			// also the term with the logarithmic derivative of the Bessel function 
			// (the term for the boundary on the exterior of the sphere) is missing because the sum turns out to be zero
			for (unsigned int l = 0; l <= m_lMax; ++l)
				H += C[l] * (m_R * ratios[l] - 1);
		}

		double Determinant() const
		{
			return H.determinant();
			//return H.fullPivLu().determinant();
		}

	protected:
		// technically the basis vectors are Ki + k, those here are only Ki
		const std::vector<Vector3D<double>>& m_basisVectors;
		const double m_R;
		const double prefactor;
		const unsigned int m_lMax;

		// see 6.35 for calculation details
		Eigen::MatrixXd OverlapInterstitial; // noted with A in the book
		Eigen::MatrixXd OverlapInterstitialHalf; // 0.5 * A, will be used to compute B

		// Call ComputeBC to fill those, needs k, so each time k is changing, call it again
		Eigen::MatrixXd B; // 0.5 * Aij * qi * qj
		std::vector<Eigen::MatrixXd> C; 

		// needs both k and E to be computed, so call ComputeHamiltonian each time those change
		Eigen::MatrixXd H;
	};

}

