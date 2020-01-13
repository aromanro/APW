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
			const double q1_lengthR = q1_length * R;
			const double q2_lengthR = q2_length * R;

			const double jl1 = SpecialFunctions::Bessel::j(l, q1_lengthR);
			const double jlp1 = SpecialFunctions::Bessel::jderiv(l, q1_lengthR) * q1_length;
			
			const double jl2 = SpecialFunctions::Bessel::j(l, q2_lengthR);
			const double jlp2 = SpecialFunctions::Bessel::jderiv(l, q2_lengthR) * q2_length;

			// I derived al and bl, the formulas are correct!

			// 6.42b
			const double al1 = jlp1 * EnergyDerivative - jl1 * BothDerivative;
			const double al2 = jlp2 * EnergyDerivative - jl2 * BothDerivative;

			// 6.42d
			const double bl1 = jl1 * RadialDerivative - jlp1 * Wavefunction;
			const double bl2 = jl2 * RadialDerivative - jlp2 * Wavefunction;

			// this is for overlap in the muffin, also correct, derived it, too

			// 6.43b
			const double sl = al1 * al2 + bl1 * bl2 * Nl; 

			// 6.44b
			const double gammal = EnergyDerivative * RadialDerivative * (jlp1 * jl2 + jl1 * jlp2)
				- (RadialDerivative * BothDerivative * jl1 * jl2 + Wavefunction * EnergyDerivative * jlp1 * jlp2);
			
			// this is an alternative, if this is used, use the commented 0.5 * (qi2 + qj2) in the Hamiltonian for the intersitial, instead of qi * qj!
			//const double gammal = 0.5 * (al1 * bl2 + al2 * bl1);

			return std::make_pair(sl, gammal);
		}
	};

	class Hamiltonian
	{
	public:
		Hamiltonian(const std::vector<Vector3D<double>>& basisVectors, double R, double cellVolume)
			: m_basisVectors(basisVectors), m_R(R), prefactor(4. * M_PI * R * R / cellVolume), prefactor2(prefactor * R * R)
		{
			// compute U, it's the same as A in APW
			// it's actually the overlap for interstitial, it works for APW so it's good for LAPW, too
			// it's the integral over the whole space for plane wave, minus the integral for the muffin
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

				U(i, i) = mtwopref * R / 3. + 1.; // 1 is from delta
			}

			S.resize(size, size);
			H.resize(size, size);
		}

		void Compute(const Vector3D<double>& k, const std::vector<Values>& vals)
		{
			const size_t m_lMax = vals.size();
			const size_t size = m_basisVectors.size();

			// I derived everything, it's correct

			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j <= i; ++j)
				{
					const Vector3D<double> qi = m_basisVectors[i] + k;
					const Vector3D<double> qj = m_basisVectors[j] + k;
					const double qiqjscalar = qi * qj;
					const double qi2 = qi * qi;
					const double qj2 = qj * qj;

					const double qilength = sqrt(qi2);
					const double qjlength = sqrt(qj2);

					const double qiqj = qilength * qjlength;
					double cosTheta = (qiqj == 0) ? 1 : qiqjscalar / qiqj;
					// some numerical issues prevented using std::legendre (the other custom implementation in Legendre::p worked fine), this solves it:
					if (isnan(cosTheta) || isinf(cosTheta) || cosTheta > 1. || cosTheta < -1) cosTheta = (cosTheta < 0 ? -1 : 1);

					// 6.43a and 6.44a
					double s = 0;
					double h = 0;
					
					for (unsigned int l = 0; l < m_lMax; ++l)
					{
						const double twolp1 = 2. * l + 1.;
						const double p = SpecialFunctions::Legendre::p(l, cosTheta);
						double sl;
						double gammal;
						std::tie(sl, gammal) = vals[l].ComputeSlGammal(l, m_R, qi, qj);

						s += twolp1 * p * sl;

						// the energy is given in Hartrees, whence the 2.
						// gammal is the matrix element for H - E, that's why the E * overlap (in the muffin) is added
						const double v = 2. * vals[l].El * sl + gammal;

						h += twolp1 * p * v;
					}
					
					// overlap for interstitial + overlap for muffin
					// this is good
					S(i, j) = U(i, j) + prefactor2 * s; 
					
					// Hamiltonian for interstitial (the same as for APW, but without a 0.5 factor due of different energy unit here) + Hamiltonian for muffin
					// the commented factor is for the alternative gamma commented in ComputeSlGammal
					H(i, j) = qiqjscalar /*0.5 * (qi2 + qj2)*/ * U(i, j) + prefactor2 * h;

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
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(H, S, Eigen::DecompositionOptions::EigenvaluesOnly);

			return es.eigenvalues();
		}


	protected:
		// technically the basis vectors are Ki + k, those here are only Ki
		const std::vector<Vector3D<double>>& m_basisVectors;
		const double m_R;
		const double prefactor;
		const double prefactor2;

		// see 6.43c for calculation details
		Eigen::MatrixXd U;

		// see 6.43a
		Eigen::MatrixXd S; 

		// see 6.44a
		Eigen::MatrixXd H;
	};

}

