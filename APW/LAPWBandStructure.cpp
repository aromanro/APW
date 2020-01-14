#include "LAPWBandStructure.h"

#include <future>

#include "LAPWBandStructure.h"
#include "LAPWHamiltonian.h"
#include "Numerov.h"
#include "Integral.h"

#include "ChemUtils.h"
#include "Pseudopotential.h"

namespace LAPW
{

	void BandStructure::NormalizeUniform(std::vector<double>& Psi, double h)
	{
		const size_t size = Psi.size();
		std::vector<double> result2(size);
		for (int i = 0; i < size; ++i)
			result2[i] = Psi[i] * Psi[i];

		const double integralForSquare = Integral::Boole(h, result2);
		const double norm = sqrt(integralForSquare);

		for (int i = 0; i < size; ++i)
			Psi[i] /= norm;
	}



	void BandStructure::NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid)
	{
		const size_t size = Psi.size();
		std::vector<double> result2(size);
		for (int i = 0; i < size; ++i)
		{
			result2[i] = Psi[i] * Psi[i];

			// if nonuniform, do the change for dr	
			const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
			result2[i] *= cnst;
		}

		const double integralForSquare = Integral::Boole(1, result2); // for nonuniform case the step is 1
		const double norm = sqrt(integralForSquare);

		for (int i = 0; i < size; ++i)
			Psi[i] /= norm;
	}

	std::vector<std::vector<double>> BandStructure::Compute(const std::atomic_bool& terminate, const Options& options)
	{
		const int numerovIntervals = 2000;// 16000;
		const int numerovGridNodes = numerovIntervals + 1;
		//const double dr = m_Rmax / numerovIntervals; // for uniform

		const size_t lMax = 5;

		// the following two are needed for the non-uniform grid computations
		const double deltaGrid = 0.005;
		const double Rp = m_Rmax / (exp(numerovIntervals * deltaGrid) - 1.);

		std::vector<std::vector<double>> res;


		APW::Potential potential;
		potential.m_potentialValues.resize(numerovGridNodes);
		for (int i = 0; i < numerovGridNodes; ++i)
		{
			//const double r = i * dr; // for uniform grid
			const double r = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = -APW::Pseudopotential::VeffCu(r) / r;
			//potential.m_potentialValues[i] = pseudopotential.Value(r);
		}


		std::vector<Values> vals(lMax + 1);
		vals[0].El = vals[1].El = vals[2].El = 0.2;
		vals[3].El = 0.3;
		vals[4].El = 0.4;
		vals[5].El = 0.6;

		const double R2 = m_Rmax * m_Rmax;

		for (unsigned int l = 0; l <= lMax; ++l)
		{
			const double El = vals[l].El;

			// compute wavefunction
			//APW::Numerov<APW::NumerovFunctionRegularGrid> numerov(potential, 0, m_Rmax, numerovGridNodes);
			APW::Numerov<APW::NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, m_Rmax, numerovGridNodes);

			std::vector<double> u = numerov.SolveSchrodingerFull(/*m_Rmax*/numerovIntervals, l, El, numerovIntervals);

			// normalize it
			//NormalizeUniform(u, dr);
			NormalizeNonUniform(u, Rp, deltaGrid);

			// for check, numerically compute the derivative of energy in a different way
			/*
			std::vector<double> un = numerov.SolveSchrodingerFull(numerovIntervals, l, El + 0.0001, numerovIntervals);
			NormalizeNonUniform(un, Rp, deltaGrid);
			for (int i = 0; i < u.size(); ++i)
				un[i] = (un[i] - u[i]) / 0.0001; // (wavefunction(E + dE) - wavefunction(E)) / dE
			*/

			const size_t size = u.size();
			const size_t lastPos = size - 1;

			vals[l].Wavefunction = u[lastPos] / m_Rmax; // Rl = u / r

			const double derivStep = numerov.function.GetDerivativeStep(numerovIntervals, 1/*dr*/);

			// compute its radial derivative
			const double up = (u[lastPos] - u[lastPos - 1]) / derivStep;
			vals[l].RadialDerivative = up / m_Rmax - u[lastPos] / R2;

			// compute the energy derivative of the wavefunction
			std::vector<double> udot = numerov.SolveGeneral(u, /*m_Rmax*/numerovIntervals, l, El, numerovIntervals);

			// the equation is inhomogeneous, add a particular solution of the homogeneous eqn, alpha * u
			// get alpha from 6.48 condition
			// this way udot is orthogonalized with u, by subtracting the projection of udot on u (that is, alpha * u)

			// fill it with u * udot
			std::vector<double> uudot(size);
			for (int i = 0; i < size; ++i)
			{
				// if nonuniform, do the change for dr
				const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
				uudot[i] = u[i] * udot[i] * cnst;
			}

			const double alpha = Integral::Boole(1/*dr*/, uudot);
			for (int i = 0; i < size; ++i)
				udot[i] -= alpha * u[i];

			// udot is correct, I get it basically the same by the two methods

			//udot = un; // use the 'numerical' derivative one, for test

			// check:
			/*
			std::vector<double> udotu(size);
			for (int i = 0; i < size; ++i)
			{
				const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
				udotu[i] = udot[i] * u[i] * cnst;
			}
			const double sum = Integral::Boole(1, udotu); // should be very close to zero
			*/

			// the formulae were derived with Rydberg atomic units, but the program uses Hartrees, so convert the derivatives to use Rydbergs for formulae
			// d/d(2E) = 0.5 d/dE
			// so the two multiplications with 0.5 that follow are for this reason

			vals[l].EnergyDerivative = 0.5 * udot[lastPos] / m_Rmax;
			
			// now derivative of both
			const double udotp = (udot[lastPos] - udot[lastPos - 1]) / derivStep;
			vals[l].BothDerivative = 0.5 * (udotp / m_Rmax - udot[lastPos] / R2);

			// 6.49
			for (int i = 0; i < size; ++i)
			{
				// if nonuniform, do the change for dr	
				const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
				udot[i] *= udot[i] * cnst;
			}
				
			// the 0.25 is for the same reason as 0.5 above
			vals[l].Nl = 0.25 * Integral::Boole(1/*dr*/, udot);

			// should be 1, see 6.50 - but if you choose to go further with relations derived for Hartree atomic units, should be 2 (this comes out of the 1/2 of the kinetic term of the Schrodinger eqn).
			//const double val = m_Rmax * m_Rmax * (vals[l].RadialDerivative * vals[l].EnergyDerivative - vals[l].Wavefunction * vals[l].BothDerivative); 		
		}


		// now solve the generalized eigenvalue problem for each k

		std::vector<std::future<void>> tasks(options.nrThreads);
		std::launch launchType = options.nrThreads == 1 ? std::launch::deferred : std::launch::async;


		res.resize(kpoints.size());

		int startPos = 0;
		int step = static_cast<int>(ceil(static_cast<double>(kpoints.size()) / options.nrThreads));
		if (step < 1) step = 1;
		int nextPos;

		for (int t = 0; t < options.nrThreads; ++t, startPos = nextPos)
		{
			if (t == options.nrThreads - 1) nextPos = kpoints.size();
			else nextPos = startPos + step;

			if (nextPos > kpoints.size()) nextPos = kpoints.size();

			tasks[t] = std::async(launchType, [this, startPos, nextPos, lMax, &vals, &res, &terminate]()->void
				{					
					Hamiltonian hamiltonian(basisVectors, m_Rmax, m_a* m_a* m_a / 4.);

					// now, loop over k points
					for (int k = startPos; k < nextPos && !terminate; ++k)
					{
						hamiltonian.Compute(kpoints[k], vals);
						const Eigen::VectorXd v = hamiltonian.GetEnergies();
						for (int i = 0; i < v.rows(); ++i)
						{
							const double val = 0.5 * v(i); // the formulae were with Rydbergs, convert back to Hartrees
							if (val > 0.8) break;

							res[k].push_back(val); 
						}

					}
				}
			);
		}

		for (auto& task : tasks)
			task.get();

		return std::move(res);
	}

}
