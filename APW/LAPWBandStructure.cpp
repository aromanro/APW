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
		std::vector<double> result2(Psi.size());
		for (int i = 0; i < result2.size(); ++i)
			result2[i] = Psi[i] * Psi[i];

		const double integralForSquare = Integral::Boole(h, result2);
		const double norm = sqrt(integralForSquare);

		for (int i = 0; i < Psi.size(); ++i)
			Psi[i] /= norm;
	}



	void BandStructure::NormalizeNonUniform(std::vector<double>& Psi, double Rp, double deltaGrid)
	{
		std::vector<double> result2(Psi.size());
		for (int i = 0; i < result2.size(); ++i)
		{
			// if nonuniform, convert the function back!!!!!
			//Psi[i] *= exp(i * deltaGrid * 0.5); // it's already converted in numerov

			result2[i] = Psi[i] * Psi[i];

			// if nonuniform, do the change for dr	
			const double cnst = Rp * deltaGrid * exp(deltaGrid * i);
			result2[i] *= cnst;
		}

		const double integralForSquare = Integral::Boole(1, result2); // for nonuniform case the step is 1
		const double norm = sqrt(integralForSquare);

		for (int i = 0; i < Psi.size(); ++i)
			Psi[i] /= norm;
	}


	std::vector<std::vector<double>> BandStructure::Compute(const std::atomic_bool& terminate, const Options& options)
	{
		const int numerovIntervals = 2000;
		const int numerovGridNodes = numerovIntervals + 1;
		const double dr = m_Rmax / numerovIntervals;

		const size_t lMax = 5;

		// the following two are needed for the non-uniform grid computations
		const double deltaGrid = 0.005;
		const double Rp = m_Rmax / (exp(numerovIntervals * deltaGrid) - 1.);

		std::vector<std::vector<double>> res;


		APW::Potential potential;
		potential.m_potentialValues.resize(numerovGridNodes);
		for (int i = 0; i < numerovGridNodes; ++i)
		{
			const double r = i * dr; // for uniform grid
			//const double r = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = -APW::Pseudopotential::VeffCu(r) / r;
			//potential.m_potentialValues[i] = pseudopotential.Value(r);
		}


		std::vector<Values> vals(lMax + 1);
		vals[0].El = vals[1].El = vals[2].El = 0.2;
		vals[3].El = 0.3;
		vals[4].El = 0.4;
		vals[5].El = 0.5;

		// TODO: must check this, probably has mistakes
		// tried it, it does not work, something is still not ok

		// I'm not yet sure which derivatives are used, the book seems wrong

		const double R2 = m_Rmax * m_Rmax;

		for (unsigned int l = 0; l <= lMax; ++l)
		{
			const double El = vals[l].El;

			// compute wavefunction
			APW::Numerov<APW::NumerovFunctionRegularGrid> numerov(potential, 0, m_Rmax, numerovGridNodes);
			//APW::Numerov<APW::NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, m_Rmax, numerovGridNodes);
			std::vector<double> u = numerov.SolveSchrodingerFull(m_Rmax/*numerovIntervals*/, l, El, numerovIntervals);

			// normalize it
			//NormalizeNonUniform(u, Rp, deltaGrid);
			NormalizeUniform(u, dr);

			const size_t lastPos = u.size() - 1;

			vals[l].Wavefunction = u[lastPos]; // / m_Rmax; // Rl = u / r


			const double derivStep = numerov.function.GetDerivativeStep(numerovIntervals, 1);

			// compute its radial derivative
			const double up = (u[lastPos] - u[lastPos - 1]) / derivStep;
			vals[l].RadialDerivative = up; // / m_Rmax - u[lastPos] / R2;

			// compute the energy derivative of the wavefunction
			std::vector<double> udot = numerov.SolveGeneral(u, m_Rmax/*numerovIntervals*/, l, El, numerovIntervals);
			// the equation is inhomogeneous, add a particular solution of the homogeneous eqn, alpha * u
			// get alpha from 6.48 condition

			// fill it with u * udot
			std::vector<double> uudot(udot.size());
			for (int i = 0; i < uudot.size(); ++i)
				uudot[i] = u[i] * udot[i];

			const double alpha = -Integral::Boole(1, uudot);

			for (int i = 0; i < udot.size(); ++i)
				udot[i] += alpha * u[i];

			vals[l].EnergyDerivative = udot[lastPos]; // / m_Rmax;
			
			// now derivative of both
			const double udotp = (udot[lastPos] - udot[lastPos - 1]) / derivStep;
			vals[l].BothDerivative = udotp; // / m_Rmax - udot[lastPos] / R2;

			// 6.49
			for (int i = 0; i < udot.size(); ++i)
				udot[i] *= udot[i] * R2;
				
			vals[l].Nl = Integral::Boole(1, udot);
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
					Hamiltonian hamiltonian(basisVectors, m_Rmax, m_a, m_a * m_a * m_a / 4.);

					// now, loop over k points
					for (int k = startPos; k < nextPos && !terminate; ++k)
					{
						hamiltonian.Compute(k, vals);
						const Eigen::VectorXd v = hamiltonian.GetEnergies();
						for (int i = 0; i < v.rows(); ++i)
							res[k].push_back(v(i));

					}
				}
			);
		}

		for (auto& task : tasks)
			task.get();

		return std::move(res);
	}

}
