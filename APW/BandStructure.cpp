#include "BandStructure.h"
#include "Hamiltonian.h"
#include "Numerov.h"

#include "ChemUtils.h"

namespace APW
{

	std::vector<std::vector<double>> BandStructure::Compute(const std::atomic_bool& terminate, const Options& options)
	{	
		const double minE = -0.2;
		const double maxE = 0.8;

		const double dE = 1E-3;
		const int numerovIntervals = 2000;
		const int numIntervals = static_cast<int>((maxE - minE) / dE);

		const int numerovGridNodes = numerovIntervals + 1;

		//const double dr = m_Rmax / numerovIntervals;

		const size_t lMax = 8;

		// the following two are needed for the non-uniform grid computations
		const double deltaGrid = 0.005;
		const double Rp = m_Rmax / (exp(numerovIntervals * deltaGrid) - 1.);

		std::vector<std::vector<double>> res;
		std::vector<std::vector<double>> ratios(numIntervals);

		// First, compute psi'/psi at muffin boundary, for each energy


		// this would load a pseudopotential file
		// there are some (local) pseudopotential files I found on GitHub, the Al file is ok for obtaining something similar with fig 6.3 from Computational Physics book
		// it's not as good as the Cu pseudopotential for APW, the match with https://www.materialscloud.org/discover/sssp/plot/efficiency/Al is worse
		// than that of Cu https://www.materialscloud.org/discover/sssp/plot/efficiency/Cu
		
		//Pseudopotential pseudopotential;
		//const bool success = pseudopotential.Load(Chemistry::ChemUtils::GetPseudopotentialFileForZ(13));
		//assert(success && pseudopotential.GetZ() == 13);

		Potential potential;
		InitializePotential(potential, numerovGridNodes, Rp, deltaGrid);

		std::vector<std::future<void>> tasks(options.nrThreads);
		std::launch launchType = options.nrThreads == 1 ? std::launch::deferred : std::launch::async;

		int startPos = 0;
		int step = static_cast<int>(ceil(static_cast<double>(numIntervals) / options.nrThreads));
		if (step < 1) step = 1; 
		int nextPos;


		for (int t = 0; t < options.nrThreads; ++t, startPos = nextPos)
		{			
			if (t == options.nrThreads - 1) nextPos = numIntervals;
			else nextPos = startPos + step;

			tasks[t] = std::async(launchType, [this, &potential, numerovGridNodes, &ratios, startPos, nextPos, minE, dE, lMax, numerovIntervals, deltaGrid, &terminate]()->void
				{
					//Numerov<NumerovFunctionRegularGrid> numerov(potential, 0, m_Rmax, numerovGridNodes);
					Numerov<NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, m_Rmax, numerovGridNodes);

					for (int posE = startPos; posE < nextPos && !terminate; ++posE)
					{
						const double E = minE + posE * dE;
						ratios[posE].resize(lMax + 1LL);
						for (int l = 0; l <= lMax && !terminate; ++l)
							// first parameter: pass m_Rmax for uniform grid, pass numerovIntervals for non-uniform (in this case the step is 1, so max radius is numerovIntervals)
							ratios[posE][l] = numerov.SolveSchrodinger(/*m_Rmax*/numerovIntervals, l, E, numerovIntervals);
					}
				}
			);
		}
		
		for (auto& task : tasks)
			task.get();

		// now compute the spectrum for each k point along the path
		ComputeBandstructure(tasks, res, ratios, numIntervals, minE, dE, lMax, terminate, options);

		return res;
	}

	void BandStructure::ComputeBandstructure(std::vector<std::future<void>>& tasks, std::vector<std::vector<double>>& res, std::vector<std::vector<double>>& ratios, int numIntervals, double minE, double dE, int lMax, const std::atomic_bool& terminate, const Options& options) const
	{
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

			tasks[t] = std::async(launchType, [this, startPos, nextPos, numIntervals, minE, dE, lMax, &ratios, &res, &terminate]()->void
				{
					Secular secular(basisVectors, m_Rmax, m_a * m_a * m_a / 4., lMax);

					// now, loop over k points
					for (int k = startPos; k < nextPos && !terminate; ++k)
					{
						secular.ComputeBC(kpoints[k]);

						// loop over all energies
						double olderDet = 0;
						double oldDet = 0;

						for (int posE = 0; posE < numIntervals && !terminate; ++posE)
						{
							const double E = minE + posE * dE;

							if (IsBlowup(ratios, posE, lMax, terminate))
							{
								oldDet = olderDet = 0;
								continue;
							}

							secular.ComputeHamiltonian(E, ratios[posE]);

							const double det = secular.Determinant();

							GetResult(res, k, E, posE, dE, det, oldDet, olderDet);

							olderDet = oldDet;
							oldDet = det;
						}
					}
				}
			);
		}

		for (auto& task : tasks)
			task.get();
	}

	void BandStructure::GetResult(std::vector<std::vector<double>>& res, int k, double E, double posE, double dE, double det, double oldDet, double olderDet) const
	{
		if (IsChangeInSign(posE, det, oldDet))
		{
			res[k].push_back(LinearInterpolation(E, dE, det, oldDet));
		}
		else if (posE > 1 && abs(oldDet) < abs(olderDet) && abs(oldDet) < abs(det) && abs(oldDet) < 1E-15 && // went over a small minimum
			((olderDet < 0 && oldDet < 0 && det < 0) || (olderDet > 0 && oldDet > 0 && det > 0))) // all have the same sign, otherwise the sign change should be detected
		{
			res[k].push_back(QuadraticInterpolation(E, dE, det, oldDet, olderDet));
		}
	}

	bool BandStructure::IsBlowup(const std::vector<std::vector<double>>& ratios, double posE, int lMax, const std::atomic_bool& terminate)
	{
		bool blowup = false;
		for (int l = 0; l <= lMax && !terminate; ++l)
		{
			if (isnan(ratios[posE][l]) || isinf(ratios[posE][l]) || abs(ratios[posE][l]) > 300)
			{
				blowup = true;
				break;
			}
		}

		return blowup;
	}

	bool BandStructure::IsChangeInSign(double posE, double det, double oldDet)
	{
		return posE > 0 && det * oldDet < 0;
	}

	double BandStructure::LinearInterpolation(double E, double dE, double det, double oldDet)
	{
		const double val = E - dE * det / (det - oldDet);

		if (!isnan(val) && !isinf(val) && val >= E - dE && val <= E)
			return val;

		return E - 0.5 * dE;		
	}

	double BandStructure::QuadraticInterpolation(double E, double dE, double det, double oldDet, double olderDet)
	{
		// quadratic interpolation:

		// write the interpolation polynomial out of the three points:
		// y(x) = y0(x-x1)(x-x2)/((x0-x1)(x0-x2)) + y1(x-x0)(x-x2)/((x1-x0)(x1-x2))+ y2(x-x0)(x-x1)/((x2-x0)(x2-x1)) 
		// points: (x0, y0) = (E - 2*dE, olderDet), (x1, y1) = (E - dE, oldDet), (x2, y2) = (E, det)
		// the minimum is where y'(x) = 0 so calculate the derivative of the polynomial, make it equal with 0, solve for x and that's the val

		// TODO: check it, I derived it too fast, might have some mistakes
		const double coef1 = olderDet * 0.5;
		const double coef2 = -oldDet;
		const double coef3 = det * 0.5;

		//Check:
		// the interpolation polynomial is:
		// y(x) = coef1 (x-x1)(x-x2) / dE^2 + coef2 (x-x0)(x-x2) / dE^2+ coef3 (x-x0)(x-x1) / dE^2 

		//const double val = E - dE;
		const double val = E - 0.5 * dE * (coef1 + 2 * coef2 + 3 * coef3) / (coef1 + coef2 + coef3);

		if (!isnan(val) && !isinf(val) && val >= E - 2 * dE && val <= E)
			return val;

		return val - dE;
	}

}