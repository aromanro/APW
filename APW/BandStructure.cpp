#include <future>

#include "BandStructure.h"
#include "Hamiltonian.h"
#include "Numerov.h"

#include "ChemUtils.h"
#include "Pseudopotential.h"

namespace APW
{

	BandStructure::BandStructure(double a, double rmax)
		: nearestNeighbors(9),
		G2{ 0, 3, 4, 8, 11, 12, 16, 19, 20, 24 },
		m_a(a), m_Rmax(rmax)
	{
		// if the passed value was zero or negative, make them touching spheres
		if (m_Rmax <= 0)
			m_Rmax = sqrt(2.) * m_a / 4.;

		basisVectors.reserve(137);
	}


	bool BandStructure::GenerateBasisVectors(unsigned int nearestNeighborsNumber)
	{
		if (nearestNeighborsNumber < 2 || nearestNeighborsNumber > 10) return false;

		nearestNeighbors = nearestNeighborsNumber - 1;
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
					
					const double vectSquared =  vect * vect;

					if (vectSquared <= G2[nearestNeighbors]) // if it's under the cutoff length, add it
						basisVectors.push_back(vect);
				}

		return true;
	}


	void BandStructure::Initialize(std::vector<std::string> path, unsigned int nrPoints,  unsigned int nearestNeighborsNumber)
	{
		kpoints.clear();
		results.clear();

		kpoints.reserve(nrPoints);
		results.reserve(nrPoints);

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


	double VeffCu(double R)
	{
		const double R2 = R * R;
		const double R3 = R2 * R;
		const double R4 = R2 * R2;

		return 29. * exp(-2.3151241717834 * pow(R, 0.81266614122432) + 2.1984250222603E-2 * pow(R, 4.2246376280056))
					- 0.15595606773483 * R - 3.1350051440417E-3 * R2 + 5.1895222293006e-2 * R3 - 2.8027608685637E-2 * R4;
	}

	// not used but might turn useful with the proper params
	// use it as -StarkloffJoannopoulos(...) / r as the one for Cu
	double StarkloffJoannopoulos(int Z, double r, double lambda, double rc)
	{
		return Z * (1. - exp(-lambda * r)) / (1. + exp(-lambda * (r - rc)));
	}
	

	// the dumbest possible
	// constant for R < Rl
	// -Z/r otherwise (Z given by the number of valence electrons)
	// don't add - or /r for this one
	double VDumb(int Z, double R, double Rl, double C)
	{
		if (R < Rl) return C;

		return -Z / R;
	}


	double VAl(double R)
	{
		return VDumb(3, R, 2.675, 1.3905 * 0.5);
	}
	
	std::vector<std::vector<double>> BandStructure::Compute(const std::atomic_bool& terminate, const Options& options)
	{	
		const double minE = -0.2;
		const double maxE = 0.8;

		const double dE = 1E-3;
		const int numerovIntervals = 2000;
		const int numIntervals = static_cast<int>((maxE - minE) / dE);

		const int numerovGridNodes = numerovIntervals + 1;

		const double dr = m_Rmax / numerovIntervals;

		const size_t lMax = 5;

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
		
		//APW::Pseudopotential pseudopotential;
		//const bool success = pseudopotential.Load(Chemistry::ChemUtils::GetPseudopotentialFileForZ(13));
		//assert(success && pseudopotential.GetZ() == 13);

		APW::Potential potential;
		potential.m_potentialValues.resize(numerovGridNodes);
		for (int i = 0; i < numerovGridNodes; ++i)
		{
			//const double r = i * dr; // for uniform grid
			const double r = Rp * (exp(i * deltaGrid) - 1.);
			potential.m_potentialValues[i] = -VeffCu(r) / r;
			//potential.m_potentialValues[i] = pseudopotential.Value(r);
		}


		std::vector<std::future<void>> tasks(options.nrThreads);
		std::launch launchType = options.nrThreads == 1 ? std::launch::deferred : std::launch::async;

		int startPos = 0;
		int step = numIntervals / options.nrThreads;
		if (step < 1) step = 1; 
		int nextPos;


		for (int t = 0; t < options.nrThreads; ++t, startPos = nextPos)
		{			
			if (t == options.nrThreads - 1) nextPos = numIntervals;
			else nextPos = startPos + step;

			tasks[t] = std::async(launchType, [this, &potential, numerovGridNodes, &ratios, startPos, nextPos, minE, dE, lMax, numerovIntervals, deltaGrid, &terminate]()->void
				{
					//APW::Numerov<APW::NumerovFunctionRegularGrid> numerov(potential, 0, m_Rmax, numerovGridNodes);
					APW::Numerov<APW::NumerovFunctionNonUniformGrid> numerov(potential, deltaGrid, m_Rmax, numerovGridNodes);

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

		res.resize(kpoints.size());
	
		startPos = 0;
		step = kpoints.size() / options.nrThreads;
		if (step < 1) step = 1;

		for (int t = 0; t < options.nrThreads; ++t, startPos = nextPos)
		{
			if (t == options.nrThreads - 1) nextPos = kpoints.size();
			else nextPos = startPos + step;

			tasks[t] = std::async(launchType, [this, startPos, nextPos, numIntervals, minE, dE, lMax, &ratios, &res, &terminate]()->void
				{
					APW::Hamiltonian hamiltonian(basisVectors, m_Rmax, m_a, m_a * m_a * m_a / 4., lMax);

					// now, loop over k points
					for (int k = startPos; k < nextPos && !terminate; ++k)
					{
						hamiltonian.ComputeBC(kpoints[k]);

						// loop over all energies
						double olderDet = 0;
						double oldDet = 0;

						for (int posE = 0; posE < numIntervals && !terminate; ++posE)
						{
							const double E = minE + posE * dE;

							bool blowup = false;
							for (int l = 0; l <= lMax && !terminate; ++l)
							{
								if (abs(ratios[posE][l]) > 300 || isnan(ratios[posE][l]) || isinf(ratios[posE][l]))
								{
									blowup = true;
									break;
								}
							}
							if (blowup)
							{
								oldDet = olderDet = 0;
								continue;
							}


							hamiltonian.ComputeHamiltonian(E, ratios[posE]);

							const double det = hamiltonian.Determinant();

							if (posE > 0 && det * oldDet < 0) // change in sign
							{
								// linear interpolation
								const double val = E - dE * det / (det - oldDet);

								res[k].push_back(val);
							}
							else if (posE > 1 && abs(oldDet) < abs(olderDet) && abs(oldDet) < abs(det) && abs(oldDet) < 1E-15 && // went over a small minimum
								((olderDet < 0 && oldDet < 0 && det < 0) || (olderDet > 0 && oldDet > 0 && det > 0))) // all have the same sign, otherwise the sign change should be detected
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

								res[k].push_back(val);
							}

							olderDet = oldDet;
							oldDet = det;
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