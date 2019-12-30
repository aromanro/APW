#include "LAPWBandStructure.h"



#include <future>

#include "BandStructure.h"
#include "Hamiltonian.h"
#include "Numerov.h"

#include "ChemUtils.h"
#include "Pseudopotential.h"

namespace LAPW
{
	// This is empty for now, I'll leave it here just in case I want to add LAPW implementation

	std::vector<std::vector<double>> BandStructure::Compute(const std::atomic_bool& terminate, const Options& options)
	{
		const double minE = -0.2;
		const double maxE = 0.8;

		const int numerovIntervals = 2000;
		const int numerovGridNodes = numerovIntervals + 1;
		const double dr = m_Rmax / numerovIntervals;

		const size_t lMax = 5;

		// the following two are needed for the non-uniform grid computations
		const double deltaGrid = 0.005;
		const double Rp = m_Rmax / (exp(numerovIntervals * deltaGrid) - 1.);

		std::vector<std::vector<double>> res;


		return res;
	}

}
