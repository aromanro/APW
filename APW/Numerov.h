#pragma once

#include <vector>
#include <algorithm> 

namespace APW {

	class Potential
	{
	public:
		inline double operator()(size_t posIndex) const { return m_potentialValues[posIndex]; };

		std::vector<double> m_potentialValues;
	};


	class NumerovFunctionRegularGrid
	{
	public:
		NumerovFunctionRegularGrid(const Potential& pot, double /*delta*/, double /*Rmax*/, size_t /*numPoints*/) : m_pot(pot) {}

		inline double GetEffectivePotential(unsigned int l, double position, size_t posIndex) const
		{
			return m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;
		}

		inline double operator()(unsigned int l, double E, double position, size_t posIndex) const
		{
			const double effectivePotential = GetEffectivePotential(l, position, posIndex);

			return  2. * (effectivePotential - E);
		}

		inline static double GetBoundaryValueFar(double position, double E)
		{
			return exp(-position * sqrt(2. * abs(E)));
		}

		inline static double GetBoundaryValueZero(double position, unsigned int l)
		{
			return pow(position, static_cast<size_t>(l) + 1);
		}

		inline static double GetMaxRadiusIndex(double E, size_t maxIndex, double stepSize)
		{
			return std::min(GetMaxRadius(E, maxIndex) / stepSize, static_cast<double>(maxIndex));
		}

		inline static double GetDerivativeStep(int /*posIndex*/, double h)
		{
			return h;
		}

		inline static double GetMaxRadius(double E, size_t /*maxIndex*/)
		{
			return 200. / sqrt(2. * abs(E));
		}

		inline static double GetWavefunctionValue(size_t /*posIndex*/, double value)
		{
			return value;
		}

		inline double GetSrcAdjustedValue(size_t /*posIndex*/, double value) const
		{
			return value;
		}

		inline static bool IsUniform()
		{
			return true;
		}

	private:
		const Potential& m_pot;
	};


	class NumerovFunctionNonUniformGrid
	{
	public:
		NumerovFunctionNonUniformGrid(const Potential& pot, double delta, double Rmax, size_t numPoints)
			: m_pot(pot), m_delta(delta)
		{
			Rp = Rmax / (exp((numPoints - 1) * delta) - 1);
			const double Rp2 = Rp * Rp;

			twodelta = 2. * m_delta;
			const double delta2 = m_delta * m_delta;

			Rp2delta2 = Rp2 * delta2;
			delta2p4 = delta2 * 0.25;
		}

		inline double GetEffectivePotential(unsigned int l, double position, size_t posIndex) const
		{
			position = GetPosition(posIndex); // the passed value is ignored, use the real one

			return m_pot(posIndex) + l * (l + 1.) / (position * position) * 0.5;
		}

		inline double operator()(unsigned int l, double E, double position, size_t posIndex) const
		{
			const double effectivePotential = GetEffectivePotential(l, position, posIndex);

			return  2. * (effectivePotential - E) * Rp2delta2 * exp(posIndex * twodelta) + delta2p4;
		}

		inline double GetBoundaryValueFar(double position, double E) const
		{
			const double realPosition = GetPosition(static_cast<int>(position));

			return exp(-realPosition * sqrt(2. * abs(E)) - static_cast<int>(position) * m_delta * 0.5);
		}

		inline double GetBoundaryValueZero(double position, unsigned int l) const
		{
			const int posInd = static_cast<int>(position);
			const double realPosition = GetPosition(posInd);

			return pow(realPosition, static_cast<size_t>(l) + 1) * exp(-static_cast<int>(position) * m_delta * 0.5);
		}


		inline double GetMaxRadiusIndex(double E, size_t maxIndex, double /*stepSize*/) const
		{
			double val = GetBoundaryValueFar(static_cast<double>(maxIndex), E);
			if (val > MaxRadiusLimit) static_cast<double>(maxIndex);

			size_t minIndex = 1;
			while (maxIndex - minIndex > 1)
			{
				const size_t midIndex = (maxIndex + minIndex) / 2;
				val = GetBoundaryValueFar(static_cast<double>(midIndex), E);
				if (val < MaxRadiusLimit)
					maxIndex = midIndex;
				else
					minIndex = midIndex;
			}

			return static_cast<double>(maxIndex);
		}

		inline double GetMaxRadius(double E, size_t maxIndex) const
		{
			double val = GetBoundaryValueFar(static_cast<double>(maxIndex), E);
			if (val > MaxRadiusLimit)
			{
				const double position = Rp * (exp(maxIndex * m_delta) - 1.);

				return position;
			}

			size_t minIndex = 1;
			while (maxIndex - minIndex > 1)
			{
				const size_t midIndex = (maxIndex + minIndex) / 2;
				val = GetBoundaryValueFar(static_cast<double>(midIndex), E);
				if (val < MaxRadiusLimit)
					maxIndex = midIndex;
				else
					minIndex = midIndex;
			}
			
			return Rp * (exp(maxIndex * m_delta) - 1.);
		}

		inline double GetDerivativeStep(int posIndex, double /*h*/) const
		{
			return Rp * exp(posIndex * m_delta) * (1. - exp(-m_delta));
		}

		inline double GetWavefunctionValue(size_t posIndex, double value) const
		{
			return exp(0.5 * posIndex * m_delta) * value;
		}

		// NOTE: This makes the 'general' solver below less general
		// as it is, it's targeted to the eqn (H-E)udot = u (6.46 in the book)
		// if you want it to be more general, use Rp * exp(posIndex * m_delta)
		// and adjust the 'source' term before passing it to the solver

		inline double GetSrcAdjustedValue(size_t posIndex, double value) const
		{
			return value * Rp2delta2 * exp(1.5 * posIndex * m_delta);
		}


		inline double GetRp() const { return Rp; }
		inline double GetDelta() const { return m_delta; }

		inline static bool IsUniform()
		{
			return false;
		}
	protected:
		inline double GetPosition(size_t posIndex) const
		{
			return Rp * (exp(posIndex * m_delta) - 1.);
		}

		const Potential& m_pot;

		const double m_delta;

		double Rp;
		double twodelta;
		double delta2p4;
		double Rp2delta2;

		static constexpr double MaxRadiusLimit = 1E-200;
	};


	template<class NumerovFunction> class Numerov
	{
	public:
		Numerov(const Potential& pot, double delta = 0, double Rmax = 0, size_t numPoints = 0) : function(pot, delta, Rmax, numPoints), h(1), h2(1), h2p12(1. / 12.) {}

		inline double SolveSchrodinger(double endPoint, unsigned int l, double E, long int steps)
		{
			if (NumerovFunction::IsUniform())
			{
				h = endPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadius(E, steps));
				steps = static_cast<long int>(endPoint / h);
			}
			else
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadiusIndex(E, steps, 1));
				steps = static_cast<long int>(endPoint);
			}

			double position = 0;
			double solution = 0;
			double wprev = 0;

			double oldSolution = 0;

			position = h;
			solution = function.GetBoundaryValueZero(position, l);
			double funcVal = function(l, E, position, 1);
			double w = (1. - h2p12 * funcVal) * solution;

			for (long int i = 2; i <= steps; ++i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);

				oldSolution = solution;
				solution = getU(w, funcVal);

				if (abs(solution) == std::numeric_limits<double>::infinity() || isnan(solution))
					return std::numeric_limits<double>::infinity();
			}
			
			const double realSolution = function.GetWavefunctionValue(steps, solution);
			const double prevSolution = function.GetWavefunctionValue(steps - 1ULL, oldSolution);

			return (realSolution - prevSolution) / (function.GetDerivativeStep(steps, h) * realSolution);
		}


		// should be needed in case of wanting to implement LAPW

		inline std::vector<double> SolveSchrodingerFull(double endPoint, unsigned int l, double E, long int steps)
		{
			const long int highLimit = steps + 1;
			std::vector<double> Psi(highLimit);

			if (NumerovFunction::IsUniform())
			{
				h = endPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadius(E, steps));
				steps = static_cast<long int>(endPoint / h);
			}
			else
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadiusIndex(E, steps, 1));
				steps = static_cast<long int>(endPoint);
			}

			double position = 0;
			double solution = 0;
			double wprev = 0;

			Psi[0] = function.GetWavefunctionValue(0, solution);

			position = h;
			solution = function.GetBoundaryValueZero(position, l);

			Psi[1] = function.GetWavefunctionValue(1, solution);

			double funcVal = function(l, E, position, 1);
			double w = (1. - h2p12 * funcVal) * solution;

			for (long int i = 2; i <= steps; ++i)
			{
				const double wnext = 2. * w - wprev + h2 * solution * funcVal;

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);

				solution = getU(w, funcVal);

				Psi[i] = function.GetWavefunctionValue(i, solution); // already adjusted for the case of the non uniform grid
			}

			return Psi;
		}


		// it's not actually the 'general' solution (for example assumes solution and src 0 at r = 0), but it adds a 'source' term
		// if 'function' would return zero, this would solve the Poisson equation with the source 'src', whence the name
		// it should be actually used to obtain the energy derivative of the wavefunction, needed in LAPW (see eq 6.46 in the Computational Physics book)
		// so the 'source' would be actually the wavefunction (but this, as above, deals with u, not with R; R = u / r)
		// to be noted that src should be multiplied by 2 because the kinetic term in Schrodinger has 1/2 in it
		// for the same reason the operator() for the 'function' has a multiplication with 2

		// NOTE: it's made less general for the non-uniform grid, targeted at eqn (H-E)udot = u (6.46 in the book)
		// see above the comment in NumerovFunctionNonUniformGrid class

		inline std::vector<double> SolveGeneral(const std::vector<double>& src, double endPoint, unsigned int l, double E, long int steps)
		{
			const long int highLimit = steps + 1;
			std::vector<double> Psi(highLimit);

			if (NumerovFunction::IsUniform())
			{
				h = endPoint / steps;
				h2 = h * h;
				h2p12 = h2 / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadius(E, steps));
				steps = static_cast<long int>(endPoint / h);
			}
			else
			{
				h = 1;
				h2 = 1;
				h2p12 = 1. / 12.;

				endPoint = std::min(endPoint, function.GetMaxRadiusIndex(E, steps, 1));
				steps = static_cast<long int>(endPoint);
			}


			double position = 0;
			double solution = 0;
			double srcVal = function.GetSrcAdjustedValue(0, 2. * src[0]);
			double wprev = h2p12 * srcVal;

			Psi[0] = function.GetWavefunctionValue(0, solution);

			position = h;
			solution = function.GetBoundaryValueZero(position, l);

			Psi[1] = function.GetWavefunctionValue(1, solution);

			double funcVal = function(l, E, position, 1);
			srcVal = function.GetSrcAdjustedValue(1, 2. * src[1]);

			double w = (1. - h2p12 * funcVal) * solution + h2p12 * srcVal;

			for (long int i = 2; i <= steps; ++i)
			{
				const double wnext = 2. * w - wprev + h2 * (solution * funcVal - srcVal);

				position = h * i;

				wprev = w;
				w = wnext;

				funcVal = function(l, E, position, i);
				srcVal = function.GetSrcAdjustedValue(i, 2. * src[i]);

				solution = getU(w, funcVal, srcVal);

				Psi[i] = function.GetWavefunctionValue(i, solution); 
			}

			return Psi;
		}

		NumerovFunction function;

	private:
		// 2.13
		inline double getU(double w, double funcVal) const
		{
			return w / (1. - h2p12 * funcVal);
		}

		inline double getU(double w, double funcVal, double srcVal) const
		{
			return (w - h2p12 * srcVal) / (1. - h2p12 * funcVal);
		}

		double h;
		double h2;
		double h2p12;
	};

}


