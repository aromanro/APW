#pragma once

// the referred formulae are from the book
// Computational Physics by J M Thijssen
// isbn: 9780521833462, https://doi.org/10.1017/CBO9781139171397


#define _USE_MATH_DEFINES
//#include <math.h>

// you can use boost for the same purpose if spherical Bessel functions are not available
#define USE_BETTER_BESSEL

#ifdef USE_BETTER_BESSEL
#include <cmath>
#endif

namespace SpecialFunctions
{
	// for now I'll let the ones implemented here

	// the recursion formulae are more clear 
	// I first implemented the code directly with them
	// I'll let them here commented out


	class Bessel
	{
	public:
		static double j(unsigned int l, double x)
		{
#ifdef USE_BETTER_BESSEL
			return std::sph_bessel(l, x);
#else
			/*
			if (0 == l) return j0(x);
			else if (1 == l) return j1(x);

			return (2. * l - 1.) / x * j(l - 1, x) - j(l - 2, x);
			*/

			if (0 == l) return sin(x) / x;
			else if (1 == l) return  sin(x) / (x * x) - cos(x) / x;

			const double sinxx = sin(x) / x;
			const double cosxx = cos(x) / x;

			double j0 = sinxx;
			double j1 = sinxx / x - cosxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double j2 = (2. * i - 1.) / x * j1 - j0;

				j0 = j1;
				j1 = j2;
			}

			return j1;
#endif
		}


		static double n(unsigned int l, double x)
		{
#ifdef USE_BETTER_BESSEL
			return std::sph_neumann(l, x);
#else
			/*
			if (0 == l) return n0(x);
			else if (1 == l) return n1(x);

			return (2. * l - 1.) / x * n(l - 1, x) - n(l - 2, x);
			*/

			if (0 == l) return -cos(x) / x;
			else if (1 == l) return -cos(x) / (x * x) - sin(x) / x;

			const double sinxx = sin(x) / x;
			const double cosxx = cos(x) / x;

			double n0 = -cosxx;
			double n1 = -cosxx / x - sinxx;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double n2 = (2. * i - 1.) / x * n1 - n0;

				n0 = n1;
				n1 = n2;
			}

			return n1;
#endif
		}

		/*
		protected:
			inline static double j0(double x) { return sin(x) / x; }
			inline static double j1(double x) { return sin(x) / (x * x) - cos(x) / x; }

			inline static double n0(double x) { return -cos(x) / x; }
			inline static double n1(double x) { return -cos(x) / (x * x) - sin(x) / x; }
		*/
	};

	// Legendre polynomials
	class Legendre
	{
	public:
		static double p(unsigned int l, double x)
		{

#ifdef USE_BETTER_LEGENDRE
			return std::legendre(l, x); // WARNING: for some reason this does not work!
#else
			if (0 == l) return 1.;
			else if (1 == l) return x;

			//return ((2. * l - 1.) * x * p(l - 1, x) - (l - 1) * p(l - 2, x)) / l;

			double p0 = 1.;
			double p1 = x;

			for (unsigned int i = 2; i <= l; ++i)
			{
				const double p2 = ((2. * i - 1.) * x * p1 - (i - 1.) * p0) / i;
				p0 = p1;
				p1 = p2;
			}

			return p1;
#endif
		}
	};

}
