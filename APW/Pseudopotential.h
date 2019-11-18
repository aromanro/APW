#pragma once

#include <string>
#include <vector>

namespace APW {

	class Pseudopotential
	{
	public:
		Pseudopotential();
		~Pseudopotential();

		bool Load(const std::string& name);

		double Value(double x) const;

		void Clear();

		unsigned int GetZ() const { return Z; }
		unsigned int GetZion() const { return Zion; }
		unsigned int ElectronsInCore() const { return Z - Zion; }
		double GetMaxRadius() const { return maxRadius; }
		bool IsValid() const { return valid; }

	protected:
		void ComputeSpline();

		double Interpolate(size_t interval, double x) const;
		size_t GetIndex(double x) const;

		unsigned int Z;
		unsigned int Zion;

		bool valid;

		double maxRadius;

		std::vector<double> position;
		std::vector<double> pseudopotential;

		std::vector<double> h;
		std::vector<double> z;
	};

}