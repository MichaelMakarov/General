#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		/// <summary>
		/// Class represents a polynom as a single variable function
		/// </summary>
		class Polynom
		{
		private:
			/// <summary>
			/// coefficients of the polynom
			/// </summary>
			std::vector<double> _data;

		public:
			Polynom(const size_t degree = 0) : _data{ std::vector<double>(degree + 1) } {}
			/// <summary>
			/// Initializing by an array of the coefficients
			/// </summary>
			/// <param name="coefficients"> - an array of the coefficients from zero to the highest degree</param>
			template<size_t size> Polynom(const double(&coefficients)[size]) : _data{ std::vector<double>(size) }
			{
				std::memcpy(_data.data(), coefficients, size * sizeof(double));
			}
			Polynom(const Polynom& p) : _data{ p._data } {}
			Polynom(Polynom&& p) noexcept : _data{ std::move(p._data) } {}
			~Polynom() = default;

			Polynom& operator = (const Polynom& p) { _data = p._data; return *this; }
			Polynom& operator = (Polynom&& p) noexcept { _data = std::move(p._data); return *this; }

			size_t degree() const { return _data.size() - 1; }
			double& operator [] (const size_t index) { return _data[index]; }
			const double& operator [] (const size_t index) const { return _data[index]; }

			double operator () (const double x) const
			{
				double result{ _data[0] };
				double mult{ x };
				for (size_t i = 1; i < _data.size(); ++i) {
					result += mult * _data[i];
					mult *= x;
				}
				return result;
			}
		};
		/// <summary>
		/// Creating a polynom as a solution of the LST
		/// </summary>
		/// <param name="X"> - a vector</param>
		/// <param name="Y"> - a vector</param>
		/// <param name="degree"> - a degree of the polynom</param>
		/// <returns></returns>
		Polynom create_polynom(const VectorDyn& X, const VectorDyn& Y, const size_t degree);
	}
}