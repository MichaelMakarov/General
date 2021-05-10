#pragma once
#include "Vector.h"
#include <functional>

namespace general
{
	namespace math
	{
		/// <summary>
		/// Class represents a polynom as a single variable function
		/// </summary>
		class PowerPolynomial
		{
		private:
			/// <summary>
			/// coefficients of the polynom
			/// </summary>
			std::vector<double> _data;

		public:
			explicit PowerPolynomial(const size_t degree = 0) : _data{ std::vector<double>(degree + 1) } {}
			/// <summary>
			/// Initializing by an array of the coefficients
			/// </summary>
			/// <param name="coefficients"> - an array of the coefficients from zero to the highest degree</param>
			template<size_t size> explicit PowerPolynomial(const double(&coefficients)[size]) : _data{ std::vector<double>(size) }
			{
				std::memcpy(_data.data(), coefficients, size * sizeof(double));
			}
			PowerPolynomial(const PowerPolynomial& p) : _data{ p._data } {}
			PowerPolynomial(PowerPolynomial&& p) noexcept : _data{ std::move(p._data) } {}
			~PowerPolynomial() = default;

			PowerPolynomial& operator = (const PowerPolynomial& p) { _data = p._data; return *this; }
			PowerPolynomial& operator = (PowerPolynomial&& p) noexcept { _data = std::move(p._data); return *this; }

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
		PowerPolynomial create_polynom(const Vector& x, const Vector& Y, const size_t degree);

		/// <summary>
		/// Class represents Legendre polynomial Pn
		/// </summary>
		class LegendrePolynomial
		{
		protected:
			size_t _degree;
			bool _normalized;

		public:
			/// <summary>
			/// Creating a Legendre polynomial
			/// </summary>
			/// <param name="degree"> - a degree</param>
			/// <param name="normalized"> - a flag whether polynomial is normalized or not</param>
			LegendrePolynomial(
				const size_t degree = 0,
				const bool normalized = false);
			LegendrePolynomial(const LegendrePolynomial& p) noexcept = default;
			LegendrePolynomial(LegendrePolynomial&& p) noexcept = default;
			~LegendrePolynomial() noexcept = default;

			LegendrePolynomial& operator = (const LegendrePolynomial& p) noexcept = default;
			LegendrePolynomial& operator = (LegendrePolynomial&& p) noexcept = default;

			size_t degree() const { return _degree; }

			bool normalized() const { return _normalized; }

			virtual double operator () (const double x) const;
		};
		/// <summary>
		/// Class represents a Legendre function Pnm
		/// </summary>
		class LegendreFunction : public LegendrePolynomial
		{
		private:
			size_t _derivation;

		public:
			/// <summary>
			/// Creating a Legendre function
			/// </summary>
			/// <param name="degree"> - a degree of the polynom</param>
			/// <param name="derivation"> - a derivation of the polynom</param>
			/// <param name="normalized"> - a flag whether the function is normalized or not</param>
			LegendreFunction(
				const size_t degree = 0,
				const size_t derivation = 0,
				const bool normalized = false) :
				LegendrePolynomial(degree, normalized),
				_derivation{ derivation }
			{}
			LegendreFunction(const LegendrePolynomial& p) noexcept;
			LegendreFunction(const LegendreFunction& f) noexcept = default;
			LegendreFunction(LegendreFunction&& f) noexcept = default;
			~LegendreFunction() noexcept = default;

			LegendreFunction& operator = (const LegendrePolynomial& p) noexcept;
			LegendreFunction& operator = (const LegendreFunction& f) noexcept = default;
			LegendreFunction& operator = (LegendreFunction&& f) noexcept = default;

			size_t derivation() const { return _derivation; }

			bool normalized() const { return _normalized; }

			double operator () (const double x) const;
		};
		/// <summary>
		/// Class represents Newtonian polynomial
		/// </summary>
		class NewtonianPolynomial
		{
		private:
			std::vector<double> _x;
			std::vector<double> _c;

			template<class Iterator>
			void calc_coefficients(Iterator first, Iterator last)
			{
				auto y{ std::vector<double>(std::distance(first, last)) };
				size_t index{ 0 };
				for (auto iter = first; iter != last; ++iter)
					y[index++] = *iter;
				for (index = 0; index < y.size(); ++index)
				{
					for (size_t n = 0; n < index; ++n) {
						y[n] /= (_x[n] - _x[index]);
						y[index] /= (_x[index] - _x[n]);
						_c[index] += y[n];
					}
					_c[index] += y[index];
				}
			}

		public:
			/// <summary>
			/// Creating a polynomial via x and y dynamic arrays
			/// </summary>
			/// <param name="x"></param>
			/// <param name="y"></param>
			NewtonianPolynomial(const std::vector<double>& x, std::vector<double>& y) :
				_x{ std::vector<double>(x.size()) }, _c{ std::vector<double>(y.size()) }
			{
				if (x.size() != y.size()) 
					throw std::invalid_argument("X and Y different dimensions!");
				std::memcpy(_x.data(), x.data(), sizeof(double) * x.size());
				calc_coefficients(y.cbegin(), y.cend());
			}
			/// <summary>
			/// Creating a polynomial via x and y static arrays
			/// </summary>
			/// <param name="x"></param>
			/// <param name="y"></param>
			template<size_t n> NewtonianPolynomial(const double(&x)[n], const double(&y)[n]) : 
				_x{ std::vector<double>(n) }, _c{ std::vector<double>(n) }
			{
				std::memcpy(_x.data(), x, sizeof(double) * n);
				calc_coefficients(y, &y[n]);
			}
			/// <summary>
			/// Creating a polynomial via x and y static arrays
			/// </summary>
			/// <param name="x"></param>
			/// <param name="y"></param>
			template<size_t n> NewtonianPolynomial(const std::array<double, n>& x, const std::array<double, n>& y) :
				_x{ std::vector<double>(n) }, _c{ std::vector<double>(n) }
			{
				std::memcpy(_x.data(), x.data(), sizeof(double) * n);
				calc_coefficients(y.cbegin(), y.cend());
			}
			NewtonianPolynomial(const NewtonianPolynomial& p) = default;
			NewtonianPolynomial(NewtonianPolynomial&& p) noexcept : _x{ std::move(p._x) }, _c{ std::move(_c) } {}
			~NewtonianPolynomial() = default;

			NewtonianPolynomial& operator = (const NewtonianPolynomial& p) = default;
			NewtonianPolynomial& operator = (NewtonianPolynomial&& p) noexcept
			{
				_x = std::move(p._x);
				_c = std::move(p._c);
				return *this;
			}
			/// <summary>
			/// Calculating a value
			/// </summary>
			/// <param name="x"></param>
			/// <returns></returns>
			double operator () (const double x) const
			{
				double mult{ 1 }, result{ 0 };
				for (size_t i = 0; i < _x.size(); ++i) {
					result += _c[i] * mult;
					mult *= (x - _x[i]);
				}
				return result;
			}
		};
	}
}