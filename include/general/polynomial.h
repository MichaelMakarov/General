#pragma once
#include "matrix.h"

namespace general {
	namespace math {

		/// <summary>
		/// Class represents a polynom as a single variable function
		/// </summary>
		template <size_t N> class power_polynomial {
		public:
			/// <summary>
			/// coefficients of the polynomial
			/// </summary>
			double data[N + 1]{};

		public:
			constexpr size_t degree() const noexcept { return N; }
			constexpr double& operator [] (const size_t index) noexcept { return data[index]; }
			constexpr const double& operator [] (const size_t index) const noexcept { return data[index]; }

			constexpr double operator () (const double x) const noexcept
			{
				double result{ data[0] };
				double mult{ x };
				for (size_t i = 1; i < N + 1; ++i) {
					result += mult * data[i];
					mult *= x;
				}
				return result;
			}
		};

		/// <summary>
		/// least squares polynomial approximation solving task sum||y - p(x)||^2 -> min
		/// </summary>
		/// <typeparam name="Container">a data container type</typeparam>
		/// <param name="x">variable values</param>
		/// <param name="y">values of the function</param>
		/// <returns>approximating polynomial</returns>
		template<size_t degree, class Container> 
		power_polynomial<degree> lstsq(const Container& x, const Container& y) {
			auto A{ matrix(x.size(), degree + 1) };
			auto B{ vector(y) };
			power_polynomial<degree> poly;
			size_t i{ 0 };
			for (const auto& elem : x) {
				A(i, 0) = 1.0;
				for (size_t k = 1; k < degree + 1; ++k)
					A(i, k) = A(i, k - 1) * elem;
				++i;
			}
			i = 0;
			for (const auto& elem : y) B[i++] = elem;
			const auto AT = transpose(A);
			const auto M{ AT * A };
			auto D{ matrix(M.rows(), M.columns()) };
			for (i = 0; i < M.rows(); ++i) D(i, i) = 1 / std::sqrt(M(i, i));
			auto c = DxA(D, AxD(inverse(DxA(D, AxD(M, D))), D)) * AT * B;
			for (i = 0; i < degree + 1; ++i) poly[i] = c[i];
			return poly;
		}
		

		/// <summary>
		/// least squares polynomial approximation solving task sum||y - p(x)||^2 -> min
		/// </summary>
		/// <param name="x">variable values</param>
		/// <param name="y">values of the function</param>
		/// <returns>approximating polynomial</returns>
		template<size_t degree, size_t size> power_polynomial<degree> lstsq(const double(&x)[size], const double(&y)[size]) {
			vec<size> B;
			power_polynomial<degree> poly;
			matrix_mxn<size, degree + 1> A;
			for (size_t i{ 0 }; i < size; ++i) {
				A(i, 0) = 1.0;
				for (size_t k{ 0 }; k < degree; ++k)
					A(i, k + 1) = A(i, k) * x[i];
				B[i] = y[i];
			} 
			const auto At = transpose(A);
			const auto M{ At * A };
			matrix_mxn<degree + 1, degree + 1> D;
			for (size_t i{ 0 }; i < M.rows(); ++i) D(i, i) = 1 / std::sqrt(M(i, i));
			auto C = DxA(D, AxD(inverse(DxA(D, AxD(M, D))), D)) * At * B;
			auto c = inverse(At * A) * At * B;
			for (size_t i{ 0 }; i < degree + 1; ++i) poly[i] = c[i];
			return poly;
		}

		/// <summary>
		/// Class represents Legendre polynomial Pn
		/// </summary>
		class legendre_polynomial {
		protected:
			size_t _degree;
			bool _normalized;

		public:
			/// <summary>
			/// Creating a Legendre polynomial
			/// </summary>
			/// <param name="degree"> - a degree</param>
			/// <param name="normalized"> - a flag whether polynomial is normalized or not</param>
			legendre_polynomial(
				const size_t degree = 0,
				const bool normalized = false);
			legendre_polynomial(const legendre_polynomial& p) noexcept = default;
			legendre_polynomial(legendre_polynomial&& p) noexcept = default;
			~legendre_polynomial() noexcept = default;

			legendre_polynomial& operator = (const legendre_polynomial& p) noexcept = default;
			legendre_polynomial& operator = (legendre_polynomial&& p) noexcept = default;

			size_t degree() const { return _degree; }

			bool normalized() const { return _normalized; }

			virtual double operator () (const double x) const;
		};
		/// <summary>
		/// Class represents a Legendre function Pnm
		/// </summary>
		class legendre_function : public legendre_polynomial
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
			legendre_function(
				const size_t degree = 0,
				const size_t derivation = 0,
				const bool normalized = false) :
				legendre_polynomial(degree, normalized),
				_derivation{ derivation }
			{}
			legendre_function(const legendre_polynomial& p) noexcept;
			legendre_function(const legendre_function& f) noexcept = default;
			legendre_function(legendre_function&& f) noexcept = default;
			~legendre_function() noexcept = default;

			legendre_function& operator = (const legendre_polynomial& p) noexcept;
			legendre_function& operator = (const legendre_function& f) noexcept = default;
			legendre_function& operator = (legendre_function&& f) noexcept = default;

			size_t derivation() const { return _derivation; }

			bool normalized() const { return _normalized; }

			double operator () (const double x) const;
		};
		/// <summary>
		/// Class represents Newtonian polynomial
		/// </summary>
		class newtonian_polynomial
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
			newtonian_polynomial(const std::vector<double>& x, std::vector<double>& y) :
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
			template<size_t n> newtonian_polynomial(const double(&x)[n], const double(&y)[n]) : 
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
			template<size_t n> newtonian_polynomial(const std::array<double, n>& x, const std::array<double, n>& y) :
				_x{ std::vector<double>(n) }, _c{ std::vector<double>(n) }
			{
				std::memcpy(_x.data(), x.data(), sizeof(double) * n);
				calc_coefficients(y.cbegin(), y.cend());
			}
			newtonian_polynomial(const newtonian_polynomial& p) = default;
			newtonian_polynomial(newtonian_polynomial&& p) noexcept : _x{ std::move(p._x) }, _c{ std::move(_c) } {}
			~newtonian_polynomial() = default;

			newtonian_polynomial& operator = (const newtonian_polynomial& p) = default;
			newtonian_polynomial& operator = (newtonian_polynomial&& p) noexcept
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