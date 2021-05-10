#include "Polynomial.h"
#include "Matrix.h"

namespace general
{
	namespace math
	{
		PowerPolynomial create_polynom(const Vector& x, const Vector& Y, const size_t degree)
		{
			if (x.size() != Y.size())
				throw std::invalid_argument("X and Y have different dimensions!");
			auto A{ Matrix(x.size(), degree + 1) };
			auto polynom{ PowerPolynomial(degree) };
			for (size_t i = 0; i < x.size(); ++i) {
				A(i, 0) = 1.0;
				for (size_t k = 1; k < degree + 1; ++k)
					A(i, k) = A(i, k - 1) * x[i];
			}
			auto AT = transpose(A);
			auto M{ AT * A };
			auto D{ Matrix(M.rows(), M.columns()) };
			for (size_t i = 0; i < M.rows(); ++i) D(i, i) = 1 / std::sqrt(M(i, i));
			auto C = DxA(D, AxD(inverse(DxA(D, AxD(M, D))), D)) * AT * Y;
			for (size_t i = 0; i < degree + 1; ++i)
				polynom[i] = C[i];
			return polynom;
		}

		LegendrePolynomial::LegendrePolynomial(
			const size_t degree,
			const bool normalized)
		{
			_normalized = normalized;
			_degree = degree;
		}

		double LegendrePolynomial::operator () (const double x) const
		{
			double poly1{ 1 };
			double poly2{ x };
			double poly;
			if (_normalized)
			{
				poly2 *= std::sqrt(3);
				poly = _degree > 0 ? poly2 : poly1;
				for (size_t n = 2; n < _degree; ++n)
				{
					poly = std::sqrt(2 * n - 1) * x * poly2 - (n - 1) / std::sqrt(2 * n - 3) * poly1;
					poly *= std::sqrt(2 * n + 1) / n;
					poly1 = poly2;
					poly2 = poly;
				}
			}
			else {
				poly = _degree > 0 ? poly2 : poly1;
				for (size_t n = 2; n < _degree; ++n)
				{
					poly = ((2 * n - 1) * x * poly2 - (n - 1) * poly1) / n;
					poly1 = poly2;
					poly2 = poly;
				}
			}
			return poly;
		}

		LegendreFunction::LegendreFunction(const LegendrePolynomial& p) noexcept
		{
			_degree = p.degree();
			_derivation = 0;
			_normalized = p.normalized();
		}

		LegendreFunction& LegendreFunction::operator=(const LegendrePolynomial& p) noexcept
		{
			_degree = p.degree();
			_derivation = 0;
			_normalized = p.normalized();
			return *this;
		}

		double LegendreFunction::operator () (const double x) const
		{
			if (_derivation > _degree) return 0;
			double funcs[2][2]{ { 1, 0 }, { x, std::sqrt(1 - x * x) } };
			double func;
			const double z{ x / funcs[1][1] };
			const double sqrt3{ std::sqrt(3) };
			if (_normalized)
			{
				funcs[1][0] *= sqrt3;
				funcs[1][1] *= sqrt3;
				func = funcs[std::min(size_t(1), _degree)][std::min(size_t(1), _derivation)];
				for (size_t n = 2; n <= _degree; ++n)
				{
					for (size_t m = 0; m <= std::min(size_t(1), _derivation); ++m)
					{
						func = std::sqrt(2 * n - 1) * x * funcs[1][m] - std::sqrt((n + m - 1) * (n - m - 1) / (2.0 * n - 3)) * funcs[0][m];
						func *= std::sqrt((2 * n + 1.0) / ((n - m) * (n + m)));
						funcs[0][m] = funcs[1][m];
						funcs[1][m] = func;
					}
				}
				for (size_t m = 2; m <= _derivation; ++m)
				{
					func = 2 * (m - 1) * z * funcs[1][1] - std::sqrt((_degree - m + 2) * (_degree + m - 1)) * funcs[1][0];
					func /= std::sqrt((_degree + m) * (_degree - m + 1));
					funcs[1][0] = funcs[1][1];
					funcs[1][1] = func;
				}
			}
			else {
				func = funcs[std::min(size_t(1), _degree)][std::min(size_t(1), _derivation)];
				for (size_t n = 2; n <= _degree; ++n)
				{
					for (size_t m = 0; m <= std::min(size_t(1), _derivation); ++m)
					{
						func = ((2 * n - 1) * x * funcs[1][m] - (n + m - 1) * funcs[0][m]) / (n - m);
						funcs[0][m] = funcs[1][m];
						funcs[1][m] = func;
					}
				}
				for (size_t m = 2; m <= _derivation; ++m)
				{
					func = 2 * (m - 1) * z * funcs[1][1] - (_degree - m + 2) * (_degree + m - 1) * funcs[1][0];
					funcs[1][0] = funcs[1][1];
					funcs[1][1] = func;
				}
			}
			return func;
		}

	}
}