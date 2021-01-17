#include "Legendre.h"
#include <cmath>

namespace general
{
	namespace math
	{
		void LegendrePolynomial::create(
			std::vector<double>& values,
			const size_t degree)
		{
			size_t size = degree + 1;
			values.resize(size);
			std::memset(values.data(), 0, size * sizeof(double));
			switch (degree)
			{
			case 0:
				values[0] = 1.0;
				break;
			case 1:
				values[1] = 1.0;
				break;
			case 2:
				values[0] = -0.5;
				values[2] = 1.5;
				break;
			case 3:
				values[1] = -1.5;
				values[3] = 2.5;
				break;
			case 4:
				values[0] = 3.0 / 8;
				values[2] = -30.0 / 8;
				values[4] = 35.0 / 8;
				break;
			case 5:
				values[1] = 15.0 / 8;
				values[3] = -70.0 / 8;
				values[5] = 63.0 / 8;
				break;
			case 6:
				values[0] = -5.0 / 16;
				values[2] = 105.0 / 16;
				values[4] = -315.0 / 16;
				values[6] = 231.0 / 16;
				break;
			case 7:
				values[1] = -35.0 / 16;
				values[3] = 315.0 / 16;
				values[5] = -693.0 / 16;
				values[7] = 429.0 / 16;
				break;
			case 8:
				values[0] = 35.0 / 128;
				values[2] = -1260.0 / 128;
				values[4] = 6930.0 / 128;
				values[6] = -12012.0 / 128;
				values[8] = 6435.0 / 128;
				break;
			case 9:
				values[1] = 315.0 / 128;
				values[3] = -4620.0 / 128;
				values[5] = 18018.0 / 128;
				values[7] = -25740.0 / 128;
				values[9] = 12155.0 / 128;
				break;
			case 10:
				values[0] = -63.0 / 256;
				values[2] = 3465.0 / 256;
				values[4] = -30030.0 / 256;
				values[6] = 90090.0 / 256;
				values[8] = -109395.0 / 256;
				values[10] = 46189.0 / 256;
				break;
			case 11:
				values[1] = -693.0 / 256;
				values[3] = 15015.0 / 256;
				values[5] = -90090.0 / 256;
				values[7] = 218790.0 / 256;
				values[9] = -230945.0 / 256;
				values[11] = 88179.0 / 256;
				break;
			case 12:
				values[0] = 231.0 / 1024;
				values[2] = -18018.0 / 1024;
				values[4] = 225225.0 / 1024;
				values[6] = -1021020.0 / 1024;
				values[8] = 2078505.0 / 1024;
				values[10] = -1939938.0 / 1024;
				values[12] = 676039.0 / 1024;
				break;
			case 13:
				values[1] = 3003.0 / 1024;
				values[3] = -90090.0 / 1024;
				values[5] = 765765.0 / 1024;
				values[7] = -2771340.0 / 1024;
				values[9] = 4849845.0 / 1024;
				values[11] = -4056234.0 / 1024;
				values[13] = 1300075.0 / 1024;
				break;
			case 14:
				values[0] = -429.0 / 2048;
				values[2] = 45045.0 / 2048;
				values[4] = -765765.0 / 2048;
				values[6] = 4849845.0 / 2048;
				values[8] = -14549535.0 / 2048;
				values[10] = 22309287.0 / 2048;
				values[12] = -16900975.0 / 2048;
				values[14] = 5014575.0 / 2048;
				break;
			case 15:
				values[1] = -6435.0 / 2048;
				values[3] = 255255.0 / 2048;
				values[5] = -2909907.0 / 2048;
				values[7] = 14549535.0 / 2048;
				values[9] = -37182145.0 / 2048;
				values[11] = 50702925.0 / 2048;
				values[13] = -35102025.0 / 2048;
				values[15] = 9694845.0 / 2048;
				break;
			case 16:
				values[0] = 6435.0 / 32768;
				values[2] = -875160.0 / 32768;
				values[4] = 19399380.0 / 32768;
				values[6] = -162954792.0 / 32768;
				values[8] = 669278610.0 / 32768;
				values[10] = -1487285800.0 / 32768;
				values[12] = 1825305300.0 / 32768;
				values[14] = -1163381400.0 / 32768;
				values[16] = 300540195.0 / 32768;
				break;
			case 17:
				values[1] = 109395.0 / 32768;
				values[3] = -5542680.0 / 32768;
				values[5] = 81477396.0 / 32768;
				values[7] = -535422888.0 / 32768;
				values[9] = 1859107250.0 / 32768;
				values[11] = -3650610600.0 / 32768;
				values[13] = 4071834900.0 / 32768;
				values[15] = -2404321560.0 / 32768;
				values[17] = 583401555.0 / 32768;
				break;
			}
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