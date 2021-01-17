#pragma once
#include <vector>

namespace general
{
	namespace math
	{
		class LegendrePolynomial
		{
		protected:
			size_t _degree;
			bool _normalized;

		private:
			void create(
				std::vector<double>& values,
				const size_t degree);

		public:
			explicit LegendrePolynomial(
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

		class LegendreFunction : public LegendrePolynomial
		{
		private:
			size_t _derivation;

		public:
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
	}
}