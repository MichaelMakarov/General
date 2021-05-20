#include "Vector.h"

namespace general
{
	namespace math 
	{
		Vec3 cross(const Vec3& f, const Vec3& s) noexcept
		{
			return Vec3({
				f[1] * s[2] - f[2] * s[1],
				f[2] * s[0] - f[0] * s[2],
				f[0] * s[1] - f[1] * s[0]
			});
		}

		Vector& Vector::operator = (Vector&& vec) noexcept
		{
			static_cast<std::vector<double>*>(this)->operator=(vec);
			return *this;
		}
		double Vector::length() const
		{
			double result{ 0 };
			for (size_t i = 0; i < this->size(); ++i) {
				result += this->operator[](i) * this->operator[](i);
			}
			return std::sqrt(result);
		}
		Vector& Vector::operator += (const Vector& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) += vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		Vector& Vector::operator -= (const Vector& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) -= vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		Vector& Vector::operator *= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i) this->operator[](i) *= value;
			return *this;
		}
		Vector& Vector::operator /= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i) this->operator[](i) /= value;
			return *this;
		}

		Vector operator + (const Vector& first, const Vector& second)
		{
			if (first.size() == second.size()) {
				auto vec{ Vector(first.size()) };
				for (size_t i = 0; i < first.size(); ++i) vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		Vector operator - (const Vector& first, const Vector& second)
		{
			if (first.size() == second.size()) {
				auto vec{ Vector(first.size()) };
				for (size_t i = 0; i < first.size(); ++i) vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		double operator * (const Vector& first, const Vector& second)
		{
			if (first.size() == second.size()) {
				double result{ 0 };
				for (size_t i = 0; i < first.size(); ++i) result += first[i] * second[i];
				return result;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		Vector operator * (const double value, const Vector& vec)
		{
			auto result{ Vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] *= value;
			return result;
		}
		Vector operator * (const Vector& vec, const double value)
		{
			auto result{ Vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] *= value;
			return result;
		}
		Vector operator / (const Vector& vec, const double value)
		{
			auto result{ Vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] /= value;
			return result;
		}
		std::ostream& operator <<(std::ostream& os, const Vector& vec)
		{
			os << "{ ";
			for (const auto& v : vec) os << v << "; ";
			os << "}";
			return os;
		}
		std::istream& operator >>(std::istream& is, Vector& vec)
		{
			for (auto& v : vec) is >> v;
			return is;
		}

		Vector Vector::ones(const size_t size)
		{
			auto vec{ Vector(size) };
			for (size_t i = 0; i < vec.size(); ++i)	vec[i] = 1.0;
			return vec;
		}

		Vector normalize(const Vector& vec)
		{
			return vec / vec.length();
		}
	}
}