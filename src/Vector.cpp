#include "Vector.h"

namespace general
{
	namespace math
	{
		Vector& Vector::operator=(Vector&& vector) noexcept
		{
			static_cast<std::vector<double>*>(this)->operator=(vector);
			return *this;
		}
		double Vector::length() const
		{
			double result{ 0.0 };
			for (size_t i = 0; i < size(); ++i)	{
				result += this->operator[](i) * this->operator[](i);
			}
			return std::sqrt(result);
		}
		Vector& Vector::operator+=(const Vector& vector)
		{
			for (size_t i = 0; i < std::min(size(), vector.size()); ++i) {
				this->operator[](i) += vector[i];
			}
			return *this;
		}
		Vector& Vector::operator-=(const Vector& vector)
		{
			for (size_t i = 0; i < std::min(size(), vector.size()); ++i) {
				this->operator[](i) += vector[i];
			}
			return *this;
		}
		Vector& Vector::operator*=(const double value)
		{
			for (size_t i = 0; i < size(); ++i) {
				this->operator[](i) *= value;
			}
			return *this;
		}
		Vector& Vector::operator/=(const double value)
		{
			for (size_t i = 0; i < size(); ++i) {
				this->operator[](i) /= value;
			}
			return *this;
		}
		Vector Vector::ones(const size_t size)
		{
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i) vector[i] = 1.0;
			return vector;
		}
		Vector operator+(const Vector& first, const Vector& second)
		{
			size_t size{ std::min(first.size(), second.size()) };
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i) {
				vector[i] = first[i] + second[i];
			}
			return vector;
		}
		Vector operator-(const Vector& first, const Vector& second)
		{
			size_t size{ std::min(first.size(), second.size()) };
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i) {
				vector[i] = first[i] - second[i];
			}
			return vector;
		}
		double operator*(const Vector& first, const Vector& second)
		{
			double result = 0.0;
			for (size_t i = 0; i < std::min(first.size(), second.size()); ++i) {
				result += first[i] * second[i];
			}
			return result;
		}
		Vector operator*(const double value, const Vector& vector)
		{
			auto result{ Vector(vector.size()) };
			for (size_t i = 0; i < result.size(); ++i) {
				result[i] = vector[i] * value;
			}
			return result;
		}
		Vector operator*(const Vector& vector, const double value)
		{
			auto result{ Vector(vector.size()) };
			for (size_t i = 0; i < result.size(); ++i) {
				result[i] = vector[i] * value;
			}
			return result;
		}
		Vector operator/(const Vector& vector, const double value)
		{
			auto result{ Vector(vector.size()) };
			for (size_t i = 0; i < result.size(); ++i) {
				result[i] = vector[i] / value;
			}
			return result;
		}
		std::ostream& operator<<(std::ostream& os, const Vector& vector)
		{
			os << "{ ";
			for (auto value : vector) {
				os << value << "; ";
			}
			os << "}";
			return os;
		}
		std::istream& operator>>(std::istream& is, Vector& vector)
		{
			for (auto& value : vector) {
				is >> value;
			}
			return is;
		}
	}
}
