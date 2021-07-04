#include "vector.h"

namespace general {
	namespace math {
		vector& vector::operator = (vector&& vec) noexcept
		{
			static_cast<std::vector<double>*>(this)->operator=(vec);
			return *this;
		}
		double vector::length() const
		{
			double result{ 0 };
			for (size_t i = 0; i < this->size(); ++i) {
				result += this->operator[](i) * this->operator[](i);
			}
			return std::sqrt(result);
		}
		vector& vector::operator += (const vector& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) += vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		vector& vector::operator -= (const vector& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) -= vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		vector& vector::operator *= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i) this->operator[](i) *= value;
			return *this;
		}
		vector& vector::operator /= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i) this->operator[](i) /= value;
			return *this;
		}

		vector operator + (const vector& first, const vector& second)
		{
			if (first.size() == second.size()) {
				auto vec{ vector(first.size()) };
				for (size_t i = 0; i < first.size(); ++i) vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		vector operator - (const vector& first, const vector& second)
		{
			if (first.size() == second.size()) {
				auto vec{ vector(first.size()) };
				for (size_t i = 0; i < first.size(); ++i) vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		double operator * (const vector& first, const vector& second)
		{
			if (first.size() == second.size()) {
				double result{ 0 };
				for (size_t i = 0; i < first.size(); ++i) result += first[i] * second[i];
				return result;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		vector operator * (const double value, const vector& vec)
		{
			auto result{ vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] *= value;
			return result;
		}
		vector operator * (const vector& vec, const double value)
		{
			auto result{ vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] *= value;
			return result;
		}
		vector operator / (const vector& vec, const double value)
		{
			auto result{ vector(vec) };
			for (size_t i = 0; i < result.size(); ++i) result[i] /= value;
			return result;
		}
		std::ostream& operator <<(std::ostream& os, const vector& vec)
		{
			os << "{ ";
			for (const auto& v : vec) os << v << "; ";
			os << "}";
			return os;
		}
		std::istream& operator >>(std::istream& is, vector& vec)
		{
			for (auto& v : vec) is >> v;
			return is;
		}

		vector vector::ones(const size_t dim)
		{
			auto vec{ vector(dim) };
			for (size_t i = 0; i < vec.size(); ++i)	vec[i] = 1.0;
			return vec;
		}

		vector normalize(const vector& vec)
		{
			return vec / vec.length();
		}
	}
}