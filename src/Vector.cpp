#include "Vector.h"

namespace general
{
	namespace math 
	{
		VectorDyn& VectorDyn::operator = (VectorDyn&& vec) noexcept
		{
			static_cast<std::vector<double>*>(this)->operator=(vec);
			return *this;
		}
		double VectorDyn::length() const
		{
			double result{ 0 };
			for (size_t i = 0; i < this->size(); ++i) {
				result += this->operator[](i) * this->operator[](i);
			}
			return std::sqrt(result);
		}
		VectorDyn& VectorDyn::operator += (const VectorDyn& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) += vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		VectorDyn& VectorDyn::operator -= (const VectorDyn& vec)
		{
			if (this->size() == vec.size()) {
				for (size_t i = 0; i < std::min(this->size(), vec.size()); ++i)
					this->operator[](i) -= vec[i];
				return *this;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		VectorDyn& VectorDyn::operator *= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i)
				this->operator[](i) *= value;
			return *this;
		}
		VectorDyn& VectorDyn::operator /= (const double value)
		{
			for (size_t i = 0; i < this->size(); ++i)
				this->operator[](i) /= value;
			return *this;
		}

		VectorDyn operator + (const VectorDyn& first, const VectorDyn& second)
		{
			if (first.size() == second.size()) {
				auto vec{ VectorDyn(first.size()) };
				for (size_t i = 0; i < first.size(); ++i)
					vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		VectorDyn operator - (const VectorDyn& first, const VectorDyn& second)
		{
			if (first.size() == second.size()) {
				auto vec{ VectorDyn(first.size()) };
				for (size_t i = 0; i < first.size(); ++i)
					vec[i] = first[i] + second[i];
				return vec;
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		double operator * (const VectorDyn& first, const VectorDyn& second)
		{
			if (first.size() == second.size()) {
				double result{ 0 };
				for (size_t i = 0; i < first.size(); ++i)
					result = first[i] * second[i];
				return std::sqrt(result);
			}
			throw std::invalid_argument("Vectors' dimensions differ!");
		}
		VectorDyn operator * (const double value, const VectorDyn& vec)
		{
			auto result{ VectorDyn(vec) };
			for (size_t i = 0; i < result.size(); ++i)
				result[i] *= value;
			return result;
		}
		VectorDyn operator * (const VectorDyn& vec, const double value)
		{
			auto result{ VectorDyn(vec) };
			for (size_t i = 0; i < result.size(); ++i)
				result[i] *= value;
			return result;
		}
		VectorDyn operator / (const VectorDyn& vec, const double value)
		{
			auto result{ VectorDyn(vec) };
			for (size_t i = 0; i < result.size(); ++i)
				result[i] /= value;
			return result;
		}
		std::ostream& operator <<(std::ostream& os, const VectorDyn& vec)
		{
			os << "{ ";
			for (const auto& v : vec)
				os << v << "; ";
			os << "}";
			return os;
		}
		std::istream& operator >>(std::istream& is, VectorDyn& vec)
		{
			for (auto& v : vec)
				is >> v;
			return is;
		}

		VectorDyn VectorDyn::ones(const size_t size)
		{
			auto vec{ VectorDyn(size) };
			for (size_t i = 0; i < vec.size(); ++i)
				vec[i] = 1.0;
			return vec;
		}

		VectorDyn normalize(const VectorDyn& vec)
		{
			return vec / vec.length();
		}
	}
}