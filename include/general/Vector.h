#pragma once
#include <array>
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
		using long_t = long long;

		template<size_t N>
		class Vector : public std::array<double, N>
		{
		public:
			Vector() noexcept : std::array<double, N>() {}
			explicit Vector(const size_t size) : std::array<double, N>(size) {}
			Vector(const std::initializer_list<double>& values) : std::array<double, N>()
			{
				if (values.size() != N) throw std::out_of_range("Invalid input list!");
				size_t i = 0;
				for (const auto& v : values)
					operator[](i++) = v;
			}
			Vector(const double(&values)[N]) : std::array<double, N>() {
				std::memcpy(data(), values, size * sizeof(double));
			}
			Vector(const Vector& vector) = default;
			Vector(Vector&& vector) noexcept : std::array<double, N>(vector) {}
			Vector(const std::array<double, N>& arr) : std::array<double, N>(arr) {}
			template<class IterType> Vector(const IterType& begin, const IterType& end) : std::array<double, N>() 
			{
				size_t n = end - begin;
				if (n != N) 
					throw std::out_of_range("Invalid number of values!");
				n = 0;
				for (auto iter = begin; iter != end; ++iter)
					operator[](n++) = *iter;
			}
			~Vector() noexcept = default;

			Vector& operator = (const Vector& vector) noexcept = default;
			Vector& operator = (Vector&& vector) noexcept {
				data() = std::move(vector.data());
				return *this;
			}

			double length() const
			{
				double value{ 0 };
				for (size_t i = 0; i < N; ++i)
					value += this->operator[](i);
				return std::sqrt(value);
			}

			Vector& operator += (const Vector& vector)
			{
				for (size_t i = 0; i < N; ++i)
					this->operator[](i) += vector[i];
				return *this;
			}
			Vector& operator -= (const Vector& vector) 
			{
				for (size_t i = 0; i < N; ++i)
					operator[](i) -= vector[i];
			}
			Vector& operator *= (const double value) 
			{
				for (size_t i = 0; i < N; ++i)
					operator[](i) *= value;
				return *this;
			}
			Vector& operator /= (const double value) 
			{
				for (size_t i = 0; i < N; ++i)
					operator[](i) /= value;
				return *this;
			}

			friend Vector operator + (const Vector& first, const Vector& second)
			{
				Vector result;
				for (size_t i = 0; i < N; ++i)
					result[i] = first[i] + second[i];
				return result;
			}
			friend Vector operator - (const Vector& first, const Vector& second)
			{
				Vector result;
				for (size_t i = 0; i < N; ++i)
					result[i] = first[i] - second[i];
				return result;
			}
			friend double operator * (const Vector& first, const Vector& second)
			{
				double result{ 0 };
				for (size_t i = 0; i < N; ++i)
					result += first[i] * second[i];
				return result;
			}
			friend Vector operator * (const double value, const Vector& vector)
			{
				Vector result;
				for (size_t i = 0; i < N; ++i)
					result[i] = value * vector[i];
				return result;
			}
			friend Vector operator * (const Vector& vector, const double value)
			{
				Vector result;
				for (size_t i = 0; i < N; ++i)
					result[i] = value * vector[i];
				return result;
			}
			friend Vector operator / (const Vector& vector, const double value)
			{
				Vector result;
				for (size_t i = 0; i < N; ++i)
					result[i] = vector[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const Vector& vector)
			{
				os << "{ ";
				for (size_t i = 0; i < N - 1; ++i)
					os << vector[i] << "; ";
				os << vector[N - 1] << "}";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, Vector& vector)
			{
				for (size_t i = 0; i < N; ++i)
					is >> vector[i];
				return is;
			}

			template<size_t size> static Vector<size> ones()
			{
				Vector<size> vec;
				for (size_t i = 0; i < size; ++i) vec[i] = 1.0;
				return vec;
			}

		};
	}
}