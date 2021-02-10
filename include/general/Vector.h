#pragma once
#include <memory>
#include <vector>
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
		using namespace general;
		using long_t = long long;

		class Vector : public std::vector<double>
		{
		public:
			Vector() : std::vector<double>() {}
			explicit Vector(const size_t size) : std::vector<double>(size) {}
			Vector(const std::initializer_list<double>& values) : std::vector<double>(values) {}
			template<size_t size> Vector(const double(&values)[size]) : std::vector<double>(size)
			{
				std::memcpy(data(), values, size * sizeof(double));
			}
			Vector(const Vector& vector) = default;
			Vector(Vector&& vector) noexcept : std::vector<double>(vector) {}
			Vector(const std::vector<double>& vector) : std::vector<double>(vector) {}
			template<class IterType> Vector(IterType& begin, IterType& end) : std::vector<double>(begin, end) {}
			~Vector() noexcept = default;

			Vector& operator = (const Vector& vector) noexcept = default;
			Vector& operator = (Vector&& vector) noexcept;

			double length() const;

			const double& operator [] (const long_t index) const
			{ 
				return static_cast<const std::vector<double>*>(this)->operator[](index < 0 ? index + size() : index); 
			}
			double& operator [] (const long_t index)
			{
				return static_cast<std::vector<double>*>(this)->operator[](index < 0 ? index + size() : index);
			}

			Vector& operator += (const Vector& vector);
			Vector& operator -= (const Vector& vector);
			Vector& operator *= (const double value);
			Vector& operator /= (const double value);

			friend Vector operator + (const Vector& first, const Vector& second);
			friend Vector operator - (const Vector& first, const Vector& second);
			friend double operator * (const Vector& first, const Vector& second);
			friend Vector operator * (const double value, const Vector& vector);
			friend Vector operator * (const Vector& vector, const double value);
			friend Vector operator / (const Vector& vector, const double value);
			friend std::ostream& operator <<(std::ostream& os, const Vector& vector);
			friend std::istream& operator >>(std::istream& is, Vector& vector);

			static Vector ones(const size_t size);

		};
	}
}