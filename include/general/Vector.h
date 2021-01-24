#pragma once
#include <memory>
#include <vector>
#include <ostream>
#include <istream>
#include "IEnumerable.h"

namespace general
{
	namespace math
	{
		using namespace general;
		using long_t = long long;

		class Vector : public IEnumerable<double>
		{
		private:
			std::unique_ptr<double[]> _pValues;
			size_t _size;

			void update();

		public:
			Vector() : _size{ 0 } { update(); }
			explicit Vector(const size_t size) : _size{ size }, _pValues{ std::make_unique<double[]>(size) } { update(); }
			Vector(const std::initializer_list<double>& values);
			template<size_t size> Vector(const double(&values)[size])
			{
				_size = size;
				_pValues = std::make_unique<double[]>(_size);
				std::memcpy(_pValues.get(), values, _size * sizeof(double));
			}
			Vector(const Vector& vector) noexcept;
			Vector(Vector&& vector) noexcept;
			Vector(const std::vector<double>& vector);
			template<class IterType> Vector(IterType& begin, IterType& end);
			~Vector() noexcept { _size = 0; _begin = _end = nullptr;  }

			Vector& operator = (const Vector& vector) noexcept;
			Vector& operator = (Vector&& vector) noexcept;

			size_t size() const { return _size; }
			double length() const;

			double& operator [] (const size_t index) const { return _pValues[index]; }
			double& operator [] (const long_t index) const { return _pValues[index < 0 ? index + _size : index]; }

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
		template<class IterType>
		inline Vector::Vector(IterType& begin, IterType& end)
		{
			_size = end - begin;
			_pValues = std::make_unique<double[]>(_size);
			size_t index{ 0 };
			for (auto& iter = begin; iter != end; iter++)
			{
				_pValues[index++] = *iter;
			}
			update();
		}
	}
}