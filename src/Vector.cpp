#include "Vector.h"

namespace general
{
	namespace math
	{
		void Vector::update()
		{
			_begin = _pValues.get();
			_end = &_pValues[_size];
		}
		Vector::Vector(const std::initializer_list<double>& values)
		{
			_size = values.size();
			_pValues = std::make_unique<double[]>(_size);
			size_t index{ 0 };
			for (auto& value : values)
			{
				_pValues[index++] = value;
			}
			update();
		}
		Vector::Vector(const Vector& vector) noexcept
		{
			_size = vector._size;
			_pValues = std::make_unique<double[]>(_size);
			std::memcpy(_pValues.get(), vector._pValues.get(), _size * sizeof(double));
			update();
		}
		Vector::Vector(Vector&& vector) noexcept
		{
			_size = vector._size;
			vector._size = 0;
			_pValues = std::move(vector._pValues);
			update();
		}
		Vector::Vector(const std::vector<double>& vector)
		{
			_size = vector.size();
			_pValues = std::make_unique<double[]>(_size);
			std::memcpy(_pValues.get(), vector.data(), _size * sizeof(double));
			update();
		}
		Vector& Vector::operator=(const Vector& vector) noexcept
		{
			_size = vector._size;
			_pValues = std::make_unique<double[]>(_size);
			std::memcpy(_pValues.get(), vector._pValues.get(), _size * sizeof(double));
			update();
			return *this;
		}
		Vector& Vector::operator=(Vector&& vector) noexcept
		{
			_size = vector._size;
			vector._size = 0;
			_pValues = std::move(vector._pValues);
			update();
			return *this;
		}
		double Vector::length() const
		{
			auto result{ 0.0 };
			for (size_t i = 0; i < _size; ++i)
			{
				result += _pValues[i] * _pValues[i];
			}
			return std::sqrt(result);
		}
		Vector& Vector::operator+=(const Vector& vector)
		{
			for (size_t i = 0; i < std::min(_size, vector._size); ++i)
			{
				_pValues[i] += vector._pValues[i];
			}
			return *this;
		}
		Vector& Vector::operator-=(const Vector& vector)
		{
			for (size_t i = 0; i < std::min(_size, vector._size); ++i)
			{
				_pValues[i] -= vector._pValues[i];
			}
			return *this;
		}
		Vector& Vector::operator*=(const double value)
		{
			for (size_t i = 0; i < _size; ++i)
			{
				_pValues[i] *= value;
			}
			return *this;
		}
		Vector& Vector::operator/=(const double value)
		{
			for (size_t i = 0; i < _size; ++i)
			{
				_pValues[i] /= value;
			}
			return *this;
		}
		Vector Vector::ones(const size_t size)
		{
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i)
			{
				vector._pValues[i] = 1.0;
			}
			return vector;
		}
		Vector operator+(const Vector& first, const Vector& second)
		{
			size_t size{ std::min(first._size, second._size) };
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i)
			{
				vector._pValues[i] = first._pValues[i] + second._pValues[i];
			}
			return vector;
		}
		Vector operator-(const Vector& first, const Vector& second)
		{
			size_t size{ std::min(first._size, second._size) };
			auto vector{ Vector(size) };
			for (size_t i = 0; i < size; ++i)
			{
				vector._pValues[i] = first._pValues[i] - second._pValues[i];
			}
			return vector;
		}
		double operator*(const Vector& first, const Vector& second)
		{
			double result = 0.0;
			for (size_t i = 0; i < std::min(first._size, second._size); ++i)
			{
				result += first._pValues[i] * second._pValues[i];
			}
			return result;
		}
		Vector operator*(const double value, const Vector& vector)
		{
			auto result{ Vector(vector._size) };
			for (size_t i = 0; i < result._size; ++i)
			{
				result._pValues[i] = vector._pValues[i] * value;
			}
			return result;
		}
		Vector operator*(const Vector& vector, const double value)
		{
			auto result{ Vector(vector._size) };
			for (size_t i = 0; i < result._size; ++i)
			{
				result._pValues[i] = vector._pValues[i] * value;
			}
			return result;
		}
		Vector operator/(const Vector& vector, const double value)
		{
			auto result{ Vector(vector._size) };
			for (size_t i = 0; i < result._size; ++i)
			{
				result._pValues[i] = vector._pValues[i] / value;
			}
			return result;
		}
		std::ostream& operator<<(std::ostream& os, const Vector& vector)
		{
			os << "{ ";
			for (auto value : vector)
			{
				os << value << "; ";
			}
			os << "}";
			return os;
		}
		std::istream& operator>>(std::istream& is, Vector& vector)
		{
			for (auto& value : vector)
			{
				is >> value;
			}
			return is;
		}
	}
}
