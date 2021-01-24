#include "Matrix.h"

namespace general
{
	namespace math
	{
		void Matrix::update()
		{
			_begin = _pValues.get();
			_end = &_pValues[_nrows * _ncols];
		}
		Matrix::Matrix(const size_t rows, const size_t columns, const std::initializer_list<double>& values)
		{
			if (values.size() != rows * columns)
				throw std::invalid_argument("Invalid dimensions!");
			_nrows = rows;
			_ncols = columns;
			_pValues = std::make_unique<double[]>(_nrows * _ncols);
			size_t index{ 0 };
			for (const auto& value : values)
				_pValues[index++] = value;
			update();
		}

		Matrix::Matrix(const size_t rows_number, const size_t columns_number, const std::vector<double> values)
		{
			size_t size{ rows_number * columns_number };
			if (values.size() != size)
				throw std::invalid_argument("Invalid dimensions!");
			_nrows = rows_number;
			_ncols = columns_number;
			_pValues = std::make_unique<double[]>(size);
			std::memcpy(_pValues.get(), values.data(), size * sizeof(double));
			update();
		}

		Matrix::Matrix(const Matrix& matrix) noexcept
		{
			_nrows = matrix._nrows;
			_ncols = matrix._ncols;
			size_t size{ _nrows * _ncols };
			_pValues = std::make_unique<double[]>(size);
			std::memcpy(_pValues.get(), matrix._pValues.get(), size * sizeof(double));
			update();
		}

		Matrix::Matrix(Matrix&& matrix) noexcept
		{
			_nrows = matrix._nrows;
			_ncols = matrix._ncols;
			_pValues = std::move(matrix._pValues);
			update();
			matrix._nrows = matrix._ncols = 0;
		}

		Matrix& Matrix::operator=(const Matrix& matrix) noexcept
		{
			_nrows = matrix._nrows;
			_ncols = matrix._ncols;
			size_t size{ _nrows * _ncols };
			_pValues = std::make_unique<double[]>(size);
			std::memcpy(_pValues.get(), matrix._pValues.get(), size * sizeof(double));
			update();
			return *this;
		}

		Matrix& Matrix::operator=(Matrix&& matrix) noexcept
		{
			_nrows = matrix._nrows;
			_ncols = matrix._ncols;
			_pValues = std::move(matrix._pValues);
			update();
			matrix._nrows = matrix._ncols = 0;
			return *this;
		}

		Vector Matrix::get_row(long_t index) const
		{
			if (index < 0) index += _nrows;
			auto vector{ Vector(_ncols) };
			for (size_t i = 0; i < _ncols; ++i)
				vector[i] = _pValues[index * _ncols + i];
			return vector;
		}

		Vector Matrix::get_column(long_t index) const
		{
			if (index < 0) index += _ncols;
			auto vector{ Vector(_nrows) };
			for (size_t i = 0; i < _nrows; ++i)
				vector[i] = _pValues[i * _ncols + index];
			return vector;
		}

		void Matrix::set_row(long_t index, const Vector& vector)
		{
			if (vector.size() == _ncols)
			{
				if (index < 0) index += _nrows;
				for (size_t i = 0; i < _ncols; ++i)
					_pValues[index * _ncols + i] = vector[i];
			}
			else {
				throw std::runtime_error("Inconsistent dimensions!");
			}
		}

		void Matrix::set_column(long_t index, const Vector& vector)
		{
			if (vector.size() == _nrows)
			{
				if (index < 0) index += _ncols;
				for (size_t i = 0; i < _nrows; ++i)
					_pValues[i * _ncols + index] = vector[i];
			}
			else {
				throw std::runtime_error("Inconsistent dimensions!");
			}
		}

		Matrix& Matrix::transpose()
		{
			auto pValues{ std::make_unique<double[]>(_nrows * _ncols) };
			for (size_t m = 0; m < _nrows; ++m)
			{
				for (size_t n = 0; n < _ncols; ++n)
					pValues[m * _ncols + n] = _pValues[n * _nrows + m];
			}
			_pValues = std::move(pValues);
			std::swap(_nrows, _ncols);
			return *this;
		}

		Matrix& Matrix::operator+=(const Matrix& matrix)
		{
			if (_nrows == matrix._nrows && _ncols == matrix._ncols)
			{
				size_t index;
				for (size_t m = 0; m < _nrows; ++m)
				{
					for (size_t n = 0; n < _ncols; ++n)
					{
						index = m * _ncols + n;
						_pValues[index] += matrix._pValues[index];
					}
				}
				return *this;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Matrix& Matrix::operator-=(const Matrix& matrix)
		{
			if (_nrows == matrix._nrows && _ncols == matrix._ncols)
			{
				size_t index;
				for (size_t m = 0; m < _nrows; ++m)
				{
					for (size_t n = 0; n < _ncols; ++n)
					{
						index = m * _ncols + n;
						_pValues[index] -= matrix._pValues[index];
					}
				}
				return *this;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Matrix& Matrix::operator*=(const double value)
		{
			for (size_t m = 0; m < _nrows; ++m)
			{
				for (size_t n = 0; n < _ncols; ++n)
					_pValues[m * _ncols + n] *= value;
			}
			return *this;
		}

		Matrix& Matrix::operator/=(const double value)
		{
			for (size_t m = 0; m < _nrows; ++m)
			{
				for (size_t n = 0; n < _ncols; ++n)
					_pValues[m * _ncols + n] *= value;
			}
			return *this;
		}

		Matrix Matrix::identity(const size_t size)
		{
			auto matrix{ Matrix(size, size) };
			for (size_t i = 0; i < size; ++i)
				matrix._pValues[i * size + i] = 1.0;
			return matrix;
		}

		Matrix operator+(const Matrix& first, const Matrix& second)
		{
			if (first._nrows == second._nrows && first._ncols == second._ncols)
			{
				auto matrix{ Matrix(first._nrows, first._ncols) };
				size_t index;
				for (size_t m = 0; m < first._nrows; ++m)
				{
					for (size_t n = 0; n < first._ncols; ++n)
					{
						index = m * first._ncols + n;
						matrix._pValues[index] = first._pValues[index] + second._pValues[index];
					}
				}
				return matrix;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Matrix operator-(const Matrix& first, const Matrix& second)
		{
			if (first._nrows == second._nrows && first._ncols == second._ncols)
			{
				auto matrix{ Matrix(first._nrows, first._ncols) };
				size_t index;
				for (size_t m = 0; m < first._nrows; ++m)
				{
					for (size_t n = 0; n < first._ncols; ++n)
					{
						index = m * first._ncols + n;
						matrix._pValues[index] = first._pValues[index] - second._pValues[index];
					}
				}
				return matrix;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Matrix operator*(const Matrix& first, const Matrix& second)
		{
			if (first._ncols == second._nrows)
			{
				auto matrix{ Matrix(first._nrows, second._ncols) };
				for (size_t m = 0; m < first._nrows; ++m)
				{
					for (size_t k = 0; k < second._nrows; ++k)
					{
						for (size_t n = 0; n < second._ncols; ++n)
							matrix._pValues[m * first._ncols + n] = 
								first._pValues[m * first._ncols + k] * second._pValues[k * second._ncols + n];
					}
				}
				return matrix;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Vector operator*(const Matrix& matrix, const Vector& vector)
		{
			if (matrix._ncols == vector.size())
			{
				auto result{ Vector(matrix._nrows) };
				for (size_t m = 0; m < matrix._nrows; ++m)
				{
					for (size_t n = 0; n < matrix._ncols; ++n)
						result[m] += matrix._pValues[m * matrix._ncols + n] * vector[n];
				}
				return result;
			}
			throw std::runtime_error("Inconsistent dimensions!");
		}

		Matrix operator*(const Matrix& matrix, const double value)
		{
			auto result{ Matrix(matrix._nrows, matrix._ncols) };
			size_t index;
			for (size_t m = 0; m < matrix._nrows; ++m)
			{
				for (size_t n = 0; n < matrix._ncols; ++n)
				{
					index = m * matrix._ncols + n;
					result._pValues[index] = matrix._pValues[index] * value;
				}
			}
			return result;
		}

		Matrix operator*(const double value, const Matrix& matrix)
		{
			auto result{ Matrix(matrix._nrows, matrix._ncols) };
			size_t index;
			for (size_t m = 0; m < matrix._nrows; ++m)
			{
				for (size_t n = 0; n < matrix._ncols; ++n)
				{
					index = m * matrix._ncols + n;
					result._pValues[index] = matrix._pValues[index] * value;
				}
			}
			return result;
		}

		Matrix operator/(const Matrix& matrix, const double value)
		{
			auto result{ Matrix(matrix._nrows, matrix._ncols) };
			size_t index;
			for (size_t m = 0; m < matrix._nrows; ++m)
			{
				for (size_t n = 0; n < matrix._ncols; ++n)
				{
					index = m * matrix._ncols + n;
					result._pValues[index] = matrix._pValues[index] / value;
				}
			}
			return result;
		}

		std::ostream& operator<<(std::ostream& os, const Matrix& matrix)
		{
			os << "{ ";
			if (matrix._pValues)
			{
				for (size_t m = 0; m < matrix._nrows; ++m)
				{
					os << "{ ";
					for (size_t n = 0; n < matrix._ncols; ++n)
						os << matrix._pValues[m * matrix._ncols + n] << "; ";
					os << "} ";
				}
			}
			os << "}";
			return os;
		}

		std::istream& operator>>(std::istream& is, Matrix& matrix)
		{
			for (size_t i = 0; i < matrix._nrows * matrix._ncols; ++i)
				is >> matrix._pValues[i];
			return is;
		}

}
}