#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		class Matrix : public IEnumerable<double>
		{
		private:
			std::unique_ptr<double[]> _pValues;
			size_t _nrows, _ncols;

			void update();

		public:
			Matrix() : _nrows{ 0 }, _ncols{ 0 } { update(); }
			Matrix(const size_t rows, const size_t columns) : 
				_nrows{ rows }, _ncols{ columns }, _pValues{ std::make_unique<double[]>(rows * columns) }
			{ update(); }
			Matrix(const size_t rows, const size_t columns, const std::initializer_list<double>& values);
			template<size_t rows, size_t columns> Matrix(const double(&values)[rows][columns])
			{
				_nrows = rows;
				_ncols = columns;
				_pValues = std::make_unique<double[]>(_nrows * _ncols);
				for (size_t m = 0; m < _nrows; ++m)
					std::memcpy(&_pValues[m * _ncols], values[m], _ncols * sizeof(double));
				update();
			}
			Matrix(const size_t rows_number, const size_t columns_number, const std::vector<double> values);
			Matrix(const Matrix& matrix) noexcept;
			Matrix(Matrix&& matrix) noexcept;
			~Matrix() noexcept { _nrows = _ncols = 0; }

			Matrix& operator = (const Matrix& matrix) noexcept;
			Matrix& operator = (Matrix&& matrix) noexcept;

			size_t rows() const { return _nrows; }
			size_t columns() const { return _ncols; }

			double& operator () (const size_t m, const size_t n) const { return _pValues[m * _ncols + n]; }
			double& operator () (long_t m, long_t n) const {
				if (m < 0) m += _nrows;
				if (n < 0) n += _ncols;
				return _pValues[m * _ncols + n];
			}
			Vector get_row(long_t index) const;
			Vector get_column(long_t index) const;
			void set_row(long_t index, const Vector& vector);
			void set_column(long_t index, const Vector& vector);

			Matrix& transpose();

			Matrix& operator += (const Matrix& matrix);
			Matrix& operator -= (const Matrix& matrix);
			Matrix& operator *= (const double value);
			Matrix& operator /= (const double value);

			friend Matrix operator + (const Matrix& first, const Matrix& second);
			friend Matrix operator - (const Matrix& first, const Matrix& second);
			friend Matrix operator * (const Matrix& first, const Matrix& second);
			friend Vector operator * (const Matrix& matrix, const Vector& vector);
			friend Matrix operator * (const Matrix& matrix, const double value);
			friend Matrix operator * (const double value, const Matrix& matrix);
			friend Matrix operator / (const Matrix& matrix, const double value);
			friend std::ostream& operator <<(std::ostream& os, const Matrix& matrix);
			friend std::istream& operator >>(std::istream& is, Matrix& matrix);

			static Matrix identity(const size_t size);
		};
	}
}