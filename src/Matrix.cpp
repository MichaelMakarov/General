#include "Matrix.h"

namespace general
{
	namespace math
	{
		Matrix::Matrix(const size_t rows, const size_t cols, const std::initializer_list<double>& list)
		{
			if (rows * cols == list.size()) {
				_rows = rows;
				_cols = cols;
				_data = std::vector<double>(list);
			}
			throw std::invalid_argument("Invalid initializer list size!");
		}
		Matrix::Matrix(const size_t rows, const size_t cols, const std::vector<double>& vec)
		{
			if (rows * cols == vec.size()) {
				_rows = rows;
				_cols = cols;
				_data = std::vector<double>(vec);
			}
			throw std::invalid_argument("Invalid vector size!");
		}
		Matrix::Matrix(Matrix&& m) noexcept
		{
			_data = std::move(m._data);
			_rows = std::move(m._rows);
			_cols = std::move(m._cols);
		}

		Matrix& Matrix::operator = (Matrix&& m) noexcept
		{
			_data = std::move(m._data);
			_rows = std::move(m._rows);
			_cols = std::move(m._cols);
			return *this;
		}

		Vector Matrix::get_row(const size_t index) const
		{
			auto row{ Vector(_cols) };
			for (size_t i = 0; i < _cols; ++i)
				row[i] = _data[index * _cols + i];
			return row;
		}
		Vector Matrix::get_column(const size_t index) const
		{
			auto column{ Vector(_rows) };
			for (size_t i = 0; i < _rows; ++i)
				column[i] = _data[i * _cols + index];
			return column;
		}
		void Matrix::set_row(const size_t index, const Vector& row)
		{
			if (_cols == row.size()) {
				for (size_t i = 0; i < _cols; ++i)
					_data[index * _cols + i] = row[i];
			}
			throw std::invalid_argument("Invalid row size!");
		}
		void Matrix::set_column(const size_t index, const Vector& column)
		{
			if (_rows == column.size()) {
				for (size_t i = 0; i < _rows; ++i)
					_data[i * _cols + index] = column[i];
			}
			throw std::invalid_argument("Invalid column size!");
		}

		Matrix& Matrix::operator += (const Matrix& m)
		{
			if (_rows == m._rows && _cols == m._cols) {
				for (size_t i = 0; i < _data.size(); ++i)
					_data[i] += m._data[i];
				return *this;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		Matrix& Matrix::operator -= (const Matrix& m)
		{
			if (_rows == m._rows && _cols == m._cols) {
				for (size_t i = 0; i < _data.size(); ++i)
					_data[i] -= m._data[i];
				return *this;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		Matrix& Matrix::operator *= (const double v)
		{
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] *= v;;
			return *this;
		}
		Matrix& Matrix::operator /= (const double v)
		{
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] /= v;;
			return *this;
		}

		Matrix operator + (const Matrix& f, const Matrix& s)
		{
			if (f._rows == s._rows && f._cols == s._cols) {
				auto m{ Matrix(f._rows, f._cols) };
				for (size_t i = 0; i < f._data.size(); ++i)
					m._data[i] = f._data[i] + s._data[i];
				return m;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		Matrix operator - (const Matrix& f, const Matrix& s)
		{
			if (f._rows == s._rows && f._cols == s._cols) {
				auto m{ Matrix(f._rows, f._cols) };
				for (size_t i = 0; i < f._data.size(); ++i)
					m._data[i] = f._data[i] - s._data[i];
				return m;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		Matrix operator * (const Matrix & f, const Matrix & s)
		{
			if (f._cols == s._rows) {
				auto result{ Matrix(f._rows, s._cols) };
				for (size_t m = 0; m < f._rows; ++m)
					for (size_t k = 0; k < f._cols; ++k)
						for (size_t n = 0; n < s._cols; ++n)
							result(m, n) += f(m, k) * s(k, n);
				return result;
			}
			throw std::invalid_argument("Inconsistent matrices!");
		}
		Vector operator * (const Matrix& matrix, const Vector& vec)
		{
			if (matrix._cols != vec.size())
				throw std::invalid_argument("Incompatible matrix and vector!");
			auto result{ Vector(matrix._rows) };
			for (size_t m = 0; m < matrix._rows; ++m)
				for (size_t n = 0; n < matrix._cols; ++n)
					result[m] += matrix(m, n) * vec[n];
			return result;
		}
		Matrix operator / (const Matrix& m, const double v)
		{
			auto result{ Matrix(m) };
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] /= v;
			return result;
		}
		Matrix operator * (const Matrix& m, const double v)
		{
			auto result{ Matrix(m) };
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] *= v;
			return result;
		}
		Matrix operator * (const double v, const Matrix& m)
		{
			auto result{ Matrix(m) };
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] *= v;
			return result;
		}

		std::ostream& operator << (std::ostream& os, const Matrix& matrix)
		{
			os << "{ ";
			for (size_t m = 0; m < matrix._rows; ++m) {
				os << "{ ";
				for (size_t n = 0; n < matrix._cols; ++n)
					os << matrix(m, n) << "; ";
				os << "} ";
			}
			os << "}";
			return os;
		}
		std::istream& operator << (std::istream& is, Matrix& m) 
		{
			for (size_t i = 0; i < m._data.size(); ++i)
				is >> m._data[i];
			return is;
		}

		static Matrix identity(const size_t size)
		{
			auto m{ Matrix(size, size) };
			for (size_t i = 0; i < size; ++i)
				m(i, i) = 1.0;
			return m;
		}

		Matrix transpose(const Matrix& matrix)
		{
			auto result{ Matrix(matrix._cols, matrix._rows) };
			for (size_t m = 0; m < matrix._rows; ++m)
				for (size_t n = 0; n < matrix._cols; ++n)
					result(n, m) = matrix(m, n);
			return result;
		}
		Matrix inverse(const Matrix& matrix)
		{
			if (matrix._rows != matrix._cols)
				throw std::invalid_argument("Matrix is not square!");
			auto result{ matrix };
			double pivot;
			for (size_t k = 0; k < result._rows; ++k) {
				pivot = result(k, k);
				if (std::abs(pivot) < 1e-20)
					throw std::runtime_error("Degenerate matrix!");
				for (size_t m = 0; m < result._rows; ++m)
					result(m, k) = -result(m, k) / pivot;
				for (size_t m = 0; m < result._rows; ++m) {
					if (m != k) {
						for (size_t n = 0; n < result._cols; ++n) {
							if (n != k)
								result(m, n) += result(k, n) * result(m, k);
						}
					}
				}
				for (size_t n = 0; n < result._cols; ++n)
					result(k, n) /= pivot;
				result(k, k) = 1 / pivot;
			}
			return result;
		}
		void LU(const Matrix& A, Matrix& L, Matrix& U)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			U = L = Matrix(A._rows, A._cols);
			double sum;
			for (size_t m = 0; m < A._rows; ++m) {
				for (size_t n = 0; n < A._cols; ++n) {
					sum = 0;
					if (m <= n) {
						for (size_t i = 0; i < A._cols; ++i)
							sum += L(m, i) * U(i, n);
						U(m, n) = A(m, n) - sum;
					}
					else {
						for (size_t i = 0; i < A._cols; ++i)
							sum += L(m, i) * U(i, n);
						L(m, n) = (A(m, n) - sum) / U(n, n);
					}
				}
				L(m, m) = 1.0;
			}
		}
		void PLU(const Matrix& A, Matrix& L, Matrix& U, std::vector<size_t>& P, size_t& k)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			k = 0;
			L = U = Matrix(A._rows, A._cols);
			P = std::vector<size_t>(A._rows);
			double buf, absv;
			size_t imax;
			for (size_t i = 0; i < P.size(); ++i) P[i] = i;
			for (size_t m = 0; m < A._rows - 1; ++m) {
				buf = std::abs(U(m, m));
				imax = m;
				for (size_t i = m + 1; i < A._rows; ++i)
					if ((absv = std::abs(U(i, m))) > buf) {
						buf = absv;
						imax = i;
					}
				if (buf < 1e-16) throw std::runtime_error("MatrixFix is degenerate!");
				if (imax != m) {
					for (size_t n = m; n < A._cols; ++n) std::swap(U(m, n), U(imax, n));
					for (size_t n = 0; n < m; ++n) std::swap(L(m, n), L(imax, n));
					std::swap(P[m], P[imax]);
					k++;
				}
				for (size_t n = m + 1; n < A._rows; ++n) {
					buf = U(n, m) / U(m, m);
					for (size_t i = m; i < A._cols; ++i) U(n, i) -= buf * U(m, i);
					L(n, m) = buf;
				}
			}
		}
		void QR(const Matrix& A, Matrix& Q, Matrix& R)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			Q = Matrix(A._rows, A._cols);
			auto vecs = std::vector<Vector>(A._cols);
			for (size_t n = 0; n < A._cols; ++n) {
				auto vec = vecs[n] = A.get_column(n);
				// Gramme-Shmidte orthogonalization
				for (size_t k = 0; k < n; ++k)
					vecs[n] -= vecs[k] * (vec * vecs[k]);
				Q.set_column(n, vecs[n] = normalize(vecs[n]));
			}
			R = transpose(Q) * A;
		}
		Vector solve(const Matrix& A, const Vector& B)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			if (A._rows != B.size())
				throw std::invalid_argument("Incompatible matrix and vector!");
			double div{ 1 }, pivot;
			auto D{ Matrix(A._rows, A._cols + 1) };
			for (size_t m = 0; m < A._rows; ++m) {
				for (size_t n = 0; n < A._cols; ++n)
					D(m, n) = A(m, n);
				D(m, A._cols) = B[m];
			}
			for (size_t k = 0; k < D._rows; ++k) {
				pivot = D(k, k);
				if (std::abs(pivot) < 1e-16)
					throw std::runtime_error("Degenerate MatrixFix!");
				for (size_t m = 0; m < D._rows; ++m) {
					if (m != k) {
						for (size_t n = 0; n < D._cols; ++n) {
							if (n != k)
								D(m, n) = div * (D(m, n) * pivot - D(m, k) * D(k, n));
						}
					}
				}
				div = 1 / pivot;
				for (size_t m = 0; m < D._rows; ++m) {
					if (m != k) D(m, k) = 0.0;
				}
			}
			return D.get_column(A._cols) / D(0, 0);
		}
		Matrix AxD(const Matrix& A, const Matrix& D)
		{
			if (A._cols != D._rows)
				throw std::invalid_argument("Incompatible matrices!");
			auto matrix{ A };
			for (size_t m = 0; m < A._rows; ++m)
				for (size_t n = 0; n < A._cols; ++n)
					matrix(m, n) *= D(n, n);
			return matrix;
		}
		Matrix DxA(const Matrix& D, const Matrix& A)
		{
			if (A._cols != D._rows)
				throw std::invalid_argument("Incompatible matrices!");
			auto matrix{ A };
			for (size_t m = 0; m < A._rows; ++m)
				for (size_t n = 0; n < A._cols; ++n)
					matrix(m, n) *= D(m, m);
			return matrix;
		}
		void SVD(const Matrix& A, Matrix& U, Matrix& s, Matrix& V)
		{
			V = Matrix(A._cols, A._cols);
			auto Sinv = s = Matrix(A._cols, A._cols);
			auto M = transpose(A) * A;
			Vector v(M._rows);
			double curr, prev, eps = 1e-16;
			for (size_t i = 0; i < M._rows; ++i) {
				curr = 2.0;
				prev = 1.0;
				for (size_t n = 0; n < M._cols; ++n) v[n] = 1.0;
				while (std::abs((curr - prev) / curr) > eps) {
					prev = curr;
					v = normalize(M * v);
					curr = v[0];
				}
				curr = v * (M * v);
				// filling the MatrixFix of the eigen vectors
				V.set_row(i, v);
				// filling the diagonal MatrixFix
				s(i, i) = std::sqrt(curr);
				for (size_t m = 0; m < M._rows; ++m)
					for (size_t n = 0; n < M._cols; ++n)
						M(m, n) -= curr * v[m] * v[n];
			}
			for (size_t i = 0; i < s._rows; ++i) Sinv(i, i) = 1 / s(i, i);
			U = A * transpose(V) * Sinv;
		}
	}
}