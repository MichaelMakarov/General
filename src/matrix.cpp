#include "matrix.h"

namespace general {
	namespace math	{

		matrix::matrix(const size_t rows, const size_t cols, const std::initializer_list<double>& list)
		{
			if (rows * cols == list.size()) {
				_rows = rows;
				_cols = cols;
				_data = std::vector<double>(list);
			}
			throw std::invalid_argument("Invalid initializer list size!");
		}
		matrix::matrix(const size_t rows, const size_t cols, const std::vector<double>& vec)
		{
			if (rows * cols == vec.size()) {
				_rows = rows;
				_cols = cols;
				_data = std::vector<double>(vec);
			}
			throw std::invalid_argument("Invalid vector size!");
		}
		matrix::matrix(matrix&& m) noexcept
		{
			_data = std::move(m._data);
			_rows = std::move(m._rows);
			_cols = std::move(m._cols);
		}

		matrix& matrix::operator = (matrix&& m) noexcept
		{
			_data = std::move(m._data);
			_rows = std::move(m._rows);
			_cols = std::move(m._cols);
			return *this;
		}

		vector matrix::get_row(const size_t index) const
		{
			auto row{ vector(_cols) };
			for (size_t i = 0; i < _cols; ++i)
				row[i] = _data[index * _cols + i];
			return row;
		}
		vector matrix::get_column(const size_t index) const
		{
			auto column{ vector(_rows) };
			for (size_t i = 0; i < _rows; ++i)
				column[i] = _data[i * _cols + index];
			return column;
		}
		void matrix::set_row(const size_t index, const vector& row)
		{
			if (_cols == row.size()) {
				for (size_t i = 0; i < _cols; ++i)
					_data[index * _cols + i] = row[i];
			}
			throw std::invalid_argument("Invalid row size!");
		}
		void matrix::set_column(const size_t index, const vector& column)
		{
			if (_rows == column.size()) {
				for (size_t i = 0; i < _rows; ++i)
					_data[i * _cols + index] = column[i];
			}
			throw std::invalid_argument("Invalid column size!");
		}

		matrix& matrix::operator += (const matrix& m)
		{
			if (_rows == m._rows && _cols == m._cols) {
				for (size_t i = 0; i < _data.size(); ++i)
					_data[i] += m._data[i];
				return *this;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		matrix& matrix::operator -= (const matrix& m)
		{
			if (_rows == m._rows && _cols == m._cols) {
				for (size_t i = 0; i < _data.size(); ++i)
					_data[i] -= m._data[i];
				return *this;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		matrix& matrix::operator *= (const double v)
		{
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] *= v;;
			return *this;
		}
		matrix& matrix::operator /= (const double v)
		{
			for (size_t i = 0; i < _data.size(); ++i)
				_data[i] /= v;;
			return *this;
		}

		matrix operator + (const matrix& f, const matrix& s)
		{
			if (f._rows == s._rows && f._cols == s._cols) {
				auto m{ matrix(f._rows, f._cols) };
				for (size_t i = 0; i < f._data.size(); ++i)
					m._data[i] = f._data[i] + s._data[i];
				return m;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		matrix operator - (const matrix& f, const matrix& s)
		{
			if (f._rows == s._rows && f._cols == s._cols) {
				auto m{ matrix(f._rows, f._cols) };
				for (size_t i = 0; i < f._data.size(); ++i)
					m._data[i] = f._data[i] - s._data[i];
				return m;
			}
			throw std::invalid_argument("Matrices' dimensions differ!");
		}
		matrix operator * (const matrix & f, const matrix & s)
		{
			if (f._cols == s._rows) {
				auto result{ matrix(f._rows, s._cols) };
				for (size_t m = 0; m < f._rows; ++m)
					for (size_t k = 0; k < f._cols; ++k)
						for (size_t n = 0; n < s._cols; ++n)
							result(m, n) += f(m, k) * s(k, n);
				return result;
			}
			throw std::invalid_argument("Inconsistent matrices!");
		}
		vector operator * (const matrix& matrix, const vector& vec)
		{
			if (matrix._cols != vec.size())
				throw std::invalid_argument("Incompatible matrix and vector!");
			auto result{ vector(matrix._rows) };
			for (size_t m = 0; m < matrix._rows; ++m)
				for (size_t n = 0; n < matrix._cols; ++n)
					result[m] += matrix(m, n) * vec[n];
			return result;
		}
		matrix operator / (const matrix& m, const double v)
		{
			auto result = matrix(m);
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] /= v;
			return result;
		}
		matrix operator * (const matrix& m, const double v)
		{
			auto result{ matrix(m) };
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] *= v;
			return result;
		}
		matrix operator * (const double v, const matrix& m)
		{
			auto result{ matrix(m) };
			for (size_t i = 0; i < result._data.size(); ++i)
				result._data[i] *= v;
			return result;
		}

		std::ostream& operator << (std::ostream& os, const matrix& matrix)
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
		std::istream& operator << (std::istream& is, matrix& m) 
		{
			for (size_t i = 0; i < m._data.size(); ++i)
				is >> m._data[i];
			return is;
		}

		static matrix identity(const size_t dim)
		{
			auto m{ matrix(dim, dim) };
			for (size_t i = 0; i < dim; ++i)
				m(i, i) = 1.0;
			return m;
		}

		matrix transpose(const matrix& init)
		{
			auto result{ matrix(init._cols, init._rows) };
			for (size_t m = 0; m < init._rows; ++m)
				for (size_t n = 0; n < init._cols; ++n)
					result(n, m) = init(m, n);
			return result;
		}
		matrix inverse(const matrix& matrix)
		{
			if (matrix._rows != matrix._cols)
				throw std::invalid_argument("Matrix is not square!");
			auto result{ matrix };
			double pivot;
			for (size_t k = 0; k < result._rows; ++k) {
				pivot = result(k, k);
				if (std::abs(pivot) < 1e-20)
					throw std::runtime_error("Degenerate matrix!");
				for (size_t m = 0; m < result._rows; ++m) result(m, k) /= -pivot;
				for (size_t m = 0; m < result._rows; ++m) {
					if (m != k) {
						for (size_t n = 0; n < result._cols; ++n) {
							if (n != k)
								result(m, n) += result(k, n) * result(m, k);
						}
					}
				}
				for (size_t n = 0; n < result._cols; ++n) result(k, n) /= pivot;
				result(k, k) = 1 / pivot;
			}
			return result;
		}
		void LU(const matrix& A, matrix& L, matrix& U)
		{
			if (A._rows != A._cols) 
				throw std::invalid_argument("Matrix is not square!");
			U = L = matrix(A._rows, A._cols);
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
		void PLU(const matrix& A, matrix& L, matrix& U, std::vector<size_t>& P, size_t& k)
		{
			if (A._rows != A._cols) 
				throw std::invalid_argument("Matrix is not square!");
			k = 0;
			L = U = matrix(A._rows, A._cols);
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
				if (buf < 1e-16) 
					throw std::runtime_error("Matrix is degenerate!");
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
		void QR(const matrix& A, matrix& Q, matrix& R)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			Q = matrix(A._rows, A._cols);
			auto vecs = std::vector<vector>(A._cols);
			for (size_t n = 0; n < A._cols; ++n) {
				auto vec = vecs[n] = A.get_column(n);
				// Gramme-Shmidte orthogonalization
				for (size_t k = 0; k < n; ++k)
					vecs[n] -= vecs[k] * (vec * vecs[k]);
				Q.set_column(n, vecs[n] = normalize(vecs[n]));
			}
			R = transpose(Q) * A;
		}
		vector solve(const matrix& A, const vector& B)
		{
			if (A._rows != A._cols)
				throw std::invalid_argument("Matrix is not square!");
			if (A._rows != B.size())
				throw std::invalid_argument("Incompatible matrix and vector!");
			double div{ 1 }, pivot;
			auto D{ matrix(A._rows, A._cols + 1) };
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
		matrix AxD(const matrix& A, const matrix& D)
		{
			if (A._cols != D._rows)
				throw std::invalid_argument("Incompatible matrices!");
			auto matrix{ A };
			for (size_t m = 0; m < A._rows; ++m)
				for (size_t n = 0; n < A._cols; ++n)
					matrix(m, n) *= D(n, n);
			return matrix;
		}
		matrix DxA(const matrix& D, const matrix& A)
		{
			if (A._cols != D._rows)
				throw std::invalid_argument("Incompatible matrices!");
			auto matrix{ A };
			for (size_t m = 0; m < A._rows; ++m)
				for (size_t n = 0; n < A._cols; ++n)
					matrix(m, n) *= D(m, m);
			return matrix;
		}
		void SVD(const matrix& A, matrix& U, matrix& s, matrix& V)
		{
			V = matrix(A._cols, A._cols);
			auto Sinv = s = matrix(A._cols, A._cols);
			auto M = transpose(A) * A;
			vector v(M._rows);
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