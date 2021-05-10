#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		template<size_t M, size_t N> 
		class MatrixMxN : public std::array<double, M * N>
		{
			static_assert(M * N > size_t(0), "MatrixFix size is equal zero!");
		public:
			MatrixMxN() : std::array<double, M * N>() {}
			MatrixMxN(const std::initializer_list<double>& values)
			{
				if (values.size() != M * N) 
					throw std::invalid_argument("Invalid args number!");
				size_t index = 0;
				for (const auto& v : values)
					this->operator[](index++) = v;
			}
			MatrixMxN(const double(&values)[M][N])
			{
				for (size_t m = 0; m < M; ++m)
					std::memcpy(&this->data()[m * N], values[m], N * sizeof(double));
			}
			MatrixMxN(const MatrixMxN& MatrixMxN) noexcept : std::array<double, M * N>(MatrixMxN) {};
			MatrixMxN(MatrixMxN&& MatrixMxN) noexcept : std::array<double, M * N>(MatrixMxN) {}
			~MatrixMxN() noexcept = default;

			MatrixMxN& operator = (const MatrixMxN& MatrixMxN) noexcept
			{
				std::memcpy(this->data(), MatrixMxN.data(), MatrixMxN.size() * sizeof(double));
				return *this;
			}
			MatrixMxN& operator = (MatrixMxN&& MatrixMxN) noexcept
			{
				this->data() = std::move(MatrixMxN.data());
				return *this;
			}

			size_t rows() const { return M; }
			size_t columns() const { return N; }

			double det() const
			{
				double det{ 1 };
				MatrixMxN<M, N> L, U;
				std::array<size_t, M> P;
				size_t k;
				try {
					PLU(*this, L, U, P, k);
					for (size_t i = 0; i < M; ++i)
						det *= U(i, i);
					if (k % 2 != 0) det = -det;
					return det;
				}
				catch (std::exception ) { return 0.0; }
			}

			const double& operator () (const size_t m, const size_t n) const { return this->operator[](m * N + n); }
			double& operator () (const size_t m, const size_t n) { return this->operator[](m* N + n); }
			
			Vec<N> get_row(size_t index) const
			{
				Vec<N> row;
				for (size_t i = 0; i < N; ++i)
					row[i] = this->operator()(index, i);
				return row;
			}
			Vec<M> get_column(size_t index) const
			{
				Vec<M> column;
				for (size_t i = 0; i < M; ++i)
					column[i] = this->operator()(i, index);
				return column;
			}
			void set_row(size_t index, const Vec<N>& vector)
			{
				for (size_t i = 0; i < N; ++i)
					this->operator()(index, i) = vector[i];
			}
			void set_column(size_t index, const Vec<M>& vector)
			{
				for (size_t i = 0; i < M; ++i)
					this->operator()(i, index) = vector[i];
			}

			MatrixMxN& operator += (const MatrixMxN& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) += matrix[i];
				return *this;
			}
			MatrixMxN& operator -= (const MatrixMxN& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) -= matrix[i];
				return *this;
			}
			MatrixMxN& operator *= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) *= value;
				return *this;
			}
			MatrixMxN& operator /= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) /= value;
				return *this;
			}

			friend MatrixMxN<M, N> operator + (const MatrixMxN<M, N>& first, const MatrixMxN<M, N>& second)
			{
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] + second[i];
				return result;
			}
			friend MatrixMxN<M, N> operator - (const MatrixMxN<M, N>& first, const MatrixMxN<M, N>& second)
			{
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] - second[i];
				return result;
			}
			template<size_t K>
			friend MatrixMxN<M, K> operator * (const MatrixMxN<M, N>& first, const MatrixMxN<N, K>& second)
			{
				MatrixMxN<M, K> result;
				for (size_t m = 0; m < M; ++m) {
					for (size_t n = 0; n < N; ++n) {
						for (size_t k = 0; k < K; ++k)
							result[m * K + k] += first[m * N + n] * second[n * K + k];
					}
				}
				return result;
			}
			friend Vec<M> operator * (const MatrixMxN<M, N>& matrix, const Vec<N>& vector)
			{
				Vec<M> result;
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result[m] += matrix[m * N + n] * vector[n];
				return result;
			}
			friend MatrixMxN<M, N> operator * (const MatrixMxN<M, N>& matrix, const double value)
			{
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend MatrixMxN<M, N> operator * (const double value, const MatrixMxN<M, N>& matrix)
			{
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend MatrixMxN<M, N> operator / (const MatrixMxN<M, N>& matrix, const double value)
			{
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < result.size(); ++i) result[i] = matrix[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const MatrixMxN<M, N>& matrix)
			{
				os << "{ ";
				for (size_t m = 0; m < M; ++m) {
					os << "{ ";
					for (size_t n = 0; n < N; ++n) os << matrix[m * N + n] << "; ";
					os << "} ";
				}
				os << "}";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, MatrixMxN& matrix)
			{
				for (size_t i = 0; i < matrix.size(); ++i) is >> matrix[i];
				return is;
			}

			static MatrixMxN<M, N> identity()
			{
				static_assert(M == N, "MatrixFix is not square!");
				MatrixMxN<M, N> result;
				for (size_t i = 0; i < M; ++i) result[i * N + i] = 1.0;
				return result;
			}
		};

		template<size_t rows, size_t cols>
		MatrixMxN<cols, rows> transpose(const MatrixMxN<rows, cols>& matrix)
		{
			MatrixMxN<cols, rows> result;
			for (size_t m = 0; m < rows; ++m) {
				for (size_t n = 0; n < cols; ++n)
					result(n, m) = matrix(m, n);
			}
			return result;
		}
		template<size_t dim>
		MatrixMxN<dim, dim> inverse(const MatrixMxN<dim, dim>& matrix)
		{
			auto result{ matrix };
			double pivot;// , det{ 1 };
			for (size_t k = 0; k < dim; ++k) {
				pivot = result(k, k);
				//det *= pivot;
				if (std::abs(pivot) < 1e-16)
					throw std::runtime_error("Degenerate matrix!");
				for (size_t m = 0; m < dim; ++m)
					result(m, k) = -result(m, k) / pivot;
				for (size_t m = 0; m < dim; ++m) {
					if (m != k) {
						for (size_t n = 0; n < dim; ++n) {
							if (n != k)
								result(m, n) += result(k, n) * result(m, k);
						}
					}
				}
				for (size_t n = 0; n < dim; ++n) result(k, n) /= pivot;
				result(k, k) = 1 / pivot;
			}
			return result;
		}
		template<size_t dim>
		void LU(const MatrixMxN<dim, dim>& A, MatrixMxN<dim, dim>& L, MatrixMxN<dim, dim>& U)
		{
			double sum;
			for (size_t m = 0; m < dim; ++m) {
				for (size_t n = 0; n < dim; ++n) {
					sum = 0;
					if (m <= n) {
						for (size_t i = 0; i < dim; ++i) sum += L(m, i) * U(i, n);
						U(m, n) = A(m, n) - sum;
					}
					else {
						for (size_t i = 0; i < dim; ++i) sum += L(m, i) * U(i, n);
						L(m, n) = (A(m, n) - sum) / U(n, n);
					}
				}
				L(m, m) = 1.0;
			}
		}
		template<size_t dim>
		void PLU(const MatrixMxN<dim, dim>& A, MatrixMxN<dim, dim>& L, MatrixMxN<dim, dim>& U, std::array<size_t, dim>& P, size_t& k)
		{
			k = 0;
			double buf, absv;
			size_t imax;
			for (size_t i = 0; i < dim; ++i) P[i] = i;
			for (size_t m = 0; m < dim - 1; ++m) {
				buf = std::abs(U(m, m));
				imax = m;
				for (size_t i = m + 1; i < dim; ++i)
					if ((absv = std::abs(U(i, m))) > buf) {
						buf = absv;
						imax = i;
					}
				if (buf < 1e-16) throw std::runtime_error("MatrixFix is degenerate!");
				if (imax != m) {
					for (size_t n = m; n < dim; ++n) std::swap(U(m, n), U(imax, n));
					for (size_t n = 0; n < m; ++n) std::swap(L(m, n), L(imax, n));
					std::swap(P[m], P[imax]);
					k++;
				}
				for (size_t n = m + 1; n < dim; ++n) {
					buf = U(n, m) / U(m, m);
					for (size_t i = m; i < dim; ++i) U(n, i) -= buf * U(m, i);
					L(n, m) = buf;
				}
			}
		}
		template<size_t dim>
		void QR(const MatrixMxN<dim, dim>& A, MatrixMxN<dim, dim>& Q, MatrixMxN<dim, dim>& R)
		{
			std::array<Vec<dim>, dim> vecs;
			for (size_t n = 0; n < dim; ++n) {
				auto vec = vecs[n] = A.get_column(n);
				// Gramme-Shmidte orthogonalization
				for (size_t k = 0; k < n; ++k)
					vecs[n] -= vecs[k] * (vec * vecs[k]);
				Q.set_column(n, vecs[n] = normalize(vecs[n]));
			}
			R = transpose(Q) * A;
		}
		template<size_t dim>
		Vec<dim> solve(const MatrixMxN<dim, dim>& A, const Vec<dim>& B)
		{
			double div{ 1 }, pivot, buf;
			MatrixMxN<dim, dim + 1> D;
			for (size_t m = 0; m < dim; ++m) {
				for (size_t n = 0; n < dim; ++n)
					D(m, n) = A(m, n);
				D(m, dim) = B[m];
			}
			for (size_t k = 0; k < dim; ++k) {
				pivot = D(k, k);
				if (std::abs(pivot) < 1e-16)
					throw std::runtime_error("Degenerate MatrixFix!");
				for (size_t m = 0; m < dim; ++m) {
					if (m != k) {
						for (size_t n = 0; n < dim + 1; ++n) {
							if (n != k)
								D(m, n) = div * (D(m, n) * pivot - D(m, k) * D(k, n));
						}
					}
				}
				div = 1 / pivot;
				for (size_t m = 0; m < dim; ++m) {
					if (m != k) D(m, k) = 0.0;
				}
			}
			return D.get_column(dim) / D[0];
		}
		template<size_t rows, size_t cols>
		MatrixMxN<rows, cols> AxD(const MatrixMxN<rows, cols>& A, const MatrixMxN<cols, cols>& D)
		{
			auto matrix{ A };
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					matrix(m, n) *= D(n, n);
			return matrix;
		}
		template<size_t rows, size_t cols>
		MatrixMxN<rows, cols> DxA(const MatrixMxN<rows, cols>& D, const MatrixMxN<rows, cols>& A)
		{
			auto matrix{ A };
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					matrix(m, n) *= D(m, m);
			return matrix;
		}
		template<size_t rows, size_t cols, size_t count = rows>
		void SVD(const MatrixMxN<rows, cols>& A, MatrixMxN<rows, count>& U, MatrixMxN<count, count>& S, MatrixMxN<count, cols>& V)
		{
			MatrixMxN<cols, cols> M = transpose(A) * A;
			Vec<cols> v;
			MatrixMxN<count, count> Sinv;
			double curr, prev, eps = 1e-16;
			for (size_t i = 0; i < count; ++i) {
				curr = 2.0;
				prev = 1.0;
				for (size_t n = 0; n < cols; ++n) v[n] = 1.0;
				while (std::abs((curr - prev) / curr) > eps) {
					prev = curr;
					v = normalize(M * v);
					curr = v[0];
				}
				curr = v * (M * v);
				// filling the MatrixFix of the eigen vectors
				V.set_row(i, v);
				// filling the diagonal MatrixFix
				S(i, i) = std::sqrt(curr);
				for (size_t m = 0; m < cols; ++m)
					for (size_t n = 0; n < cols; ++n)
						M(m, n) -= curr * v[m] * v[n];
			}
			for (size_t i = 0; i < count; ++i) Sinv(i, i) = 1 / S(i, i);
			U = A * transpose(V) * Sinv;
		}


		class Matrix
		{
		private:
			std::vector<double> _data;
			size_t _rows, _cols;

		public:
			Matrix() noexcept = default;
			Matrix(const size_t rows, const size_t cols) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) } 
			{}
			Matrix(const size_t rows, const size_t cols, const std::initializer_list<double>& list);
			template<size_t rows, size_t cols> Matrix(const double(&arr)[rows][cols]) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) }
			{
				for (size_t m = 0; m < rows; ++m) {
					std::memcpy(_data.data() + m * cols, arr[m], sizeof(double)* cols);
				}
			}
			Matrix(const size_t rows, const size_t cols, const std::vector<double>& vec);
			template<size_t rows, size_t cols> Matrix(const MatrixMxN<rows, cols>& m) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) }
			{
				std::memcpy(_data.data(), m.data(), sizeof(double) * rows * cols);
			}
			Matrix(const Matrix& m) noexcept = default;
			Matrix(Matrix&& m) noexcept;
			~Matrix() = default;

			Matrix& operator = (const Matrix& m) = default;
			Matrix& operator = (Matrix&& m) noexcept;

			size_t rows() const { return _rows; }
			size_t columns() const { return _cols; }

			const double& operator () (const size_t m, const size_t n) const { return _data[m * _cols + n]; }
			double& operator () (const size_t m, const size_t n) { return _data[m * _cols + n]; }

			Vector get_row(const size_t index) const;
			Vector get_column(const size_t index) const;
			void set_row(const size_t index, const Vector& row);
			void set_column(const size_t index, const Vector& column);

			Matrix& operator += (const Matrix& m);
			Matrix& operator -= (const Matrix& m);
			Matrix& operator *= (const double v);
			Matrix& operator /= (const double v);

			friend Matrix operator + (const Matrix& f, const Matrix& s);
			friend Matrix operator - (const Matrix& f, const Matrix& s);
			friend Matrix operator * (const Matrix& f, const Matrix& s);
			friend Vector operator * (const Matrix& m, const Vector& v);
			friend Matrix operator / (const Matrix& m, const double s);
			friend Matrix operator * (const Matrix& m, const double v);
			friend Matrix operator * (const double v, const Matrix& m);

			friend std::ostream& operator << (std::ostream& os, const Matrix& m);
			friend std::istream& operator << (std::istream& is, Matrix& m);

			static Matrix identity(const size_t size);

			friend Matrix transpose(const Matrix& matrix);
			friend Matrix inverse(const Matrix& matrix);
			friend void LU(const Matrix& A, Matrix& L, Matrix& U);
			friend void PLU(const Matrix& A, Matrix& L, Matrix& U, std::vector<size_t>& P, size_t& k);
			friend void QR(const Matrix& A, Matrix& Q, Matrix& R);
			friend Vector solve(const Matrix& A, const Vector& B);
			friend Matrix AxD(const Matrix& A, const Matrix& D);
			friend Matrix DxA(const Matrix& D, const Matrix& A);
			friend void SVD(const Matrix& A, Matrix& U, Matrix& S, Matrix& V);
		};
	}
}