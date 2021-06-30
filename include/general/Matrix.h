#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		constexpr double abs(const double val) noexcept {
			if (val < 0.0) return -val;
			else return val;
		}

		template<size_t M, size_t N> 
		class MatrixMxN
		{
			static_assert(M * N > size_t(0), "MatrixFix size is equal zero!");
		public:
			double elems[M][N]{};
		public:
			constexpr size_t rows() const noexcept { return M; }
			constexpr size_t columns() const noexcept { return N; }

			double** data() { return elems; }
			const double** data() const { return elems; }

			[[nodiscard]] constexpr double& operator() (const size_t m, const size_t n) 
			{ 
				return elems[m][n];
			}
			[[nodiscard]] constexpr const double& operator() (const size_t m, const size_t n) const
			{
				return elems[m][n];
			}

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
			
			constexpr Vec<N> get_row(const size_t index) const
			{
				Vec<N> row;
				for (size_t i = 0; i < N; ++i) row[i] = this->operator()(index, i);
				return row;
			}
			constexpr Vec<M> get_column(const size_t index) const
			{
				Vec<M> column;
				for (size_t i = 0; i < M; ++i) column[i] = this->operator()(i, index);
				return column;
			}
			void set_row(const size_t index, const Vec<N>& vector)
			{
				for (size_t i = 0; i < N; ++i) this->operator()(index, i) = vector[i];
			}
			void set_column(const size_t index, const Vec<M>& vector)
			{
				for (size_t i = 0; i < M; ++i) this->operator()(i, index) = vector[i];
			}

			MatrixMxN& operator += (const MatrixMxN& matrix)
			{
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						this->operator()(m, n) += matrix(m, n);
				return *this;
			}
			MatrixMxN& operator -= (const MatrixMxN& matrix)
			{
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						this->operator()(m, n) -= matrix(m, n);
				return *this;
			}
			MatrixMxN& operator *= (const double value)
			{
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						this->operator()(m, n) *= value;
				return *this;
			}
			MatrixMxN& operator /= (const double value)
			{
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						this->operator()(m, n) /= value;
				return *this;
			}

			friend constexpr MatrixMxN<M, N> operator + (const MatrixMxN<M, N>& first, const MatrixMxN<M, N>& second)
			{
				MatrixMxN<M, N> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result(m, n) = first(m, n) + second(m, n);
				return result;
			}
			friend constexpr MatrixMxN<M, N> operator - (const MatrixMxN<M, N>& first, const MatrixMxN<M, N>& second)
			{
				MatrixMxN<M, N> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result(m, n) = first(m, n) - second(m, n);
				return result;
			}
			template<size_t K>
			friend constexpr MatrixMxN<M, K> operator * (const MatrixMxN<M, N>& first, const MatrixMxN<N, K>& second)
			{
				MatrixMxN<M, K> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						for (size_t k = 0; k < K; ++k)
							result(m, k) += first(m, n) * second(n, k);
				return result;
			}
			friend constexpr Vec<M> operator * (const MatrixMxN<M, N>& matrix, const Vec<N>& vector)
			{
				Vec<M> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result[m] += matrix(m, n) * vector[n];
				return result;
			}
			friend constexpr MatrixMxN<M, N> operator * (const MatrixMxN<M, N>& matrix, const double value)
			{
				MatrixMxN<M, N> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result(m, n) = matrix(m, n) * value;
				return result;
			}
			friend constexpr MatrixMxN<M, N> operator * (const double value, const MatrixMxN<M, N>& matrix)
			{
				MatrixMxN<M, N> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result(m, n) = matrix(m, n) * value;
				return result;
			}
			friend constexpr MatrixMxN<M, N> operator / (const MatrixMxN<M, N>& matrix, const double value)
			{
				MatrixMxN<M, N> result{};
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result(m, n) = matrix(m, n) / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const MatrixMxN<M, N>& matrix)
			{
				os << "{ ";
				for (size_t m = 0; m < M; ++m) {
					os << "{ ";
					for (size_t n = 0; n < N; ++n) os << matrix(m, n) << "; ";
					os << "} ";
				}
				os << "}";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, MatrixMxN<M, N>& matrix)
			{
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						is >> matrix(m, n);
				return is;
			}

			constexpr static MatrixMxN<M, N> identity()
			{
				static_assert(M == N, "Matrix is not square!");
				MatrixMxN<M, N> result{};
				for (size_t i = 0; i < M; ++i) result(i, i) = 1.0;
				return result;
			}
		};

		template<size_t M, size_t N>
		constexpr MatrixMxN<N, M> transpose(const MatrixMxN<M, N>& matrix)
		{
			MatrixMxN<N, M> result;
			for (size_t m = 0; m < M; ++m) {
				for (size_t n = 0; n < N; ++n)
					result(n, m) = matrix(m, n);
			}
			return result;
		}
		template<size_t dim>
		constexpr MatrixMxN<dim, dim> inverse(const MatrixMxN<dim, dim>& matrix)
		{
			const double zero{ 1e-20 };
			auto result{ matrix };
			double pivot;// , det{ 1 };
			size_t descent[dim];
			size_t k;
			for (size_t m = 0; m < dim; ++m) descent[m] = m;
			for (size_t i = 0; i < dim; ++i) {
				k = descent[i];
				pivot = result(k, k);
				//det *= pivot;
				if (abs(pivot) < zero) {
					if (i < dim - 1) {
						std::swap(descent[i], descent[i + 1]);
						k = descent[i];
						pivot = result(k, k);
					}
					else throw std::runtime_error("Degenerate matrix!");
				}
				for (size_t m = 0; m < dim; ++m) result(m, k) /= -pivot;
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
				buf = abs(U(m, m));
				imax = m;
				for (size_t i = m + 1; i < dim; ++i)
					if ((absv = std::abs(U(i, m))) > buf) {
						buf = absv;
						imax = i;
					}
				if (buf < 1e-20) throw std::runtime_error("Matrix is degenerate!");
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
			MatrixMxN<dim, dim + 1> D{};
			size_t descent[dim];
			for (size_t m = 0; m < dim; ++m) {
				for (size_t n = 0; n < dim; ++n)
					D(m, n) = A(m, n);
				D(m, dim) = B[m];
				descent[m] = m;
			}
			for (size_t k : descent) {
				pivot = D(k, k);
				if (abs(pivot) < 1e-20) {
					if (k != dim - 1) {
						std::swap(descent[k], descent[k + 1]);
						++k;
						pivot = result(k, k);
					}
					else throw std::runtime_error("Degenerate matrix!");
				}
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
		constexpr MatrixMxN<rows, cols> AxD(const MatrixMxN<rows, cols>& A, const MatrixMxN<cols, cols>& D)
		{
			auto matrix{ A };
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					matrix(m, n) *= D(n, n);
			return matrix;
		}
		template<size_t rows, size_t cols>
		constexpr MatrixMxN<rows, cols> DxA(const MatrixMxN<rows, cols>& D, const MatrixMxN<rows, cols>& A)
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
			Vec<cols> v{};
			MatrixMxN<count, count> Sinv{};
			double curr, prev, eps = 1e-16;
			for (size_t i = 0; i < count; ++i) {
				curr = 2.0;
				prev = 1.0;
				for (size_t n = 0; n < cols; ++n) v[n] = 1.0;
				while (abs((curr - prev) / curr) > eps) {
					prev = curr;
					v = normalize(M * v);
					curr = v[0];
				}
				curr = v * (M * v);
				// filling the matrix of the eigen vectors
				V.set_row(i, v);
				// filling the diagonal matrix
				S(i, i) = std::sqrt(curr);
				for (size_t m = 0; m < cols; ++m)
					for (size_t n = 0; n < cols; ++n)
						M(m, n) -= curr * v[m] * v[n];
			}
			for (size_t i = 0; i < count; ++i) Sinv(i, i) = 1 / S(i, i);
			U = A * transpose(V) * Sinv;
		}

		using Matrix3x3 = MatrixMxN<3, 3>;

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