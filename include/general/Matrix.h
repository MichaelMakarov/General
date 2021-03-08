#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		template<size_t M, size_t N> requires NotZero<M> && NotZero<N>
		class MatrixFix : public std::array<double, M * N>
		{
			static_assert(M * N > size_t(0), "MatrixFix size is equal zero!");
		public:
			MatrixFix() : std::array<double, M * N>() {}
			MatrixFix(const std::initializer_list<double>& values)
			{
				if (values.size() != M * N) 
					throw std::invalid_argument("Invalid args number!");
				size_t index = 0;
				for (const auto& v : values)
					this->operator[](index++) = v;
			}
			MatrixFix(const double(&values)[M][N])
			{
				for (size_t m = 0; m < M; ++m)
					std::memcpy(&this->data()[m * N], values[m], N * sizeof(double));
			}
			MatrixFix(const MatrixFix& MatrixFix) noexcept : std::array<double, M * N>(MatrixFix) {};
			MatrixFix(MatrixFix&& MatrixFix) noexcept : std::array<double, M * N>(MatrixFix) {}
			~MatrixFix() noexcept = default;

			MatrixFix& operator = (const MatrixFix& MatrixFix) noexcept
			{
				std::memcpy(this->data(), MatrixFix.data(), MatrixFix.size() * sizeof(double));
				return *this;
			}
			MatrixFix& operator = (MatrixFix&& MatrixFix) noexcept
			{
				this->data() = std::move(MatrixFix.data());
				return *this;
			}

			size_t rows() const { return M; }
			size_t columns() const { return N; }

			double det() const
			{
				double det{ 1 };
				MatrixFix<M, N> L, U;
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
			
			VectorFix<N> get_row(size_t index) const
			{
				VectorFix<N> row;
				for (size_t i = 0; i < N; ++i)
					row[i] = this->operator()(index, i);
				return row;
			}
			VectorFix<M> get_column(size_t index) const
			{
				VectorFix<M> column;
				for (size_t i = 0; i < M; ++i)
					column[i] = this->operator()(i, index);
				return column;
			}
			void set_row(size_t index, const VectorFix<N>& vector)
			{
				for (size_t i = 0; i < N; ++i)
					this->operator()(index, i) = vector[i];
			}
			void set_column(size_t index, const VectorFix<M>& vector)
			{
				for (size_t i = 0; i < M; ++i)
					this->operator()(i, index) = vector[i];
			}

			MatrixFix& operator += (const MatrixFix& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) += matrix[i];
				return *this;
			}
			MatrixFix& operator -= (const MatrixFix& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) -= matrix[i];
				return *this;
			}
			MatrixFix& operator *= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) *= value;
				return *this;
			}
			MatrixFix& operator /= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) /= value;
				return *this;
			}

			friend MatrixFix<M, N> operator + (const MatrixFix<M, N>& first, const MatrixFix<M, N>& second)
			{
				MatrixFix<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] + second[i];
				return result;
			}
			friend MatrixFix<M, N> operator - (const MatrixFix<M, N>& first, const MatrixFix<M, N>& second)
			{
				MatrixFix<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] - second[i];
				return result;
			}
			template<size_t K>
			friend MatrixFix<M, K> operator * (const MatrixFix<M, N>& first, const MatrixFix<N, K>& second)
			{
				MatrixFix<M, K> result;
				for (size_t m = 0; m < M; ++m) {
					for (size_t n = 0; n < N; ++n) {
						for (size_t k = 0; k < K; ++k)
							result[m * K + k] += first[m * N + n] * second[n * K + k];
					}
				}
				return result;
			}
			friend VectorFix<M> operator * (const MatrixFix<M, N>& matrix, const VectorFix<N>& vector)
			{
				VectorFix<M> result;
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result[m] += matrix[m * N + n] * vector[n];
				return result;
			}
			friend MatrixFix<M, N> operator * (const MatrixFix<M, N>& matrix, const double value)
			{
				MatrixFix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend MatrixFix<M, N> operator * (const double value, const MatrixFix<M, N>& matrix)
			{
				MatrixFix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend MatrixFix<M, N> operator / (const MatrixFix<M, N>& matrix, const double value)
			{
				MatrixFix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const MatrixFix<M, N>& matrix)
			{
				os << "{ ";
				for (size_t m = 0; m < M; ++m)
				{
					os << "{ ";
					for (size_t n = 0; n < N; ++n)
						os << matrix[m * N + n] << "; ";
					os << "} ";
				}
				os << "}";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, MatrixFix& matrix)
			{
				for (size_t i = 0; i < matrix.size(); ++i)
					is >> matrix[i];
				return is;
			}

			static MatrixFix<M, N> identity()
			{
				static_assert(M == N, "MatrixFix is not square!");
				MatrixFix<M, N> result;
				for (size_t i = 0; i < M; ++i)
					result[i * N + i] = 1.0;
				return result;
			}
		};

		template<size_t rows, size_t cols>
		MatrixFix<cols, rows> transpose(const MatrixFix<rows, cols>& matrix)
		{
			MatrixFix<cols, rows> result;
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					result(n, m) = matrix(m, n);
			return result;
		}
		template<size_t dim>
		MatrixFix<dim, dim> inverse(const MatrixFix<dim, dim>& matrix)
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
				for (size_t n = 0; n < dim; ++n)
					result(k, n) /= pivot;
				result(k, k) = 1 / pivot;
			}
			return result;
		}
		template<size_t dim>
		void LU(const MatrixFix<dim, dim>& A, MatrixFix<dim, dim>& L, MatrixFix<dim, dim>& U)
		{
			double sum;
			for (size_t m = 0; m < dim; ++m) {
				for (size_t n = 0; n < dim; ++n) {
					sum = 0;
					if (m <= n) {
						for (size_t i = 0; i < dim; ++i)
							sum += L(m, i) * U(i, n);
						U(m, n) = A(m, n) - sum;
					}
					else {
						for (size_t i = 0; i < dim; ++i)
							sum += L(m, i) * U(i, n);
						L(m, n) = (A(m, n) - sum) / U(n, n);
					}
				}
				L(m, m) = 1.0;
			}
		}
		template<size_t dim>
		void PLU(const MatrixFix<dim, dim>& A, MatrixFix<dim, dim>& L, MatrixFix<dim, dim>& U, std::array<size_t, dim>& P, size_t& k)
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
		void QR(const MatrixFix<dim, dim>& A, MatrixFix<dim, dim>& Q, MatrixFix<dim, dim>& R)
		{
			std::array<VectorFix<dim>, dim> vecs;
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
		VectorFix<dim> solve(const MatrixFix<dim, dim>& A, const VectorFix<dim>& B)
		{
			double div{ 1 }, pivot, buf;
			MatrixFix<dim, dim + 1> D;
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
		MatrixFix<rows, cols> AxD(const MatrixFix<rows, cols>& A, const MatrixFix<cols, cols>& D)
		{
			auto matrix{ A };
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					matrix(m, n) *= D(n, n);
			return matrix;
		}
		template<size_t rows, size_t cols>
		MatrixFix<rows, cols> DxA(const MatrixFix<rows, cols>& D, const MatrixFix<rows, cols>& A)
		{
			auto matrix{ A };
			for (size_t m = 0; m < rows; ++m)
				for (size_t n = 0; n < cols; ++n)
					matrix(m, n) *= D(m, m);
			return matrix;
		}
		template<size_t rows, size_t cols, size_t count = rows>
		void SVD(const MatrixFix<rows, cols>& A, MatrixFix<rows, count>& U, MatrixFix<count, count>& S, MatrixFix<count, cols>& V)
		{
			MatrixFix<cols, cols> M = transpose(A) * A;
			VectorFix<cols> v;
			MatrixFix<count, count> Sinv;
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


		class MatrixDyn
		{
		private:
			std::vector<double> _data;
			size_t _rows, _cols;

		public:
			MatrixDyn() noexcept = default;
			MatrixDyn(const size_t rows, const size_t cols) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) } 
			{}
			MatrixDyn(const size_t rows, const size_t cols, const std::initializer_list<double>& list);
			template<size_t rows, size_t cols> MatrixDyn(const double(&arr)[rows][cols]) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) }
			{
				for (size_t m = 0; m < rows; ++m) {
					std::memcpy(_data.data() + m * cols, arr[m], sizeof(double)* cols);
				}
			}
			MatrixDyn(const size_t rows, const size_t cols, const std::vector<double>& vec);
			template<size_t rows, size_t cols> MatrixDyn(const MatrixFix<rows, cols>& m) : 
				_rows{ rows }, _cols{ cols }, _data{ std::vector<double>(rows * cols) }
			{
				std::memcpy(_data.data(), m.data(), sizeof(double) * rows * cols);
			}
			MatrixDyn(const MatrixDyn& m) noexcept = default;
			MatrixDyn(MatrixDyn&& m) noexcept;
			~MatrixDyn() = default;

			MatrixDyn& operator = (const MatrixDyn& m) = default;
			MatrixDyn& operator = (MatrixDyn&& m) noexcept;

			size_t rows() const { return _rows; }
			size_t columns() const { return _cols; }

			const double& operator () (const size_t m, const size_t n) const { return _data[m * _cols + n]; }
			double& operator () (const size_t m, const size_t n) { return _data[m * _cols + n]; }

			VectorDyn get_row(const size_t index) const;
			VectorDyn get_column(const size_t index) const;
			void set_row(const size_t index, const VectorDyn& row);
			void set_column(const size_t index, const VectorDyn& column);

			MatrixDyn& operator += (const MatrixDyn& m);
			MatrixDyn& operator -= (const MatrixDyn& m);
			MatrixDyn& operator *= (const double v);
			MatrixDyn& operator /= (const double v);

			friend MatrixDyn operator + (const MatrixDyn& f, const MatrixDyn& s);
			friend MatrixDyn operator - (const MatrixDyn& f, const MatrixDyn& s);
			friend MatrixDyn operator * (const MatrixDyn& f, const MatrixDyn& s);
			friend VectorDyn operator * (const MatrixDyn& m, const VectorDyn& v);
			friend MatrixDyn operator / (const MatrixDyn& m, const double s);
			friend MatrixDyn operator * (const MatrixDyn& m, const double v);
			friend MatrixDyn operator * (const double v, const MatrixDyn& m);

			friend std::ostream& operator << (std::ostream& os, const MatrixDyn& m);
			friend std::istream& operator << (std::istream& is, MatrixDyn& m);

			static MatrixDyn identity(const size_t size);

			friend MatrixDyn transpose(const MatrixDyn& matrix);
			friend MatrixDyn inverse(const MatrixDyn& matrix);
			friend void LU(const MatrixDyn& A, MatrixDyn& L, MatrixDyn& U);
			friend void PLU(const MatrixDyn& A, MatrixDyn& L, MatrixDyn& U, std::vector<size_t>& P, size_t& k);
			friend void QR(const MatrixDyn& A, MatrixDyn& Q, MatrixDyn& R);
			friend VectorDyn solve(const MatrixDyn& A, const VectorDyn& B);
			friend MatrixDyn AxD(const MatrixDyn& A, const MatrixDyn& D);
			friend MatrixDyn DxA(const MatrixDyn& D, const MatrixDyn& A);
			friend void SVD(const MatrixDyn& A, MatrixDyn& U, MatrixDyn& S, MatrixDyn& V);
		};
	}
}