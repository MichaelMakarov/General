#pragma once
#include "Vector.h"

namespace general
{
	namespace math
	{
		template<size_t M, size_t N>
		class Matrix : public std::array<double, M * N>
		{
		public:
			Matrix() : std::array<double, M * N>() {}
			Matrix(const std::initializer_list<double>& values)
			{
				if (values.size() != M * N) 
					throw std::invalid_argument("Invalid args number!");
				size_t index = 0;
				for (const auto& v : values)
					this->operator[](index++) = v;
			}
			Matrix(const double(&values)[M][N])
			{
				for (size_t m = 0; m < M; ++m)
					std::memcpy(&data()[m * N], values[m], N * sizeof(double));
			}
			Matrix(const Matrix& matrix) noexcept : std::array<double, M * N>(matrix) {};
			Matrix(Matrix&& matrix) noexcept : std::array<double, M * N>(matrix) {}
			~Matrix() noexcept = default;

			Matrix& operator = (const Matrix& matrix) noexcept
			{
				std::memcpy(this->data(), matrix.data(), matrix.size() * sizeof(double));
				return *this;
			}
			Matrix& operator = (Matrix&& matrix) noexcept
			{
				this->data() = std::move(matrix.data());
				return *this;
			}

			size_t rows() const { return M; }
			size_t columns() const { return N; }

			const double& operator () (const size_t m, const size_t n) const { return this->operator[](m * N + n); }
			double& operator () (const size_t m, const size_t n) { return this->operator[](m* N + n); }
			Vector<N> get_row(size_t index) const
			{
				Vector<N> row;
				for (size_t i = 0; i < N; ++i)
					row[i] = this->operator()(index, i);
				return row;
			}
			Vector<M> get_column(size_t index) const
			{
				Vector<M> column;
				for (size_t i = 0; i < M; ++i)
					column[i] = this->operator()(i, index);
				return column;
			}
			void set_row(size_t index, const Vector<N>& vector)
			{
				for (size_t i = 0; i < N; ++i)
					this->operator()(index, i) = vector[i];
			}
			void set_column(size_t index, const Vector<M>& vector)
			{
				for (size_t i = 0; i < M; ++i)
					this->operator()(i, index) = vector[i];
			}

			Matrix& operator += (const Matrix& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) += matrix[i];
				return *this;
			}
			Matrix& operator -= (const Matrix& matrix)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) -= matrix[i];
				return *this;
			}
			Matrix& operator *= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) *= value;
				return *this;
			}
			Matrix& operator /= (const double value)
			{
				for (size_t i = 0; i < this->size(); ++i)
					this->operator[](i) /= value;
				return *this;
			}

			friend Matrix<M, N> operator + (const Matrix<M, N>& first, const Matrix<M, N>& second)
			{
				Matrix<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] + second[i];
				return result;
			}
			friend Matrix<M, N> operator - (const Matrix<M, N>& first, const Matrix<M, N>& second)
			{
				Matrix<M, N> result;
				for (size_t i = 0; i < first.size(); ++i)
					result[i] = first[i] - second[i];
				return result;
			}
			template<size_t K>
			friend Matrix<M, N> operator * (const Matrix<M, K>& first, const Matrix<K, N>& second)
			{
				Matrix<M, N> result;
				for (size_t m = 0; m < M; ++m)
				{
					for (size_t k = 0; k < K; ++k)
					{
						for (size_t n = 0; n < N; ++n)
							result[m * K + n] =	first._pValues[m * K + k] * second._pValues[k * N + n];
					}
				}
				return result;
			}
			friend Vector<M> operator * (const Matrix<M, N>& matrix, const Vector<N>& vector)
			{
				Vector<M> result;
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						result[m] += matrix[m * N + n] * vector[n];
				return result;
			}
			friend Matrix<M, N> operator * (const Matrix<M, N>& matrix, const double value)
			{
				Matrix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend Matrix<M, N> operator * (const double value, const Matrix<M, N>& matrix)
			{
				Matrix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] * value;
				return result;
			}
			friend Matrix<M, N> operator / (const Matrix<M, N>& matrix, const double value)
			{
				Matrix<M, N> result;
				for (size_t i = 0; i < result.size(); ++i)
					result[i] = matrix[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const Matrix<M, N>& matrix)
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
			friend std::istream& operator >>(std::istream& is, Matrix& matrix)
			{
				for (size_t i = 0; i < matrix.size(); ++i)
					is >> matrix[i];
				return is;
			}

			static Matrix<M, N> identity()
			{
				static_assert(M == N, "Matrix is not square!");
				Matrix<M, N> result;
				for (size_t i = 0; i < M; ++i)
					result[i * N + i] = 1.0;
				return result;
			}
		};
	}
}