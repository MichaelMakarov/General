#pragma once
#include <array>
#include <vector>
#include <ostream>
#include <istream>
#include <type_traits>
#include <concepts>

namespace general
{
	namespace math
	{
		template<size_t T> 
		concept NotZero = requires { { T > 0 }; };

		template<size_t N> requires NotZero<N>
		class VectorFix : public std::array<double, N>
		{
			static_assert(N > size_t(0), "VectorFix size is equal zero!");
		public:
			VectorFix() noexcept : std::array<double, N>() {}
			explicit VectorFix(const size_t size) : std::array<double, N>(size) {}
			VectorFix(const std::initializer_list<double>& values) : std::array<double, N>()
			{
				if (values.size() != N) throw std::out_of_range("Invalid input list!");
				size_t i = 0;
				for (const auto& v : values)
					this->operator[](i++) = v;
			}
			VectorFix(const double(&values)[N]) : std::array<double, N>() {
				std::memcpy(this->data(), values, N * sizeof(double));
			}
			VectorFix(const VectorFix& VectorFix) = default;
			VectorFix(VectorFix&& VectorFix) noexcept : std::array<double, N>(VectorFix) {}
			VectorFix(const std::array<double, N>& arr) : std::array<double, N>(arr) {}
			template<class IterType> VectorFix(const IterType& begin, const IterType& end) : std::array<double, N>() 
			{
				size_t n = end - begin;
				if (n != N) 
					throw std::out_of_range("Invalid number of values!");
				n = 0;
				for (auto iter = begin; iter != end; ++iter)
					this->operator[](n++) = *iter;
			}
			~VectorFix() noexcept = default;

			VectorFix& operator = (const VectorFix& VectorFix) noexcept = default;
			VectorFix& operator = (VectorFix&& VectorFix) noexcept {
				this->data() = std::move(VectorFix.data());
				return *this;
			}

			double length() const
			{
				double value{ 0 };
				for (size_t i = 0; i < N; ++i)
					value += this->operator[](i);
				return std::sqrt(value);
			}

			VectorFix& operator += (const VectorFix& VectorFix)
			{
				for (size_t i = 0; i < N; ++i)
					this->operator[](i) += VectorFix[i];
				return *this;
			}
			VectorFix& operator -= (const VectorFix& VectorFix) 
			{
				for (size_t i = 0; i < N; ++i)
					this->operator[](i) -= VectorFix[i];
			}
			VectorFix& operator *= (const double value) 
			{
				for (size_t i = 0; i < N; ++i)
					this->operator[](i) *= value;
				return *this;
			}
			VectorFix& operator /= (const double value) 
			{
				for (size_t i = 0; i < N; ++i)
					this->operator[](i) /= value;
				return *this;
			}

			friend VectorFix operator + (const VectorFix& first, const VectorFix& second)
			{
				VectorFix result;
				for (size_t i = 0; i < N; ++i)
					result[i] = first[i] + second[i];
				return result;
			}
			friend VectorFix operator - (const VectorFix& first, const VectorFix& second)
			{
				VectorFix result;
				for (size_t i = 0; i < N; ++i)
					result[i] = first[i] - second[i];
				return result;
			}
			friend double operator * (const VectorFix& first, const VectorFix& second)
			{
				double result{ 0 };
				for (size_t i = 0; i < N; ++i)
					result += first[i] * second[i];
				return result;
			}
			friend VectorFix operator * (const double value, const VectorFix& vec)
			{
				VectorFix result;
				for (size_t i = 0; i < N; ++i)
					result[i] = value * vec[i];
				return result;
			}
			friend VectorFix operator * (const VectorFix& vec, const double value)
			{
				VectorFix result;
				for (size_t i = 0; i < N; ++i)
					result[i] = value * vec[i];
				return result;
			}
			friend VectorFix operator / (const VectorFix& vec, const double value)
			{
				VectorFix result;
				for (size_t i = 0; i < N; ++i)
					result[i] = vec[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const VectorFix& vec)
			{
				os << "{ ";
				for (size_t i = 0; i < N - 1; ++i)
					os << vec[i] << "; ";
				os << vec[N - 1] << "}";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, VectorFix& vec)
			{
				for (size_t i = 0; i < N; ++i)
					is >> vec[i];
				return is;
			}

			template<size_t size> static VectorFix<size> ones()
			{
				VectorFix<size> vec;
				for (size_t i = 0; i < size; ++i) vec[i] = 1.0;
				return vec;
			}
		};

		template<size_t size>
		VectorFix<size> normalize(const VectorFix<size>& vec)
		{
			return vec / vec.length();
		}


		class VectorDyn : public std::vector<double>
		{
		public:
			VectorDyn() : std::vector<double>() {}
			VectorDyn(const size_t size) : std::vector<double>(size) {}
			//VectorDyn(const std::initializer_list<double>& list) : std::vector<double>(list) {}
			template<size_t size> VectorDyn(const double(&arr)[size]) : std::vector<double>(size) 
			{
				std::memcpy(this->data(), arr, sizeof(double) * size);
			}
			VectorDyn(const std::vector<double>& vec) : std::vector<double>(vec) {}
			template<size_t size> VectorDyn(const VectorFix<size>& vec) : std::vector<double>(size)
			{
				std::memcpy(this->data(), vec.data(), sizeof(double) * size);
			}
			VectorDyn(const VectorDyn& vec) : std::vector<double>(vec) {}
			VectorDyn(VectorDyn&& vec) noexcept : std::vector<double>(vec) {}
			~VectorDyn() = default;

			VectorDyn& operator = (const VectorDyn& vec) = default;
			VectorDyn& operator = (VectorDyn&& vec) noexcept;

			double length() const;

			VectorDyn& operator += (const VectorDyn& vec);
			VectorDyn& operator -= (const VectorDyn& vec);
			VectorDyn& operator *= (const double value);
			VectorDyn& operator /= (const double value);

			friend VectorDyn operator + (const VectorDyn& first, const VectorDyn& second);
			friend VectorDyn operator - (const VectorDyn& first, const VectorDyn& second);
			friend double operator * (const VectorDyn& first, const VectorDyn& second);
			friend VectorDyn operator * (const double value, const VectorDyn& vec);
			friend VectorDyn operator * (const VectorDyn& vec, const double value);
			friend VectorDyn operator / (const VectorDyn& vec, const double value);
			friend std::ostream& operator <<(std::ostream& os, const VectorDyn& vec);
			friend std::istream& operator >>(std::istream& is, VectorDyn& vec);

			static VectorDyn ones(const size_t size);

		};

		VectorDyn normalize(const VectorDyn& vec);
	}
}