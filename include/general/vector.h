#pragma once
#include <array>
#include <vector>
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
		template<size_t N> class vec
		{
			static_assert(N > size_t(0), "Size is equal zero!");
		public:
			double elems[N]{};
		public:
			[[nodiscard]] constexpr double& operator[](const size_t index) {
				return elems[index];
			}
			[[nodiscard]] constexpr const double& operator[](const size_t index) const {
				return elems[index];
			}

			constexpr size_t size() const { return N; }

			double* data() { return elems; }
			const double* data() const { return elems; }

			double length() const noexcept {
				double value{ 0 };
				for (size_t i = 0; i < N; ++i) value += this->operator[](i) * this->operator[](i);
				return std::sqrt(value);
			}

			vec& operator += (const vec& v) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) += v[i];
				return *this;
			}
			vec& operator -= (const vec& v) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) -= v[i];
				return *this;
			}
			vec& operator *= (const double value) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) *= value;
				return *this;
			}
			vec& operator /= (const double value) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) /= value;
				return *this;
			}

			friend constexpr vec operator + (const vec& first, const vec& second) noexcept
			{
				vec result{};
				for (size_t i = 0; i < N; ++i) result[i] = first[i] + second[i];
				return result;
			}
			friend constexpr vec operator - (const vec& first, const vec& second) noexcept
			{
				vec result{};
				for (size_t i = 0; i < N; ++i) result[i] = first[i] - second[i];
				return result;
			}
			friend constexpr double operator * (const vec& first, const vec& second) noexcept
			{
				double result{ 0 };
				for (size_t i = 0; i < N; ++i) result += first[i] * second[i];
				return result;
			}
			friend constexpr vec operator * (const double value, const vec& v) noexcept
			{
				vec result{};
				for (size_t i = 0; i < N; ++i) result[i] = value * v[i];
				return result;
			}
			friend constexpr vec operator * (const vec& v, const double value) noexcept
			{
				vec result{};
				for (size_t i = 0; i < N; ++i) result[i] = value * v[i];
				return result;
			}
			friend constexpr vec operator / (const vec& v, const double value) noexcept
			{
				vec result{};
				for (size_t i = 0; i < N; ++i) result[i] = v[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const vec& v)
			{
				os << "{ ";
				for (size_t i = 0; i < N - 1; ++i) os << v[i] << "; ";
				os << v[N - 1] << " }";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, vec& v)
			{
				for (size_t i = 0; i < N; ++i) is >> v[i];
				return is;
			}
			constexpr static vec<N> ones()
			{
				vec<N> vec{};
				for (size_t i = 0; i < N; ++i) vec[i] = 1.0;
				return vec;
			}

		};

		template<size_t size> vec<size> normalize(const vec<size>& vec) noexcept
		{
			return vec / vec.length();
		}

		template<size_t begin, size_t end, size_t size> constexpr vec<end - begin + 1> slice(const vec<size> v)
		{
			static_assert(end < size && begin < end && (end - begin) < size, "Indices out of bounds!");
			vec<end - begin + 1> vec;
			for (size_t i{ begin }; i <= end; ++i) vec[i - begin] = v[i];
			return vec;
		}

		template<size_t dim1, size_t dim2> vec<dim1 + dim2> constexpr unite(const vec<dim1>& first, const vec<dim2>& second)
		{
			vec<dim1 + dim2> result;
			for (size_t i{}; i < dim1; ++i) result[i] = first[i];
			for (size_t i{}; i < dim2; ++i) result[i + dim1] = second[i];
			return result;
		}

		/// <summary>
		/// 2d vector (x, y)
		/// </summary>
		using vec2 = vec<2>;
		/// <summary>
		/// 3d vector (x, y, z)
		/// </summary>
		using vec3 = vec<3>;
		/// <summary>
		/// vector multiplication
		/// </summary>
		/// <param name="f">is left vector</param>
		/// <param name="s">is right vector</param>
		/// <returns>vector</returns>
		constexpr vec3 cross(const vec3& f, const vec3& s) noexcept
		{
			return vec3{
				f[1] * s[2] - f[2] * s[1],
				f[2] * s[0] - f[0] * s[2],
				f[0] * s[1] - f[1] * s[0]
			};
		}


		class vector : public std::vector<double>
		{
		public:
			vector() : std::vector<double>() {}
			vector(const size_t size) : std::vector<double>(size) {}
			vector(const std::initializer_list<double>& list) : std::vector<double>(list) {}
			vector(const std::vector<double>& vec) : std::vector<double>(vec) {}
			vector(const vector& vec) : std::vector<double>(vec) {}
			vector(vector&& vec) noexcept : std::vector<double>(vec) {}
			~vector() = default;

			vector& operator = (const vector& vec) = default;
			vector& operator = (vector&& vec) noexcept;

			double length() const;

			vector& operator += (const vector& vec);
			vector& operator -= (const vector& vec);
			vector& operator *= (const double value);
			vector& operator /= (const double value);

			friend vector operator + (const vector& first, const vector& second);
			friend vector operator - (const vector& first, const vector& second);
			friend double operator * (const vector& first, const vector& second);
			friend vector operator * (const double value, const vector& vec);
			friend vector operator * (const vector& vec, const double value);
			friend vector operator / (const vector& vec, const double value);
			friend std::ostream& operator <<(std::ostream& os, const vector& vec);
			friend std::istream& operator >>(std::istream& is, vector& vec);

			static vector ones(const size_t size);

		};

		vector normalize(const vector& vec);
	}
}