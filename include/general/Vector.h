#pragma once
#include <array>
#include <vector>
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
		template<size_t N>
		class Vec : public std::array<double, N>
		{
			static_assert(N > size_t(0), "Size is equal zero!");
		public:
			Vec() noexcept : std::array<double, N>() {}
			Vec(const double(&arr)[N]) : std::array<double, N>() {
				std::memcpy(this->data(), arr, N * sizeof(double));
			}
			Vec(const Vec& v) noexcept = default;
			Vec(Vec&& v) noexcept : std::array<double, N>(v) {}
			Vec(const std::array<double, N>& arr) noexcept : std::array<double, N>(arr) {}
			~Vec() = default;

			Vec& operator = (const Vec& v) noexcept
			{
				static_cast<std::array<double, N>*>(this)->operator=(v);
				return *this;
			}
			Vec& operator = (Vec&& v) noexcept {
				static_cast<std::array<double, N>*>(this)->operator=(v);
				return *this;
			}

			double length() const noexcept
			{
				double value{ 0 };
				for (size_t i = 0; i < N; ++i) value += this->operator[](i) * this->operator[](i);
				return std::sqrt(value);
			}

			Vec& operator += (const Vec& v) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) += v[i];
				return *this;
			}
			Vec& operator -= (const Vec& v) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) -= v[i];
				return *this;
			}
			Vec& operator *= (const double value) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) *= value;
				return *this;
			}
			Vec& operator /= (const double value) noexcept
			{
				for (size_t i = 0; i < N; ++i) this->operator[](i) /= value;
				return *this;
			}

			friend Vec operator + (const Vec& first, const Vec& second) noexcept
			{
				Vec result;
				for (size_t i = 0; i < N; ++i) result[i] = first[i] + second[i];
				return result;
			}
			friend Vec operator - (const Vec& first, const Vec& second) noexcept
			{
				Vec result;
				for (size_t i = 0; i < N; ++i) result[i] = first[i] - second[i];
				return result;
			}
			friend double operator * (const Vec& first, const Vec& second) noexcept
			{
				double result{ 0 };
				for (size_t i = 0; i < N; ++i) result += first[i] * second[i];
				return result;
			}
			friend Vec operator * (const double value, const Vec& vec) noexcept
			{
				Vec result;
				for (size_t i = 0; i < N; ++i) result[i] = value * vec[i];
				return result;
			}
			friend Vec operator * (const Vec& vec, const double value) noexcept
			{
				Vec result;
				for (size_t i = 0; i < N; ++i) result[i] = value * vec[i];
				return result;
			}
			friend Vec operator / (const Vec& vec, const double value) noexcept
			{
				Vec result;
				for (size_t i = 0; i < N; ++i) result[i] = vec[i] / value;
				return result;
			}
			friend std::ostream& operator <<(std::ostream& os, const Vec& vec)
			{
				os << "{ ";
				for (size_t i = 0; i < N - 1; ++i) os << vec[i] << "; ";
				os << vec[N - 1] << " }";
				return os;
			}
			friend std::istream& operator >>(std::istream& is, Vec& vec)
			{
				for (size_t i = 0; i < N; ++i) is >> vec[i];
				return is;
			}

			static Vec<N> ones()
			{
				Vec<N> vec;
				for (size_t i = 0; i < N; ++i) vec[i] = 1.0;
				return vec;
			}

		};

		template<size_t size> Vec<size> normalize(const Vec<size>& vec) noexcept
		{
			return vec / vec.length();
		}

		template<size_t begin, size_t end, size_t size> Vec<end - begin + 1> slice(const Vec<size> v)
		{
			static_assert(end < size && begin < end && (end - begin) < size, "Indices out of bounds!");
			Vec<end - begin + 1> vec;
			for (size_t i = begin; i <= end; ++i) vec[i] = v[i];
			return vec;
		}

		/// <summary>
		/// 2d vector (x, y)
		/// </summary>
		using Vec2 = Vec<2>;
		/// <summary>
		/// 3d vector (x, y, z)
		/// </summary>
		using Vec3 = Vec<3>;
		/// <summary>
		/// Vector multiplication
		/// </summary>
		/// <param name="f">is left vector</param>
		/// <param name="s">is right vector</param>
		/// <returns>vector</returns>
		Vec3 cross(const Vec3& f, const Vec3& s) noexcept;


		class Vector : public std::vector<double>
		{
		public:
			Vector() : std::vector<double>() {}
			Vector(const size_t size) : std::vector<double>(size) {}
			Vector(const std::initializer_list<double>& list) : std::vector<double>(list) {}
			Vector(const std::vector<double>& vec) : std::vector<double>(vec) {}
			Vector(const Vector& vec) : std::vector<double>(vec) {}
			Vector(Vector&& vec) noexcept : std::vector<double>(vec) {}
			~Vector() = default;

			Vector& operator = (const Vector& vec) = default;
			Vector& operator = (Vector&& vec) noexcept;

			double length() const;

			Vector& operator += (const Vector& vec);
			Vector& operator -= (const Vector& vec);
			Vector& operator *= (const double value);
			Vector& operator /= (const double value);

			friend Vector operator + (const Vector& first, const Vector& second);
			friend Vector operator - (const Vector& first, const Vector& second);
			friend double operator * (const Vector& first, const Vector& second);
			friend Vector operator * (const double value, const Vector& vec);
			friend Vector operator * (const Vector& vec, const double value);
			friend Vector operator / (const Vector& vec, const double value);
			friend std::ostream& operator <<(std::ostream& os, const Vector& vec);
			friend std::istream& operator >>(std::istream& is, Vector& vec);

			static Vector ones(const size_t size);

		};

		Vector normalize(const Vector& vec);
	}
}