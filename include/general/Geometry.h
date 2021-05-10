#pragma once
#include <ostream>
#include <istream>
#include "Matrix.h"

namespace general
{
	namespace math
	{
		class Vec2 : public Vec<2>
		{
		public:
			Vec2() = default;
			Vec2(const double x, const double y)
			{
				this->operator[](0) = x;
				this->operator[](1) = y;
			}
			Vec2(const Vec<2>& v) : Vec<2>(v) {}
			Vec2(const Vec2& v) : Vec<2>(v) {}
			Vec2(Vec2&& v) noexcept : Vec<2>(v) {}
			~Vec2() = default;

			Vec2& operator = (const Vec<2>& v)
			{
				static_cast<Vec<2>*>(this)->operator=(v);
				return *this;
			}
			Vec2& operator = (const Vec2& v)
			{
				static_cast<Vec<2>*>(this)->operator=(v);
				return *this;
			}
			Vec2& operator = (const Vec2&& v) noexcept
			{
				static_cast<Vec<2>*>(this)->operator=(v);
				return *this;
			}

			const double& x() const noexcept { return this->operator[](0); }
			const double& y() const noexcept { return this->operator[](1); }
			inline double& x() noexcept { return this->operator[](0); }
			inline double& y() noexcept { return this->operator[](1); }
		};

		class Vec3 : public Vec<3>
		{
		public:
			Vec3() = default;
			Vec3(const double x, const double y, const double z)
			{
				this->operator[](0) = x;
				this->operator[](1) = y;
				this->operator[](2) = z;
			}
			Vec3(const Vec<3>& v) : Vec<3>(v) {}
			Vec3(const Vec3& v) : Vec<3>(v) {}
			Vec3(Vec3&& v) noexcept : Vec<3>(v) {}
			~Vec3() = default;

			Vec3& operator = (const Vec<3>& v)
			{
				static_cast<Vec<3>*>(this)->operator=(v);
				return *this;
			}
			Vec3& operator = (const Vec3& v)
			{
				static_cast<Vec<3>*>(this)->operator=(v);
				return *this;
			}
			Vec3& operator = (const Vec3&& v) noexcept
			{
				static_cast<Vec<3>*>(this)->operator=(v);
				return *this;
			}

			constexpr const double& x() const noexcept { return this->operator[](0); }
			constexpr const double& y() const noexcept { return this->operator[](1); }
			constexpr const double& z() const noexcept { return this->operator[](2); }
			inline double& x() noexcept { return this->operator[](0); }
			inline double& y() noexcept { return this->operator[](1); }
			inline double& z() noexcept { return this->operator[](2); }
		};

		Vec3 cross(const Vec3& f, const Vec3& s);

		using Matrix3x3 = MatrixMxN<3, 3>;


		// Struct implements the vector of 6th dimension.
		// Contains 3 position and 3 velocity values.
		struct PV
		{
		private:
			const double& get_by_index(const size_t index) const;
		public:
			Vec3 pos, vel;

			PV() noexcept = default;
			PV(
				const double px,
				const double py,
				const double pz,
				const double vx,
				const double vy,
				const double vz) :
				pos{ Vec3(px, py, pz) }, vel{ Vec3(vx, vy, vz) }
			{}
			PV(
				const Vec3& position,
				const Vec3& velocity) :
				pos{ position }, vel{ velocity }
			{}
			PV(const PV& pv) noexcept = default;
			PV(PV&& pv) noexcept;
			~PV() noexcept = default;

			PV& operator = (const PV& pv) noexcept = default;
			PV& operator = (PV&& pv) noexcept;

			double& operator [] (const size_t index);
			const double& operator [] (const size_t index) const;

			PV& operator += (const PV& pv);
			PV& operator -= (const PV& pv);
			PV& operator *= (const double m);
			PV& operator /= (const double v);

			friend PV operator + (const PV& f, const PV& s);
			friend PV operator - (const PV& f, const PV& s);
			friend PV operator * (const double m, const PV& pv);
			friend PV operator * (const PV& pv, const double m);
			friend PV operator / (const PV& pv, const double v);

			friend std::ostream& operator << (std::ostream& o, const PV& pv);
		};
	}
}