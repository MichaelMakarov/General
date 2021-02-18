#pragma once
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
		struct Vec2
		{
			double X, Y;
			Vec2() : X{ 0 }, Y{ 0 } {}
			Vec2(
				const double x, 
				const double y) : X{ x }, Y{ y }
			{}
			Vec2(const Vec2& v) : X{ v.X }, Y{ v.Y } 
			{}
			Vec2(Vec2&& v) noexcept;
			~Vec2() noexcept = default;

			Vec2& operator = (const Vec2& v) noexcept = default;
			Vec2& operator = (Vec2&& v) noexcept;

			double length() const;

			Vec2& operator += (const Vec2& v);
			Vec2& operator -= (const Vec2& v);
			Vec2& operator /= (const double n);
			Vec2& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& os, const Vec2& v);
			friend std::istream& operator >> (std::istream& is, Vec2& v);

			friend Vec2 operator + (const Vec2& f, const Vec2& s);
			friend Vec2 operator - (const Vec2& f, const Vec2& s);
			friend Vec2 operator * (const Vec2& v, const double n);
			friend Vec2 operator / (const Vec2& v, const double n);
			friend Vec2 operator * (const double n, const Vec2& v);
			friend double operator * (const Vec2& f, const Vec2& s);
		};

		// Struct represents point in the orthogonal coordinate system.
		// Implements adding, substraction, multiplying, can be used for vector implementing.
		struct Vec3
		{
			double X, Y, Z;
			Vec3() noexcept : X(0), Y(0), Z(0) {}
			Vec3(
				const double x,
				const double y,
				const double z) : X(x), Y(y), Z(z)
			{}
			Vec3(const Vec3& xyz) noexcept = default;
			Vec3(Vec3&& xyz) noexcept;
			~Vec3() noexcept = default;

			Vec3& operator = (const Vec3& xyz) noexcept = default;
			Vec3& operator = (Vec3&& xyz) noexcept;

			double length() const;

			Vec3& operator += (const Vec3& v);
			Vec3& operator -= (const Vec3& v);
			Vec3& operator /= (const double n);
			Vec3& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& os, const Vec3& v);
			friend std::istream& operator >> (std::istream& is, Vec3& v);

			friend Vec3 operator + (const Vec3& f, const Vec3& s);
			friend Vec3 operator - (const Vec3& f, const Vec3& s);
			friend Vec3 operator * (const Vec3& v, const double n);
			friend Vec3 operator / (const Vec3& v, const double n);
			friend Vec3 operator * (const double n, const Vec3& v);
			friend double operator * (const Vec3& f, const Vec3& s);

			static Vec3 cross(const Vec3& f, const Vec3& s);
		};

		// Struct implements the vector of 6th dimension.
		// Contains 3 position and 3 velocity values.
		struct PV
		{
			Vec3 Pos, Vel;

			PV() noexcept = default;
			PV(
				const double px,
				const double py,
				const double pz,
				const double vx,
				const double vy,
				const double vz) :
				Pos{ Vec3(px, py, pz) }, Vel{ Vec3(vx, vy, vz) }
			{}
			PV(
				const Vec3& position,
				const Vec3& velocity) :
				Pos{ position }, Vel{ velocity }
			{}
			PV(const PV& pv) noexcept = default;
			PV(PV&& pv) noexcept;
			~PV() noexcept = default;

			PV& operator = (const PV& pv) noexcept = default;
			PV& operator = (PV&& pv) noexcept;

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