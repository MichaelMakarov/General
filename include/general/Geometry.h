#pragma once
#include <ostream>

namespace general
{
	namespace math
	{
		// Struct represents point in the orthogonal coordinate system.
		// Implements adding, substraction, multiplying, can be used for vector implementing.
		struct Vec3
		{
			double X, Y, Z;
			Vec3() : X(0), Y(0), Z(0) {}
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

			friend std::ostream& operator << (std::ostream& o, const Vec3& v);

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
			double P1, P2, P3, V1, V2, V3;

			PV() : P1(0), P2(0), P3(0), V1(0), V2(0), V3(0)
			{}
			PV(
				const double px,
				const double py,
				const double pz,
				const double vx,
				const double vy,
				const double vz) :
				P1{ px }, P2{ py }, P3{ pz },
				V1{ vx }, V2{ vy }, V3{ vz }
			{}
			PV(
				const Vec3& position,
				const Vec3& velocity) :
				P1{ position.X },
				P2{ position.Y },
				P3{ position.Z },
				V1{ velocity.X },
				V2{ velocity.Y },
				V3{ velocity.Z }
			{}
			PV(const PV& pv) noexcept = default;
			PV(PV&& pv) noexcept;
			~PV() = default;

			PV& operator = (const PV& pv) noexcept = default;
			PV& operator = (PV&& pv) noexcept;

			PV& operator += (const PV& pv);
			PV& operator -= (const PV& pv);
			PV& operator *= (const double m);

			friend PV operator + (const PV& f, const PV& s);
			friend PV operator - (const PV& f, const PV& s);
			friend PV operator * (const double m, const PV& pv);
			friend PV operator * (const PV& pv, const double m);

			friend std::ostream& operator << (std::ostream& o, const PV& pv);
		};
	}
}