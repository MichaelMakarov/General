#pragma once
#include <istream>
#include <ostream>
#include "Geometry.h"

namespace general
{
	namespace math
	{
		class Quaternion
		{
		public:
			double s;
			Vec3 v;

		public:
			Quaternion() noexcept : s{ 1 }, v{ Vec3(0, 0, 0) } {}
			Quaternion(const double s0, const double x, const double y, const double z) noexcept :
				s{ s0 }, v{ Vec3(x, y, z) }
			{}
			Quaternion(const double s0, const Vec3& vec) noexcept : s{ s0 }, v{ vec } {}
			Quaternion(const Quaternion& q) noexcept : s{ q.s }, v{ q.v } {}
			Quaternion(Quaternion&& q) noexcept : s{ std::move(q.s) }, v{ std::move(q.v) }{}
			~Quaternion() = default;

			Quaternion& operator = (const Quaternion& q) = default;
			Quaternion& operator = (Quaternion&& q) noexcept;

			double mod() const;

			Quaternion& operator += (const Quaternion& q);
			Quaternion& operator -= (const Quaternion& q);
			Quaternion& operator *= (const Quaternion& q);
			Quaternion& operator *= (const double n);
			Quaternion& operator /= (const double n);

			friend Quaternion operator + (const Quaternion& f, const Quaternion& s);
			friend Quaternion operator - (const Quaternion& f, const Quaternion& s);
			friend Quaternion operator * (const Quaternion& f, const Quaternion& s);
			friend Quaternion operator * (const Quaternion& q, const double n);
			friend Quaternion operator / (const Quaternion& q, const double n);

			friend std::ostream& operator << (std::ostream& os, const Quaternion& q);
			friend std::istream& operator >> (std::istream& is, Quaternion& q);

			static Quaternion conj(const Quaternion& q);
			static double dot(const Quaternion& f, const Quaternion& s);
			static double angle(const Quaternion& f, const Quaternion& s);
			static Quaternion inv(const Quaternion& q);
		};
	}
}