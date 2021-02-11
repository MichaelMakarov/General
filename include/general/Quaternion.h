#pragma once
#include <istream>
#include <ostream>

namespace general
{
	namespace math
	{
		class Quaternion
		{
		public:
			double S, X, Y, Z;

		public:
			Quaternion() noexcept : S{ 1 }, X{ 0 }, Y{ 0 }, Z{ 0 } {}
			Quaternion(const double s, const double x, const double y, const double z) noexcept :
				S{ s }, X{ x }, Y{ y }, Z{ z }
			{}
			Quaternion(const Quaternion& q) noexcept : S{ q.S }, X{ q.X }, Y{ q.Y }, Z{ q.Z } {}
			Quaternion(Quaternion&& q) noexcept : S{ q.S }, X{ q.X }, Y{ q.Y }, Z{ q.Z }
			{
				q.S = q.X = q.Y = q.Z = 0;
			}
			~Quaternion() = default;

			Quaternion& operator = (const Quaternion& q) noexcept;
			Quaternion& operator = (Quaternion&& q) noexcept;

			double module() const;

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