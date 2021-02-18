#include "Quaternion.h"
#include <cmath>

namespace general
{
	namespace math
	{
		Quaternion& Quaternion::operator=(const Quaternion& q) noexcept
		{
			S = q.S;
			X = q.X;
			Y = q.Y;
			Z = q.Z;
			return *this;
		}
		Quaternion& Quaternion::operator=(Quaternion&& q) noexcept
		{
			S = q.S;
			X = q.X;
			Y = q.Y;
			Z = q.Z;
			q.S = q.X = q.Y = q.Z = 0;
			return *this;
		}
		double Quaternion::mod() const
		{
			return std::sqrt(S * S + X * X + Y * Y + Z * Z);
		}
		Quaternion& Quaternion::operator+=(const Quaternion& q)
		{
			S += q.S;
			X += q.X;
			Y += q.Y;
			Z += q.Z;
			return *this;
		}
		Quaternion& Quaternion::operator-=(const Quaternion& q)
		{
			S -= q.S;
			X -= q.X;
			Y -= q.Y;
			Z -= q.Z;
			return *this;
		}
		Quaternion& Quaternion::operator*=(const Quaternion& q)
		{
			const double s = S, x = X, y = Y, z = Z;
			S = s * q.S - x * q.X - y * q.Y - z * q.Z;
			X = x * q.S + s * q.X + y * q.Z - z * q.Y;
			Y = y * q.S + s * q.Y + z * q.X - x * q.Z;
			Z = z * q.S + s * q.Z + x * q.Y - y * q.X;
			return *this;
		}
		Quaternion& Quaternion::operator*=(const double n)
		{
			S *= n;
			X *= n;
			Y *= n;
			Z *= n;
			return *this;
		}
		Quaternion& Quaternion::operator/=(const double n)
		{
			S /= n;
			X /= n;
			Y /= n;
			Z /= n;
			return *this;
		}
		Quaternion Quaternion::conj(const Quaternion& q)
		{
			return Quaternion(q.S, -q.X, -q.Y, -q.Z);
		}
		double Quaternion::dot(const Quaternion& f, const Quaternion& s)
		{
			return f.S * s.S + f.X * s.X + f.Y * s.Y + f.Z * s.Z;
		}
		double Quaternion::angle(const Quaternion& f, const Quaternion& s)
		{
			return dot(f, s) / f.mod() / s.mod();
		}
		Quaternion Quaternion::inv(const Quaternion& q)
		{
			return Quaternion::conj(q) / (q.S * q.S + q.X * q.X + q.Y * q.Y + q.Z * q.Z);
		}
		Quaternion operator+(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(f.S + s.S, f.X + s.X, f.Y + s.Y, f.Z + s.Z);
		}
		Quaternion operator-(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(f.S - s.S, f.X - s.X, f.Y - s.Y, f.Z - s.Z);
		}
		Quaternion operator*(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(
				f.S * s.S - f.X * s.X - f.Y * s.Y - f.Z * s.Z,
				f.X * s.S + f.S * s.X + f.Y * s.Z - f.Z * s.Y,
				f.Y * s.S + f.S * s.Y + f.Z * s.X - f.X * s.Z,
				f.Z * s.S + f.S * s.Z + f.X * s.Y - f.Y * s.X
			);
		}
		Quaternion operator*(const Quaternion& q, const double n)
		{
			return Quaternion(q.S * n, q.X * n, q.Y * n, q.Z * n);
		}
		Quaternion operator/(const Quaternion& q, const double n)
		{
			return Quaternion(q.S / n, q.X / n, q.Y / n, q.Z / n);
		}
		std::ostream& operator<<(std::ostream& os, const Quaternion& q)
		{
			os << "{ " << q.S << "; " << q.X << "; " << q.Y << "; " << q.Z << " }";
			return os;
		}
		std::istream& operator>>(std::istream& is, Quaternion& q)
		{
			is >> q.S >> q.X >> q.Y >> q.Z;
			return is;
		}
	}
}