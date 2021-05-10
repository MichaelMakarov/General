#include "Quaternion.h"
#include <cmath>

namespace general
{
	namespace math
	{
		Quaternion& Quaternion::operator=(Quaternion&& q) noexcept
		{
			s = std::move(q.s);
			v = std::move(q.v);
			return *this;
		}
		double Quaternion::mod() const
		{
			return std::sqrt(s * s + v * v);
		}
		Quaternion& Quaternion::operator+=(const Quaternion& q)
		{
			s += q.s;
			v += q.v;
			return *this;
		}
		Quaternion& Quaternion::operator-=(const Quaternion& q)
		{
			s -= q.s;
			v -= q.v;
			return *this;
		}
		Quaternion& Quaternion::operator*=(const Quaternion& q)
		{
			const double d = s, x = v.x(), y = v.y(), z = v.z();
			s = d * q.s - v * q.v;
			v.x() = d * q.v.x() + q.s * x + y * q.v.z() - z * q.v.y();
			v.y() = d * q.v.y() + q.s * y + z * q.v.x() - x * q.v.z();
			v.z() = d * q.v.z() + q.s * z + x * q.v.y() - y * q.v.x();
			return *this;
		}
		Quaternion& Quaternion::operator*=(const double n)
		{
			s *= n;
			v *= n;
			return *this;
		}
		Quaternion& Quaternion::operator/=(const double n)
		{
			s /= n;
			v /= n;
			return *this;
		}
		Quaternion Quaternion::conj(const Quaternion& q)
		{
			return Quaternion(q.s, (-1) * q.v);
		}
		double Quaternion::dot(const Quaternion& f, const Quaternion& s)
		{
			return f.s * s.s + f.v * s.v;
		}
		double Quaternion::angle(const Quaternion& f, const Quaternion& s)
		{
			return dot(f, s) / f.mod() / s.mod();
		}
		Quaternion Quaternion::inv(const Quaternion& q)
		{
			return Quaternion::conj(q) / (q.s * q.s + q.v * q.v);
		}
		Quaternion operator+(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(f.s + s.s, f.v + s.v);
		}
		Quaternion operator-(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(f.s - s.s, f.v - s.v);
		}
		Quaternion operator*(const Quaternion& f, const Quaternion& s)
		{
			return Quaternion(
				f.s * s.s - f.v * s.v,
				f.s * s.v + s.s * f.v + cross(f.v, s.v)
			);
		}
		Quaternion operator*(const Quaternion& q, const double n)
		{
			return Quaternion(q.s * n, q.v * n);
		}
		Quaternion operator/(const Quaternion& q, const double n)
		{
			return Quaternion(q.s / n, q.v / n);
		}
		std::ostream& operator<<(std::ostream& os, const Quaternion& q)
		{
			os << "{ " << q.s << "; " << q.v.x() << "; " << q.v.y() << "; " << q.v.z() << " }";
			return os;
		}
		std::istream& operator>>(std::istream& is, Quaternion& q)
		{
			is >> q.s >> q.v.x() >> q.v.y() >> q.v.z();
			return is;
		}
	}
}