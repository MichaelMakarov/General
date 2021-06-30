#include "Quaternion.h"
#include <cmath>

namespace general
{
	namespace math
	{
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
			const double d = s, x = v[0], y = v[1], z = v[2];
			s = d * q.s - v * q.v;
			v = s * q.v + q.s * v + cross(v, q.v);
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
		double angle(const Quaternion& f, const Quaternion& s)
		{
			return dot(f, s) / f.mod() / s.mod();
		}
		
		std::ostream& operator<<(std::ostream& os, const Quaternion& q)
		{
			os << "{ " << q.s << "; " << q.v[0] << "; " << q.v[1] << "; " << q.v[2] << " }";
			return os;
		}
		std::istream& operator>>(std::istream& is, Quaternion& q)
		{
			is >> q.s >> q.v[0] >> q.v[1] >> q.v[2];
			return is;
		}
	}
}