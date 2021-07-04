#include "quaternion.h"
#include <cmath>

namespace general {
	namespace math {
		double quaternion::mod() const
		{
			return std::sqrt(s * s + v * v);
		}
		quaternion& quaternion::operator+=(const quaternion& q)
		{
			s += q.s;
			v += q.v;
			return *this;
		}
		quaternion& quaternion::operator-=(const quaternion& q)
		{
			s -= q.s;
			v -= q.v;
			return *this;
		}
		quaternion& quaternion::operator*=(const quaternion& q)
		{
			const double d = s, x = v[0], y = v[1], z = v[2];
			s = d * q.s - v * q.v;
			v = s * q.v + q.s * v + cross(v, q.v);
			return *this;
		}
		quaternion& quaternion::operator*=(const double n)
		{
			s *= n;
			v *= n;
			return *this;
		}
		quaternion& quaternion::operator/=(const double n)
		{
			s /= n;
			v /= n;
			return *this;
		}
		double angle(const quaternion& f, const quaternion& s)
		{
			return dot(f, s) / f.mod() / s.mod();
		}
		
		std::ostream& operator<<(std::ostream& os, const quaternion& q)
		{
			os << "{ " << q.s << "; " << q.v[0] << "; " << q.v[1] << "; " << q.v[2] << " }";
			return os;
		}
		std::istream& operator>>(std::istream& is, quaternion& q)
		{
			is >> q.s >> q.v[0] >> q.v[1] >> q.v[2];
			return is;
		}
	}
}