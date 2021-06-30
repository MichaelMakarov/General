#include "Mathematics.h"

namespace general
{
	namespace math
	{
        Vec3 rotate_vector(const Vec3& v, const Vec3& axis, const double angle)
        {
            return rotate_vector(v, Quaternion(std::cos(angle * 0.5), std::sin(angle * 0.5) * axis));
        }
        Quaternion quaternion_from_matrix(const Matrix3x3& m)
        {
            auto q = Quaternion{ 0.25, 0.25, 0.25, 0.25 };
            if ((q.s = 0.5 * std::sqrt(1 + m(0, 0) + m(1, 1) + m(2, 2))) > 1e-6) {
                q.v[0] *= (m(2, 1) - m(1, 2)) / q.s;
                q.v[1] *= (m(0, 2) - m(2, 0)) / q.s;
                q.v[2] *= (m(1, 0) - m(0, 1)) / q.s;
            }
            else if ((q.v[0] = 0.5 * std::sqrt(1 + m(0, 0) - m(1, 1) - m(2, 2))) > 1e-6) {
                q.s *= (m(2, 1) - m(1, 2)) / q.v[0];
                q.v[1] *= (m(1, 0) + m(0, 1)) / q.v[0];
                q.v[2] *= (m(2, 0) + m(0, 2)) / q.v[0];
            }
            else if ((q.v[1] = 0.5 * std::sqrt(1 - m(0, 0) + m(1, 1) - m(2, 2))) > 1e-6) {
                q.s *= (m(0, 2) - m(2, 0)) / q.v[1];
                q.v[0] *= (m(0, 1) + m(1, 0)) / q.v[1];
                q.v[2] *= (m(2, 1) + m(1, 2)) / q.v[1];
            }
            else {
                q.v[2] = 0.5 * std::sqrt(1 - m(0, 0) - m(1, 1) + m(2, 2));
                q.s *= (m(1, 0) - m(0, 1)) / q.v[2];
                q.v[0] *= (m(2, 0) + m(0, 2)) / q.v[2];
                q.v[1] *= (m(2, 1) + m(1, 2)) / q.v[2];
            }
            return q;
        }
        Quaternion quaternion_from_eulerangles(const Vec3& angles)
        {
            const double cost_2{ std::cos(angles[1] * 0.5) };
            const double sint_2{ std::sin(angles[1] * 0.5) };
            const double ppp_2{ (angles[0] + angles[2]) * 0.5 };
            const double pmp_2{ (angles[0] - angles[2]) * 0.5 };
            return Quaternion{
                cost_2 * std::cos(ppp_2),
                sint_2 * std::cos(pmp_2),
                sint_2 * std::sin(pmp_2),
                cost_2 * std::sin(ppp_2)
            };
        }
        Vec3 eulerangles_from_quaternion(const Quaternion& q)
        {
            const double u{ std::atan2(q.v[2], q.s) }; // (psi + phi) / 2
            const double v{ std::atan2(q.v[1], q.v[0]) }; // (psi - phi) / 2
            return Vec3{
                u + v,
                2 * std::atan(std::sqrt((q.v[0] * q.v[0] + q.v[1] * q.v[1]) / (q.s * q.s + q.v[2] * q.v[2]))),
                u - v
            };
        }
	}
}