#include "Mathematics.h"

namespace general
{
	namespace math
	{
        Vec3 rotate_vector(const Vec3& v, const Quaternion& q)
        {
            auto r = q * Quaternion(0, v) * Quaternion::inv(q);
            return r.v;
        }
        Vec3 rotate_vector(const Vec3& v, const Vec3& axis, const double angle)
        {
            return rotate_vector(v, Quaternion(std::cos(angle * 0.5), std::sin(angle * 0.5) * axis));
        }
        Matrix3x3 matrix_from_quaternion(const Quaternion& q)
        {
            return Matrix3x3({
                q.s * q.s + q.v.x() * q.v.x() - q.v.y() * q.v.y() - q.v.z() * q.v.z(),
                2 * (q.v.x() * q.v.y() - q.s * q.v.z()),
                2 * (q.s * q.v.y() + q.v.x() * q.v.z()),
                2 * (q.s * q.v.z() + q.v.x() * q.v.y()),
                q.s * q.s - q.v.x() * q.v.x() + q.v.y() * q.v.y() - q.v.z() * q.v.z(),
                2 * (q.v.y() * q.v.z() - q.s * q.v.x()),
                2 * (q.v.x() * q.v.z() - q.s * q.v.y()),
                2 * (q.s * q.v.x() + q.v.y() * q.v.z()),
                q.s * q.s - q.v.x() * q.v.x() - q.v.y() * q.v.y() + q.v.z() * q.v.z()
            });
        }
        Quaternion quaternion_from_matrix(const Matrix3x3& m)
        {
            auto q{ Quaternion(0.25, 0.25, 0.25, 0.25) };
            if ((q.s = 0.5 * std::sqrt(1 + m(0, 0) + m(1, 1) + m(2, 2))) > 1e-6) {
                q.v.x() *= (m(2, 1) - m(1, 2)) / q.s;
                q.v.y() *= (m(0, 2) - m(2, 0)) / q.s;
                q.v.z() *= (m(1, 0) - m(0, 1)) / q.s;
            }
            else if ((q.v.x() = 0.5 * std::sqrt(1 + m(0, 0) - m(1, 1) - m(2, 2))) > 1e-6) {
                q.s *= (m(2, 1) - m(1, 2)) / q.v.x();
                q.v.y() *= (m(1, 0) + m(0, 1)) / q.v.x();
                q.v.z() *= (m(2, 0) + m(0, 2)) / q.v.x();
            }
            else if ((q.v.y() = 0.5 * std::sqrt(1 - m(0, 0) + m(1, 1) - m(2, 2))) > 1e-6) {
                q.s *= (m(0, 2) - m(2, 0)) / q.v.y();
                q.v.x() *= (m(0, 1) + m(1, 0)) / q.v.y();
                q.v.z() *= (m(2, 1) + m(1, 2)) / q.v.y();
            }
            else {
                q.v.z() = 0.5 * std::sqrt(1 - m(0, 0) - m(1, 1) + m(2, 2));
                q.s *= (m(1, 0) - m(0, 1)) / q.v.z();
                q.v.x() *= (m(2, 0) + m(0, 2)) / q.v.z();
                q.v.y() *= (m(2, 1) + m(1, 2)) / q.v.z();
            }
            return q;
        }
        Quaternion quaternion_from_eulerangles(const Vec3& angles)
        {
            const double cost_2{ std::cos(angles.y() * 0.5) };
            const double sint_2{ std::sin(angles.y() * 0.5) };
            const double ppp_2{ (angles.x() + angles.z()) * 0.5 };
            const double pmp_2{ (angles.x() - angles.z()) * 0.5 };
            return Quaternion(
                cost_2 * std::cos(ppp_2),
                sint_2 * std::cos(pmp_2),
                sint_2 * std::sin(pmp_2),
                cost_2 * std::sin(ppp_2)
            );
        }
        Vec3 eulerangles_from_quaternion(const Quaternion& q)
        {
            const double u{ std::atan2(q.v.z(), q.s) }; // (psi + phi) / 2
            const double v{ std::atan2(q.v.y(), q.v.x()) }; // (psi - phi) / 2
            return Vec3(
                u + v,
                2 * std::atan(std::sqrt((q.v.x() * q.v.x() + q.v.y() * q.v.y()) / (q.s * q.s + q.v.z() * q.v.z()))),
                u - v
            );
        }
	}
}