#include "Mathematics.h"

namespace general
{
	namespace math
	{
        long double factorial(const size_t x)
        {
            static const long double values[41]
            {
                1,
                1,
                2,
                6,
                24,
                120,
                720,
                5040,
                40320,
                362880,
                3628800,
                39916800,
                479001600,
                6227020800,
                87178291200,
                1307674368000,
                20922789888000,
                355687428096000,
                6402373705728000,
                121645100408832000,
                2432902008176640000.0,
                5.109094217170944e+19L,
                1.1240007277776077e+21L,
                2.5852016738884978e+22L,
                6.2044840173323941e+23L,
                1.5511210043330986e+25L,
                4.0329146112660565e+26L,
                1.0888869450418352e+28L,
                3.0488834461171384e+29L,
                8.8417619937397008e+30L,
                2.6525285981219103e+32L,
                8.2228386541779224e+33L,
                2.6313083693369352e+35L,
                8.6833176188118859e+36L,
                2.9523279903960412e+38L,
                1.0628380765425749e+40L,
                3.8262170755532692e+41L,
                1.4157003179547095e+43L,
                5.3796612082278958e+44L,
                2.0980678712088794e+46L,
                8.3922714848355173e+45L
            };
            if (x < 40) return values[x];
            throw std::invalid_argument("x value is unsupported!");
        }
        geometry::XYZ rotate_vector(const geometry::XYZ& v, const Quaternion& q)
        {
            auto r = q * Quaternion(0, v.X, v.Y, v.Z) * Quaternion::inv(q);
            return geometry::XYZ(r.X, r.Y, r.Z);
        }
        geometry::XYZ rotate_vector(const geometry::XYZ& v, const geometry::XYZ& axis, const double angle)
        {
            const double cost{ std::cos(angle * 0.5) };
            const double sint{ std::sin(angle * 0.5) };
            return rotate_vector(v, Quaternion(cost, sint * axis.X, sint * axis.Y, sint * axis.Z));
        }
        Matrix3x3 matrix_from_quaternion(const Quaternion& q)
        {
            return Matrix3x3(
                q.S * q.S + q.X * q.X - q.Y * q.Y - q.Z * q.Z,
                2 * (q.X * q.Y - q.S * q.Z),
                2 * (q.S * q.Y + q.X * q.Z),
                2 * (q.S * q.Z + q.X * q.Y),
                q.S * q.S - q.X * q.X + q.Y * q.Y - q.Z * q.Z,
                2 * (q.Y * q.Z - q.S * q.X),
                2 * (q.X * q.Z - q.S * q.Y),
                2 * (q.S * q.X + q.Y * q.Z),
                q.S * q.S - q.X * q.X - q.Y * q.Y + q.Z * q.Z
            );
        }
        Quaternion quaternion_from_matrix(const Matrix3x3& m)
        {
            auto q{ Quaternion(0.25, 0.25, 0.25, 0.25) };
            if ((q.S = 0.5 * std::sqrt(1 + m(0, 0) + m(1, 1) + m(2, 2))) > 1e-6) {
                q.X *= (m(2, 1) - m(1, 2)) / q.S;
                q.Y *= (m(0, 2) - m(2, 0)) / q.S;
                q.Z *= (m(1, 0) - m(0, 1)) / q.S;
            }
            else if ((q.X = 0.5 * std::sqrt(1 + m(0, 0) - m(1, 1) - m(2, 2))) > 1e-6) {
                q.S *= (m(2, 1) - m(1, 2)) / q.X;
                q.Y *= (m(1, 0) + m(0, 1)) / q.X;
                q.Z *= (m(2, 0) + m(0, 2)) / q.X;
            }
            else if ((q.Y = 0.5 * std::sqrt(1 - m(0, 0) + m(1, 1) - m(2, 2))) > 1e-6) {
                q.S *= (m(0, 2) - m(2, 0)) / q.Y;
                q.X *= (m(0, 1) + m(1, 0)) / q.Y;
                q.Z *= (m(2, 1) + m(1, 2)) / q.Y;
            }
            else {
                q.Z = 0.5 * std::sqrt(1 - m(0, 0) - m(1, 1) + m(2, 2));
                q.S *= (m(1, 0) - m(0, 1)) / q.Z;
                q.X *= (m(2, 0) + m(0, 2)) / q.Z;
                q.Y *= (m(2, 1) + m(1, 2)) / q.Z;
            }
            return q;
        }
	}
}