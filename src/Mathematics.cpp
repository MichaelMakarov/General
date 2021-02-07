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
	}
}