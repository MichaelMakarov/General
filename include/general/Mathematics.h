#include "GeneralConstants.h"
#include <istream>
#include "Geometry.h"
#include "Quaternion.h"

namespace general
{
	namespace math
	{
        // Degreese to radians conversion
        inline double deg_to_rad(const double degrees)
        {
            return degrees * PI / 180.0;
        }
        // Radians to degreese conversion
        inline double rad_to_deg(const double radians)
        {
            return radians * 180.0 / PI;
        }
        // Factorial
        long double factorial(const size_t x);
        // A vector transforming by rotation using quaternion.
        // v - the initial vector,
        // q - the quaternion of the rotation.
        geometry::XYZ rotate_vector(const geometry::XYZ& v, const Quaternion& q);
        // A vector transforming by rotation using angle and the axis of the rotation.
        // v - the initial vector,
        // axis - the direction of the axis of rotation,
        // angle - the angle of rotation in radians.
        geometry::XYZ rotate_vector(const geometry::XYZ& v, const geometry::XYZ& axis, const double angle);
	}
}