#include "GeneralConstants.h"
#include "Quaternion.h"
#include "Matrix3x3.h"

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
        // Radians to angle seconds conversion
        inline double rad_to_sec(const double radians)
        {
            return rad_to_deg(radians) * SEC_PER_ROUND;
        }
        // Angle seconds to radians conversion
        inline double sec_to_rad(const double seconds)
        {
            return seconds * RAD_PER_SEC;
        }
        // Factorial
        long double factorial(const size_t x);
        // A vector transforming by rotation using quaternion.
        // v - the initial vector,
        // q - the quaternion of the rotation.
        Vec3 rotate_vector(const Vec3& v, const Quaternion& q);
        // A vector transforming by rotation using angle and the axis of the rotation.
        // v - the initial vector,
        // axis - the direction of the axis of rotation,
        // angle - the angle of rotation in radians.
        Vec3 rotate_vector(const Vec3& v, const Vec3& axis, const double angle);
        // Matrix of the cosines calculation using a quaternion
        Matrix3x3 matrix_from_quaternion(const Quaternion& q);
        // Quaternion calculation using the matrix of the cosines
        Quaternion quaternion_from_matrix(const Matrix3x3& m);
        // Quaternion calculation using Euler angles.
        // angles - (psi, teta, phi) vector.
        Quaternion quaternion_from_eulerangles(const Vec3& angles);
        // Euler angled calculation using quaternion.
        // Returns (psi, teta, phi) vector.
        Vec3 eulerangles_from_quaternion(const Quaternion& q);
	}
}