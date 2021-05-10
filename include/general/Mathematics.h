#pragma once
#include "GeneralConstants.h"
#include "Quaternion.h"

namespace general
{
	namespace math
	{
        /// <summary>
        /// day part conversion to seconds
        /// </summary>
        /// <param name="t">a part of a day</param>
        /// <returns>seconds</returns>
        inline constexpr double daypart_to_sec(const double t)
        {
            return t * time::SEC_PER_DAY;
        }
        /// <summary>
        /// degrees conversion to radians
        /// </summary>
        /// <param name="degrees">angle in degrees</param>
        /// <returns>angle in radians</returns>
        inline constexpr double deg_to_rad(const double degrees)
        {
            return degrees * PI / 180.0;
        }
        /// <summary>
        /// radians conversion to degrees
        /// </summary>
        /// <param name="radians">angle in radians</param>
        /// <returns>angle in degrees</returns>
        inline constexpr double rad_to_deg(const double radians)
        {
            return radians * 180.0 / PI;
        }
        // Radians to angle seconds conversion
        inline constexpr double rad_to_sec(const double radians)
        {
            return rad_to_deg(radians) * SEC_PER_ROUND;
        }
        /// <summary>
        /// angular seconds conversion to radians
        /// </summary>
        /// <param name="seconds">angular seconds</param>
        /// <returns>corresponding angle in radians</returns>
        inline constexpr double sec_to_rad(const double seconds)
        {
            return seconds * RAD_PER_SEC;
        }
        /// <summary>
        /// adjust angle in radians to interval [0; 2 * pi]
        /// </summary>
        /// <param name="angle"> - angle in radians</param>
        /// <returns>angle in radians</returns>
        inline double rad_to_2pi(const double angle)
        {
            if (angle > 0.0) return angle - std::floor(angle / PI2) * PI2;
            else return angle + std::ceil(-angle / PI2) * PI2;
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