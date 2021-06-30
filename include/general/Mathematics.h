#pragma once
#include "GeneralConstants.h"
#include "Quaternion.h"
#include "Matrix.h"

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
        inline constexpr double deg_to_rad(const double degrees) noexcept {
            return degrees * PI / 180.0;
        }
        /// <summary>
        /// radians conversion to degrees
        /// </summary>
        /// <param name="radians">angle in radians</param>
        /// <returns>angle in degrees</returns>
        inline constexpr double rad_to_deg(const double radians) noexcept {
            return radians * 180.0 / PI;
        }
        // Radians to angle seconds conversion
        inline constexpr double rad_to_sec(const double radians) noexcept {
            return rad_to_deg(radians) * SEC_PER_ROUND;
        }
        /// <summary>
        /// angular seconds conversion to radians
        /// </summary>
        /// <param name="seconds">angular seconds</param>
        /// <returns>corresponding angle in radians</returns>
        inline constexpr double sec_to_rad(const double seconds) noexcept {
            return seconds * RAD_PER_SEC;
        }
        /// <summary>
        /// adjust angle in radians to interval [0; 2 * pi]
        /// </summary>
        /// <param name="angle"> - angle in radians</param>
        /// <returns>angle in radians</returns>
        constexpr inline double rad_to_2pi(const double angle) noexcept {
            const auto mult = static_cast<long long>(angle / PI2);
            if (angle > 0.0) return angle - mult * PI2;
            else return angle - mult * PI2 + PI2;
        }
        // A vector transforming by rotation using quaternion.
        // v - the initial vector,
        // q - the quaternion of the rotation.
        constexpr Vec3 rotate_vector(const Vec3& v, const Quaternion& q)
        {
            auto r = q * Quaternion(0, v) * inv(q);
            return r.v;
        }
        // A vector transforming by rotation using angle and the axis of the rotation.
        // v - the initial vector,
        // axis - the direction of the axis of rotation,
        // angle - the angle of rotation in radians.
        Vec3 rotate_vector(const Vec3& v, const Vec3& axis, const double angle);
        // Matrix of the cosines calculation using a quaternion
        inline constexpr Matrix3x3 matrix_from_quaternion(const Quaternion& q) noexcept
        {
            return Matrix3x3{
                { { q.s * q.s + q.v[0] * q.v[0] - q.v[1] * q.v[1] - q.v[2] * q.v[2],
                    2 * (q.v[0] * q.v[1] - q.s * q.v[2]),
                    2 * (q.s * q.v[1] + q.v[0] * q.v[2]) },
                  { 2 * (q.s * q.v[2] + q.v[0] * q.v[1]),
                    q.s * q.s - q.v[0] * q.v[0] + q.v[1] * q.v[1] - q.v[2] * q.v[2],
                    2 * (q.v[1] * q.v[2] - q.s * q.v[0]) },
                  { 2 * (q.v[0] * q.v[2] - q.s * q.v[1]),
                    2 * (q.s * q.v[0] + q.v[1] * q.v[2]),
                    q.s * q.s - q.v[0] * q.v[0] - q.v[1] * q.v[1] + q.v[2] * q.v[2] } }
            };
        }
        // Quaternion calculation using the matrix of the cosines
        Quaternion quaternion_from_matrix(const Matrix3x3& m);
        // Quaternion calculation using Euler angles.
        // angles - (psi, teta, phi) vector.
        Quaternion quaternion_from_eulerangles(const Vec3& angles);
        // Euler angled calculation using quaternion.
        // Returns (psi, teta, phi) vector.
        Vec3 eulerangles_from_quaternion(const Quaternion& q);

        /// <summary>
        /// discrete fourier transform from signal to spectrum
        /// </summary>
        /// <typeparam name="Iter">is data container iterator</typeparam>
        /// <param name="first">is a first element iterator</param>
        /// <param name="last">is a last element iterator</param>
        /// <returns>spectral values</returns>
        template<class Iter> std::vector<Complex> dft(const Iter& first, const Iter& last) {
            size_t num = std::distance(first, last);
            auto x = std::vector<Complex>(num);
            double arg;
            size_t k;
            for (size_t n{}; n < num; ++n) {
                k = 0;
                for (auto iter = first; iter != last; ++iter, ++k) {
                    arg = PI2 * n * k / num;
                    x[n] += (*iter) * Complex { std::cos(arg), -std::sin(arg) };
                }
            }
            return x;
        }
        
        /*void _ifft(
            const std::vector<double>& values, const std::vector<double>& x,
            const size_t lind, const size_t rind) {
            if (rind - lind > 1) {
                const size_t mind{ (lind + rind) / 2 };
                _ifft(values, x, lind, mind);
                _ifft(values, x, mind, rind);
            }
            else {
                for ()
            }
            
        };*/
        
        /// <summary>
        /// inverse discrete fourier transform from spectrum to signal
        /// </summary>
        /// <param name="x"> - spectral values</param>
        /// <returns>signal values</returns>
        inline std::vector<Complex> idft(const std::vector<Complex>& x) {
            auto values = std::vector<Complex>(x.size());
            const double size = static_cast<double>(values.size());
            double arg;
            for (size_t k{}; k < x.size(); ++k) {
                for (size_t n{}; n < x.size(); ++n) {
                    arg = PI2 * n * k / size;
                    values[k] += x[n] * Complex{ std::cos(arg), std::sin(arg) };
                }
                values[k] /= size;
            }
            return values;
        }


        inline constexpr double absolute(const double value) noexcept {
            if (value < 0.0) return -value;
            else return value;
        }
        inline constexpr int sign(const double value) noexcept {
            if (value == 0.0) return 0;
            else if (value < 0.0) return -1;
            else return 1;
        }
        inline constexpr long long roundval(const double value) noexcept {
            const auto integer = static_cast<long long>(value);
            const auto residual = value - integer;
            if (absolute(residual) < 0.5) return integer;
            else return integer + sign(residual);
        }
	}
}