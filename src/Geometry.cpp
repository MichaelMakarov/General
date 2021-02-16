#include "Geometry.h"
#include <cmath>

namespace general
{
    namespace math
    {
        Vec3::Vec3(Vec3&& xyz) noexcept
        {
            X = xyz.X;
            Y = xyz.Y;
            Z = xyz.Z;
            xyz.X = xyz.Y = xyz.Z = 0;
        }

        Vec3& Vec3::operator=(Vec3&& xyz) noexcept
        {
            X = xyz.X;
            Y = xyz.Y;
            Z = xyz.Z;
            xyz.X = xyz.Y = xyz.Z = 0;
            return *this;

        }
        double Vec3::length() const
        {
            return std::sqrt(X * X + Y * Y + Z * Z);
        }

        Vec3& Vec3::operator += (const Vec3& v)
        {
            X += v.X;
            Y += v.Y;
            Z += v.Z;
            return *this;
        }
        Vec3& Vec3::operator -= (const Vec3& v)
        {
            X -= v.X;
            Y -= v.Y;
            Z -= v.Z;
            return *this;
        }
        Vec3& Vec3::operator *= (const double n)
        {
            X *= n;
            Y *= n;
            Z *= n;
            return *this;
        }
        Vec3 Vec3::cross(const Vec3& f, const Vec3& s)
        {
            return Vec3(
                f.Y * s.Z - f.Z * s.Y,
                f.Z * s.X - f.X * s.Z,
                f.X * s.Y - f.Y * s.X
            );
        }
        Vec3& Vec3::operator /= (const double n)
        {
            X /= n;
            Y /= n;
            Z /= n;
            return *this;
        }

        Vec3 operator + (const Vec3& f, const Vec3& s)
        {
            return Vec3(f.X + s.X, f.Y + s.Y, f.Z + s.Z);
        }
        Vec3 operator - (const Vec3& f, const Vec3& s)
        {
            return Vec3(f.X - s.X, f.Y - s.Y, f.Z - s.Z);
        }
        Vec3 operator * (const Vec3& v, const double n)
        {
            return Vec3(v.X * n, v.Y * n, v.Z * n);
        }
        Vec3 operator / (const Vec3& v, const double n)
        {
            return Vec3(v.X / n, v.Y / n, v.Z / n);
        }
        Vec3 operator * (const double n, const Vec3& v)
        {
            return Vec3(v.X * n, v.Y * n, v.Z * n);
        }

        double operator*(const Vec3& f, const Vec3& s)
        {
            return f.X * s.X + f.Y * s.Y + f.Z * s.Z;
        }

        std::ostream& operator << (std::ostream& o, const Vec3& v)
        {
            o << "{ " << v.X << "; " << v.Y << "; " << v.Z << " }";
            return o;
        }

        std::istream& operator>>(std::istream& is, Vec3& v)
        {
            is >> v.X >> v.Y;
            return is;
        }

        std::ostream& operator << (std::ostream& o, const PV& pv)
        {
            o << pv.P1 << "; " << pv.P2 << "; " << pv.P3 << "; " <<
                pv.V1 << "; " << pv.V2 << "; " << pv.V3;
            return o;
        }

        PV::PV(PV&& pv) noexcept
        {
            P1 = pv.P1;
            P2 = pv.P2;
            P3 = pv.P3;
            V1 = pv.V1;
            V2 = pv.V2;
            V3 = pv.V3;
            pv.P1 = pv.P2 = pv.P3 = pv.V1 = pv.V2 = pv.V3 = 0;
        }

        PV& PV::operator=(PV&& pv) noexcept
        {
            P1 = pv.P1;
            P2 = pv.P2;
            P3 = pv.P3;
            V1 = pv.V1;
            V2 = pv.V2;
            V3 = pv.V3;
            pv.P1 = pv.P2 = pv.P3 = pv.V1 = pv.V2 = pv.V3 = 0;
            return *this;
        }

        PV& PV::operator += (const PV& pv)
        {
            P1 += pv.P1;
            P2 += pv.P2;
            P3 += pv.P3;
            V1 += pv.V1;
            V2 += pv.V2;
            V3 += pv.V3;
            return *this;
        }
        PV& PV::operator -= (const PV& pv)
        {
            P1 -= pv.P1;
            P2 -= pv.P2;
            P3 -= pv.P3;
            V1 -= pv.V1;
            V2 -= pv.V2;
            V3 -= pv.V3;
            return *this;
        }
        PV& PV::operator *= (const double m)
        {
            P1 *= m;
            P2 *= m;
            P3 *= m;
            V1 *= m;
            V2 *= m;
            V3 *= m;
            return *this;
        }

        PV operator + (const PV& f, const PV& s)
        {
            return PV(
                f.P1 + s.P1,
                f.P2 + s.P2,
                f.P3 + s.P3,
                f.V1 + s.V1,
                f.V2 + s.V2,
                f.V3 + s.V3
            );
        }
        PV operator - (const PV& f, const PV& s)
        {
            return PV(
                f.P1 - s.P1,
                f.P2 - s.P2,
                f.P3 - s.P3,
                f.V1 - s.V1,
                f.V2 - s.V2,
                f.V3 - s.V3
            );
        }
        PV operator * (const PV& pv, const double m)
        {
            return PV(
                pv.P1 * m,
                pv.P2 * m,
                pv.P3 * m,
                pv.V1 * m,
                pv.V2 * m,
                pv.V3 * m
            );
        }
        PV operator * (const double m, const PV& pv)
        {
            return PV(
                pv.P1 * m,
                pv.P2 * m,
                pv.P3 * m,
                pv.V1 * m,
                pv.V2 * m,
                pv.V3 * m
            );
        }
        Vec2::Vec2(Vec2&& v) noexcept
        {
            X = v.X;
            Y = v.Y;
            v.X = v.Y = 0;
        }
        Vec2& Vec2::operator=(Vec2&& v) noexcept
        {
            X = v.X;
            Y = v.Y;
            v.X = v.Y = 0;
            return *this;
        }
        double Vec2::length() const
        {
            return std::sqrt(X * X + Y * Y);
        }
        Vec2& Vec2::operator+=(const Vec2& v)
        {
            X += v.X;
            Y += v.Y;
            return *this;
        }
        Vec2& Vec2::operator-=(const Vec2& v)
        {
            X -= v.X;
            Y -= v.Y;
            return *this;
        }
        Vec2& Vec2::operator/=(const double n)
        {
            X /= n;
            Y /= n;
            return *this;
        }
        Vec2& Vec2::operator*=(const double n)
        {
            X *= n;
            Y *= n;
            return *this;
        }
        std::ostream& operator<<(std::ostream& o, const Vec2& v)
        {
            o << "{ " << v.X << "; " << v.Y << "}";
            return o;
        }
        std::istream& operator>>(std::istream& is, Vec2& v)
        {
            is >> v.X >> v.Y;
            return is;
        }
        Vec2 operator+(const Vec2& f, const Vec2& s)
        {
            return Vec2(f.X + s.X, f.Y + s.Y);
        }
        Vec2 operator-(const Vec2& f, const Vec2& s)
        {
            return Vec2(f.X - s.X, f.Y - s.Y);
        }
        Vec2 operator*(const Vec2& v, const double n)
        {
            return Vec2(v.X * n, v.Y * n);
        }
        Vec2 operator/(const Vec2& v, const double n)
        {
            return Vec2(v.X / n, v.Y / n);
        }
        Vec2 operator*(const double n, const Vec2& v)
        {
            return Vec2(v.X * n, v.Y * n);
        }
        double operator*(const Vec2& f, const Vec2& s)
        {
            return f.X * s.X + f.Y * s.Y;
        }
}
}