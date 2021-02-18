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
            o << "pos: " << pv.Pos << "; vel: " << pv.Vel << ";";
            return o;
        }

        PV::PV(PV&& pv) noexcept
        {
            Pos = std::move(pv.Pos);
            Vel = std::move(pv.Vel);
        }

        PV& PV::operator=(PV&& pv) noexcept
        {
            Pos = std::move(pv.Pos);
            Vel = std::move(pv.Vel);
            return *this;
        }

        PV& PV::operator += (const PV& pv)
        {
            Pos += pv.Pos;
            Vel += pv.Vel;
            return *this;
        }
        PV& PV::operator -= (const PV& pv)
        {
            Pos -= pv.Pos;
            Vel -= pv.Vel;
            return *this;
        }
        PV& PV::operator *= (const double m)
        {
            Pos *= m;
            Vel *= m;
            return *this;
        }
        PV& PV::operator /= (const double v)
        {
            Pos /= v;
            Vel /= v;
            return *this;
        }

        PV operator + (const PV& f, const PV& s)
        {
            return PV(f.Pos + s.Pos, f.Vel + s.Vel);
        }
        PV operator - (const PV& f, const PV& s)
        {
            return PV(f.Pos - s.Pos, f.Vel - s.Vel);
        }
        PV operator * (const PV& pv, const double m)
        {
            return PV(pv.Pos * m, pv.Vel * m);
        }
        PV operator * (const double m, const PV& pv)
        {
            return PV(pv.Pos * m, pv.Vel * m);
        }
        PV operator / (const PV& pv, const double v)
        {
            return PV(pv.Pos / v, pv.Vel / v);
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