#include "Geometry.h"

namespace general
{
    namespace math
    {

        Vec3 cross(const Vec3& f, const Vec3& s)
        {
            return Vec3(
                f[1] * s[2] - f[2] * s[1],
                f[2] * s[0] - f[0] * s[2],
                f[0] * s[1] - f[1] * s[0]
            );
        }


        std::ostream& operator << (std::ostream& o, const PV& pv)
        {
            o << "pos: " << pv.pos << "; vel: " << pv.vel;
            return o;
        }

        PV::PV(PV&& pv) noexcept
        {
            pos = std::move(pv.pos);
            vel = std::move(pv.vel);
        }

        PV& PV::operator=(PV&& pv) noexcept
        {
            pos = std::move(pv.pos);
            vel = std::move(pv.vel);
            return *this;
        }
        const double& PV::get_by_index(const size_t index) const
        {
            switch (index) {
            case 0: return pos.x();
            case 1: return pos.y();
            case 2: return pos.z();
            case 3: return vel.x();
            case 4: return vel.y();
            case 5: return vel.z();
            } throw std::out_of_range("Index is out of range!");
        }
        double& PV::operator [] (const size_t index)
        {
            return const_cast<double&>(get_by_index(index));
        }
        const double& PV::operator [] (const size_t index) const
        {
            return get_by_index(index);
        }

        PV& PV::operator += (const PV& pv)
        {
            pos += pv.pos;
            vel += pv.vel;
            return *this;
        }
        PV& PV::operator -= (const PV& pv)
        {
            pos -= pv.pos;
            vel -= pv.vel;
            return *this;
        }
        PV& PV::operator *= (const double m)
        {
            pos *= m;
            vel *= m;
            return *this;
        }
        PV& PV::operator /= (const double v)
        {
            pos /= v;
            vel /= v;
            return *this;
        }

        PV operator + (const PV& f, const PV& s)
        {
            return PV(f.pos + s.pos, f.vel + s.vel);
        }
        PV operator - (const PV& f, const PV& s)
        {
            return PV(f.pos - s.pos, f.vel - s.vel);
        }
        PV operator * (const PV& pv, const double m)
        {
            return PV(pv.pos * m, pv.vel * m);
        }
        PV operator * (const double m, const PV& pv)
        {
            return PV(pv.pos * m, pv.vel * m);
        }
        PV operator / (const PV& pv, const double v)
        {
            return PV(pv.pos / v, pv.vel / v);
        }

    }
}