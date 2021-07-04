#pragma once
#include "general_constants.h"
#include "vector.h"
#include <istream>
#include <ostream>

namespace general {
	namespace math {

		class complex {
		public:
			double x{}, y{};

		public:
			double mod() const { return std::sqrt(x * x + y * y); }
			double arg() const { return std::atan2(y, x); }

			complex() = default;
			complex(const complex& c) = default;
			complex(complex&& c) noexcept = default;
			constexpr complex(const double re, const double im = 0.0) noexcept : x{ re }, y{ im } {}

			complex& operator = (const complex& c) = default;
			complex& operator = (complex&& c) noexcept = default;

			complex& operator = (const double v) noexcept {
				x = v;
				return *this;
			}

			complex& operator += (const complex& c) noexcept {
				x += c.x;
				y += c.y;
				return *this;
			}
			complex& operator -= (const complex& c) noexcept {
				x -= c.x;
				y -= c.y;
				return *this;
			}
			complex& operator *= (const double c) noexcept {
				x *= c;
				y *= c;
				return *this;
			}
			complex& operator /= (const double c) {
				x /= c;
				y /= c;
				return *this;
			}
			complex& operator *= (const complex& c) noexcept {
				const double xv{ x }, yv{ y };
				x = xv * c.x - yv * c.y;
				y = xv * c.y + yv * c.x;
				return *this;
			}
			complex& operator /= (const complex& c) {
				const double xv{ x }, yv{ y }, mod = c.x * c.x + c.y * c.y;
				x = (xv * c.x + yv * c.y) / mod;
				y = (yv * c.y - xv * c.y) / mod;
				return *this;
			}

			friend constexpr complex operator - (const complex& c) noexcept {
				return complex{ -c.x, -c.y };
			}

			friend constexpr complex operator + (const complex& f, const complex& s) noexcept {
				return complex{ f.x + s.x, f.y + s.y };
			}
			friend constexpr complex operator - (const complex& f, const complex& s) noexcept {
				return complex{ f.x - s.x, f.y - s.y };
			}
			friend constexpr complex operator * (const complex& c, const double v) noexcept {
				return complex{ c.x * v, c.y * v };
			}
			friend constexpr complex operator / (const complex& c, const double v) {
				return complex{ c.x / v, c.y / v };
			}
			friend constexpr complex operator * (const double v, const complex& c) noexcept {
				return complex{ c.x * v, c.y * v };
			}
			friend constexpr complex operator * (const complex& f, const complex& s) noexcept {
				return complex{ f.x * s.x - f.y * s.y, f.x * s.y + f.y * s.x };
			}
			friend constexpr complex operator / (const complex& f, const complex& s) {
				const double mod{ s.x * s.x + s.y * s.y };
				return complex{ (f.x * s.x + f.y * s.y) / mod, (f.y * s.x - f.x * s.y) / mod };
			}
			friend constexpr complex operator / (const double v, const complex& c) {
				const double mod{ c.x * c.x + c.y * c.y };
				return complex{ c.x * v / mod, -c.y * v / mod };
			}
			friend complex sqrt(const complex& c) {
				const double mod = std::sqrt(c.mod()), arg = c.arg() * 0.5;
				return complex{ mod * std::cos(arg), mod * std::sin(arg) };
			}
			friend complex pow(const complex& c, const double v) {
				const double mod = std::pow(c.mod(), v), arg = c.arg() * v;
				return complex{ mod * std::cos(arg), mod * std::sin(arg) };
			}
			friend constexpr complex conj(const complex& c) noexcept {
				return complex{ c.x, -c.y };
			}
			friend complex exp(const complex& c) {
				const double amp{ std::exp(c.x) };
				return complex{ amp * std::cos(c.y), amp * std::sin(c.y) };
			}
			friend complex log(const complex& c) {
				return complex{ std::log(c.mod()), c.arg() };
			}
			friend constexpr complex zhukovskiy(const complex& c) {
				return (c + 1.0 / c) * 0.5;
			}
			friend constexpr complex zhukconj(const complex& c) {
				return (c - 1.0 / c) * 0.5;
			}
			friend constexpr complex dlo(const complex& c) {
				return (c - 1.0) / (c + 1.0);
			}
			friend complex cosh(const complex& c) {
				return zhukovskiy(exp(c));
			}
			friend complex sinh(const complex& c) {
				return zhukconj(exp(c));
			}
			friend complex tanh(const complex& c) {
				return dlo(exp(c * 2.0));
			}
			friend complex cos(const complex& c) {
				return zhukovskiy(exp(i() * c));
			}
			friend complex sin(const complex& c) {
				return zhukconj(exp(i() * c));
			}
			friend complex tan(const complex& c) {
				return dlo(exp(2.0 * i() * c)) * (-i());
			}
			friend complex acos(const complex& c) {
				return -i() * log(c + sqrt(c * c - 1));
			}
			friend complex asin(const complex& c) {
				return PI1_2 - acos(c);
			}
			friend complex atan(const complex& c) {
				return -0.5 * i() * log((i() + c) / (i() - c));
			}

			friend std::ostream& operator << (std::ostream& ostr, const complex& c) {
				ostr << "( " << c.x << "; " << c.y << " )";
				return ostr;
			}
			friend std::istream& operator >> (std::istream& istr, complex& c) {
				istr >> c.x >> c.y;
				return istr;
			}

			constexpr static complex i() noexcept {
				return complex{ 0, 1 };
			}
		};

		class quaternion {
		public:
			double s;
			vec3 v;

		public:
			double mod() const;

			quaternion& operator += (const quaternion& q);
			quaternion& operator -= (const quaternion& q);
			quaternion& operator *= (const quaternion& q);
			quaternion& operator *= (const double n);
			quaternion& operator /= (const double n);

			friend constexpr quaternion operator + (const quaternion& f, const quaternion& s) {
				return quaternion{ f.s + s.s, f.v + s.v };
			}
			friend constexpr quaternion operator-(const quaternion& f, const quaternion& s) {
				return quaternion{ f.s - s.s, f.v - s.v };
			}
			friend constexpr quaternion operator*(const quaternion& f, const quaternion& s) {
				return quaternion{
					f.s * s.s - f.v * s.v,
					f.s * s.v + s.s * f.v + cross(f.v, s.v)
				};
			}
			friend constexpr quaternion operator*(const quaternion& q, const double n) {
				return quaternion{ q.s * n, q.v * n };
			}
			friend constexpr quaternion operator/(const quaternion& q, const double n) {
				return quaternion{ q.s / n, q.v / n };
			}

			friend std::ostream& operator << (std::ostream& os, const quaternion& q);
			friend std::istream& operator >> (std::istream& is, quaternion& q);

			friend constexpr quaternion conj(const quaternion& q) {
				return quaternion{ q.s, -q.v[0], -q.v[1], -q.v[2] };
			}
			friend constexpr double dot(const quaternion& f, const quaternion& s) {
				return f.s * s.s + f.v * s.v;
			}
			friend double angle(const quaternion& f, const quaternion& s);
			friend constexpr quaternion inv(const quaternion& q) {
				return conj(q) / (q.s * q.s + q.v * q.v);
			}
		};
	}
}