#pragma once
#include <istream>
#include <ostream>
#include "Vector.h"
#include "GeneralConstants.h"

namespace general
{
	namespace math
	{
		class Complex
		{
		public:
			double x{}, y{};

		public:
			double mod() const { return std::sqrt(x * x + y * y); }
			double arg() const { return std::atan2(y, x); }

			Complex() = default;
			Complex(const Complex& c) = default;
			Complex(Complex&& c) noexcept = default;
			constexpr Complex(const double re, const double im = 0.0) noexcept : x{ re }, y{ im } {}

			Complex& operator = (const Complex& c) = default;
			Complex& operator = (Complex&& c) noexcept = default;

			Complex& operator = (const double v) noexcept {
				x = v;
				return *this;
			}

			Complex& operator += (const Complex& c) noexcept {
				x += c.x;
				y += c.y;
				return *this;
			}
			Complex& operator -= (const Complex& c) noexcept {
				x -= c.x;
				y -= c.y;
				return *this;
			}
			Complex& operator *= (const double c) noexcept {
				x *= c;
				y *= c;
				return *this;
			}
			Complex& operator /= (const double c) {
				x /= c;
				y /= c;
				return *this;
			}
			Complex& operator *= (const Complex& c) noexcept {
				const double xv{ x }, yv{ y };
				x = xv * c.x - yv * c.y;
				y = xv * c.y + yv * c.x;
				return *this;
			}
			Complex& operator /= (const Complex& c) {
				const double xv{ x }, yv{ y }, mod = c.x * c.x + c.y * c.y;
				x = (xv * c.x + yv * c.y) / mod;
				y = (yv * c.y - xv * c.y) / mod;
				return *this;
			}

			friend constexpr Complex operator - (const Complex& c) noexcept {
				return Complex{ -c.x, -c.y };
			}

			friend constexpr Complex operator + (const Complex& f, const Complex& s) noexcept {
				return Complex{ f.x + s.x, f.y + s.y };
			}
			friend constexpr Complex operator - (const Complex& f, const Complex& s) noexcept {
				return Complex{ f.x - s.x, f.y - s.y };
			}
			friend constexpr Complex operator * (const Complex& c, const double v) noexcept {
				return Complex{ c.x * v, c.y * v };
			}
			friend constexpr Complex operator / (const Complex& c, const double v) {
				return Complex{ c.x / v, c.y / v };
			}
			friend constexpr Complex operator * (const double v, const Complex& c) noexcept {
				return Complex{ c.x * v, c.y * v };
			}
			friend constexpr Complex operator * (const Complex& f, const Complex& s) noexcept {
				return Complex{ f.x * s.x - f.y * s.y, f.x * s.y + f.y * s.x };
			}
			friend constexpr Complex operator / (const Complex& f, const Complex& s) {
				const double mod{ s.x * s.x + s.y * s.y };
				return Complex{ (f.x * s.x + f.y * s.y) / mod, (f.y * s.x - f.x * s.y) / mod };
			}
			friend constexpr Complex operator / (const double v, const Complex& c) {
				const double mod{ c.x * c.x + c.y * c.y };
				return Complex{ c.x * v / mod, -c.y * v / mod };
			}
			friend Complex sqrt(const Complex& c) {
				const double mod = std::sqrt(c.mod()), arg = c.arg() * 0.5;
				return Complex{ mod * std::cos(arg), mod * std::sin(arg) };
			}
			friend Complex pow(const Complex& c, const double v) {
				const double mod = std::pow(c.mod(), v), arg = c.arg() * v;
				return Complex{ mod * std::cos(arg), mod * std::sin(arg) };
			}
			friend constexpr Complex conj(const Complex& c) noexcept {
				return Complex{ c.x, -c.y };
			}
			friend Complex exp(const Complex& c) {
				const double amp{ std::exp(c.x) };
				return Complex{ amp * std::cos(c.y), amp * std::sin(c.y) };
			}
			friend Complex log(const Complex& c) {
				return Complex{ std::log(c.mod()), c.arg() };
			}
			friend constexpr Complex zhukovskiy(const Complex& c) {
				return (c + 1.0 / c) * 0.5;
			}
			friend constexpr Complex zhukconj(const Complex& c) {
				return (c - 1.0 / c) * 0.5;
			}
			friend constexpr Complex dlo(const Complex& c) {
				return (c - 1.0) / (c + 1.0);
			}
			friend Complex cosh(const Complex& c) {
				return zhukovskiy(exp(c));
			}
			friend Complex sinh(const Complex& c) {
				return zhukconj(exp(c));
			}
			friend Complex tanh(const Complex& c) {
				return dlo(exp(c * 2.0));
			}
			friend Complex cos(const Complex& c) {
				return zhukovskiy(exp(i() * c));
			}
			friend Complex sin(const Complex& c) {
				return zhukconj(exp(i() * c));
			}
			friend Complex tan(const Complex& c) {
				return dlo(exp(2.0 * i() * c)) * (-i());
			}
			friend Complex acos(const Complex& c) {
				return -i() * log(c + sqrt(c * c - 1));
			}
			friend Complex asin(const Complex& c) {
				return PI1_2 - acos(c);
			}
			friend Complex atan(const Complex& c) {
				return -0.5 * i() * log((i() + c) / (i() - c));
			}

			friend std::ostream& operator << (std::ostream& ostr, const Complex& c) {
				ostr << "( " << c.x << "; " << c.y << " )";
				return ostr;
			}
			friend std::istream& operator >> (std::istream& istr, Complex& c) {
				istr >> c.x >> c.y;
				return istr;
			}

			constexpr static Complex i() noexcept {
				return Complex{ 0, 1 };
			}
		};

		class Quaternion
		{
		public:
			double s;
			Vec3 v;

		public:
			double mod() const;

			Quaternion& operator += (const Quaternion& q);
			Quaternion& operator -= (const Quaternion& q);
			Quaternion& operator *= (const Quaternion& q);
			Quaternion& operator *= (const double n);
			Quaternion& operator /= (const double n);

			friend constexpr Quaternion operator + (const Quaternion& f, const Quaternion& s)
			{
				return Quaternion{ f.s + s.s, f.v + s.v };
			}
			friend constexpr Quaternion operator-(const Quaternion& f, const Quaternion& s)
			{
				return Quaternion{ f.s - s.s, f.v - s.v };
			}
			friend constexpr Quaternion operator*(const Quaternion& f, const Quaternion& s)
			{
				return Quaternion{
					f.s * s.s - f.v * s.v,
					f.s * s.v + s.s * f.v + cross(f.v, s.v)
				};
			}
			friend constexpr Quaternion operator*(const Quaternion& q, const double n)
			{
				return Quaternion{ q.s * n, q.v * n };
			}
			friend constexpr Quaternion operator/(const Quaternion& q, const double n)
			{
				return Quaternion{ q.s / n, q.v / n };
			}

			friend std::ostream& operator << (std::ostream& os, const Quaternion& q);
			friend std::istream& operator >> (std::istream& is, Quaternion& q);

			friend constexpr Quaternion conj(const Quaternion& q)
			{
				return Quaternion{ q.s, -q.v[0], -q.v[1], -q.v[2] };
			}
			friend constexpr double dot(const Quaternion& f, const Quaternion& s)
			{
				return f.s * s.s + f.v * s.v;
			}
			friend double angle(const Quaternion& f, const Quaternion& s);
			friend constexpr Quaternion inv(const Quaternion& q)
			{
				return conj(q) / (q.s * q.s + q.v * q.v);
			}
		};
	}
}