#include "Times.h"
#include "GeneralConstants.h"
#include "Polynomial.h"
#include "Matrix.h"
#include "Quaternion.h"
#include "Mathematics.h"
#include <iostream>

using namespace general::time;
using namespace general::math;

void test_datetime();
void test_polynomials();
void test_matrix3();
void test_vector();
void test_matrix();
void test_quaternion();
void test_complex();



int main()
{
	test_datetime();
	test_polynomials();
	test_matrix3();
	test_vector();
	test_matrix();
	test_quaternion();
	test_complex();
	return 0;
}

void test_datetime()
{
	std::cout << "\n...Time tests...\n";

	constexpr DateTime dt0;
	constexpr auto dt1{ DateTime(2000, 1, 3, 23, 23, 23, 423) };
	constexpr auto dt2{ DateTime(2000, 1, 3, 0, 0, 1, 0) };
	constexpr auto dt3{ DateTime(2013, 12, 23, 1, 56, 34, 100) };
	constexpr auto dt4{ DateTime(2013, 12, 23, 11, 5, 7, 0) };
	constexpr auto dt5{ DateTime(2013, 12, 23, 20, 23, 19, 0) };

	std::cout << "dt0: " << dt0 << std::endl;
	std::cout << "dt1: " << dt1 << std::endl;
	std::cout << "dt2: " << dt2 << std::endl;
	std::cout << "dt3: " << dt3 << std::endl;
	std::cout << "dt4: " << dt4 << std::endl;
	std::cout << "dt5: " << dt5 << std::endl;

	constexpr auto jd0{ JD(JD2000) };
	constexpr auto jd1{ JD(dt1) };
	constexpr auto jd2{ JD(dt2) };
	constexpr auto jd3{ JD(dt3) };
	constexpr auto jd4{ JD(dt4) };
	constexpr auto jd5{ JD(dt5) };
	constexpr auto jd6{ jd5 + 4 * SEC_PER_HOUR };
	constexpr auto jd7{ jd5 - 4 * SEC_PER_HOUR };

	std::cout << "JD0 = " << jd0 << std::endl;
	std::cout << "JD1 = " << jd1 << std::endl;
	std::cout << "JD2 = " << jd2 << std::endl;
	std::cout << "JD3 = " << jd3 << std::endl;
	std::cout << "JD4 = " << jd4 << std::endl;
	std::cout << "JD5 = " << jd5 << std::endl;
	std::cout << "JD6 = " << jd6 << std::endl;
	std::cout << "JD7 = " << jd7 << std::endl;

	std::cout << "JD0 to DateTime: " << jd_to_datetime(jd0) << std::endl;
	std::cout << "JD1 to DateTime: " << jd_to_datetime(jd1) << std::endl;
	std::cout << "JD2 to DateTime: " << jd_to_datetime(jd2) << std::endl;
	std::cout << "JD3 to DateTime: " << jd_to_datetime(jd3) << std::endl;
	std::cout << "JD4 to DateTime: " << jd_to_datetime(jd4) << std::endl;
	std::cout << "JD5 to DateTime: " << jd_to_datetime(jd5) << std::endl;
	std::cout << "JD6 to DateTime: " << jd_to_datetime(jd6) << std::endl;
	std::cout << "JD6 to DateTime: " << jd_to_datetime(jd7) << std::endl;

	std::string dt6 = "2010/12/11 23:34:56.231", dt7 = "13/11/2021 13:1:3.111111";
	std::cout << dt6 << " parsed: " << datetime_from_str(dt6) << std::endl;
	std::cout << dt7 << " parsed: " << datetime_from_str(dt7, "d/M/y h:m:s.f") << std::endl;
}

void test_legendre()
{
	std::cout << "\n...Legendre polynoms tests...\n";

	auto p1{ LegendrePolynomial(3) };
	auto p2{ LegendrePolynomial(10) };
	auto p3{ LegendrePolynomial(15) };
	
	double v1 = -0.5, v2 = 0.5;

	std::cout << "P3 degree: " << p1.degree() << "; P10 degree: " << p2.degree() << std::endl;

	std::cout << "P3(" << v1 << ") = " << p1(v1) << std::endl;
	std::cout << "P10(" << v1 << ") = " << p2(v1) << std::endl;
	std::cout << "P15(" << v1 << ") = " << p3(v1) << std::endl;
	std::cout << "P3(" << v2 << ") = " << p1(v2) << std::endl;
	std::cout << "P10(" << v2 << ") = " << p2(v2) << std::endl;
	std::cout << "P15(" << v2 << ") = " << p3(v2) << std::endl;

	std::cout << "...Legendre functions tests...\n";

	auto f1{ LegendreFunction(3, 1) };
	auto f2{ LegendreFunction(10, 5) };
	auto f3{ LegendreFunction(15, 0) };
	auto f4{ LegendreFunction(3, 2) };
	auto f5{ LegendreFunction(4, 4) };

	std::cout << "P3_1(" << v1 << ") = " << f1(v1) << std::endl;
	std::cout << "P10_5(" << v1 << ") = " << f2(v1) << std::endl;
	std::cout << "P15_0(" << v1 << ") = " << f3(v1) << std::endl;
	std::cout << "P3_2(" << v1 << ") = " << f4(v1) << std::endl;
	std::cout << "P3_1(" << v2 << ") = " << f1(v2) << std::endl;
	std::cout << "P10_5(" << v2 << ") = " << f2(v2) << std::endl;
	std::cout << "P15_0(" << v2 << ") = " << f3(v2) << std::endl;
	std::cout << "P3_2(" << v2 << ") = " << f4(v2) << std::endl;
	std::cout << "P4_4(" << v2 << ") = " << f5(v2) << std::endl;
}

void test_power()
{
	std::cout << "\n...Power polynomial test...\n";

	const double x[]{ 1.0, 2.0, 3.0, 4.0, 5.0 };
	const double y[]{ 0.5, -0.6, 0.52, 3.4, 8.4999 };
	auto poly = lstsq<2>(x, y);
	std::cout << "differences: ";
	for (size_t i = 0; i < 5; ++i) {
		std::cout << "y - p(" << x[i] << ") = " << y[i] - poly(x[i]) << "; ";
	}
	std::cout << std::endl;
}

void test_newtonian()
{
	std::cout << "\n...Newtonian polynomial test...\n";

	const double x[]{ 1.0, 2.0, 3.0, 4.0, 5.0 };
	const double y[]{ 0.5, -0.6, 0.52, 3.4, 8.4999 };
	auto poly{ NewtonianPolynomial(x, y) };
	double arg{ x[0] };
	while (arg <= x[4]) {
		std::cout << "P(" << arg << ") = " << poly(arg) << "; ";
		arg += 0.1;
	}
	std::cout << std::endl;
}

void test_polynomials()
{
	test_legendre();
	test_power();
	test_newtonian();
}

void test_matrix3()
{
	std::cout << "\n...Matrix tests...\n";

	constexpr const double a[3]{ 0.212340538, 0.590533136, 0.911412040 };
	constexpr Matrix3x3 m1;
	constexpr Matrix3x3 m2{ { { 1, 2, 3 }, { 1, 2, 4 }, { 3, 2, 1 } } };
	constexpr auto m3{ inverse(m2) };
	constexpr Matrix3x3 m4{ Matrix3x3::identity() };
	constexpr Matrix3x3 m5{ 
		{ { a[0], a[0] * a[0], a[0] * a[0] * a[0] },
		  { a[1], a[1] * a[1], a[1] * a[1] * a[1] },
		  { a[2], a[2] * a[2], a[2] * a[2] * a[2] } }
	};
	constexpr auto m6{ inverse(m5) };
	constexpr Matrix3x3 m7{ { {2, 1, 3 }, { 1, 3, -3 }, { -2, 4, 4 } } };
	constexpr auto m8 = inverse(m7);
	constexpr auto m9 = transpose(m7);
	constexpr auto m10 = m8 * m7;
	constexpr auto m11 = m5 * m6;
	constexpr auto m12 = m3 * m2;
	
	std::cout << "M1: " << m1 << std::endl;
	std::cout << "M2: " << m2 << std::endl;
	std::cout << "M4: " << m4 << std::endl;
	std::cout << "M3 = M2^-1: " << m3 << std::endl;
	std::cout << "M2 * M3: " << m2 * m3 << std::endl;
	std::cout << "M2 + M4: " << m2 + m4 << std::endl;
	std::cout << "M2 - M4: " << m2 - m4 << std::endl;
	std::cout << "M2 * 3.5: " << m2 * 3.5 << std::endl;
	std::cout << "M5: " << m5 << std::endl;
	std::cout << "M6: " << m6 << std::endl;
	std::cout << "M5 * M6: " << m5 * m6 << std::endl;
	std::cout << "M7: " << m7 << std::endl;
	std::cout << "M8: " << m8 << std::endl;
	std::cout << "M9: " << m9 << std::endl;
	std::cout << "M10: " << m10 << std::endl;
}

void test_vector()
{
	std::cout << "\n...Vector tests...\n";

	constexpr Vec<5> vec1;
	constexpr auto vec2{ Vec<10>() };
	constexpr auto vec3 = Vec<7>{ 2.3, -0.1, 9.8, 12.0, 4.3, -1.5, 0.9 };
	constexpr Vec<9> vec4{ 5, 5, 5, 5, 5, 1, 1, 1, 2 };
	constexpr Vec3 vec7{ 3.4, 0.9, 1.2 };
	Vec<4> vec8{ 4, 5, 6, 7 };
	Vec<4> vec9{}, vec10{};
	vec9 = std::move(vec8);
	vec10 = vec9;

	std::cout << "empty vector of size 5: " << vec1 << std::endl;
	std::cout << "empty vector of size 10: " << vec2 << std::endl;
	std::cout << "vector 3: " << vec3 << std::endl;
	std::cout << "vector 4: " << vec4 << std::endl;
	std::cout << "vector 7: " << vec7 << std::endl;

	std::cout << "vector 2 length = " << vec2.length() << std::endl;
	std::cout << "vector 3 length = " << vec3.length() << std::endl;

}

void test_matrix()
{
	std::cout << "\n...Matrix tests...\n";

	constexpr MatrixMxN<3, 3> m1;
	constexpr MatrixMxN<5, 4> m2;
	constexpr MatrixMxN<3, 3> m3{ { { 4.0, 5.0, -1.0 }, { -5.6, 10.0, 2.34 }, { -0.31, 3.33, -9.0 } } };
	constexpr MatrixMxN<2, 5> m4{ { { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 0 } } };
	constexpr auto m5{ m3 };
	constexpr Vec<5> v{ 2, -1, 0, 0, 3 };
	const double n{ 0.5 };
	constexpr auto m6 = MatrixMxN<4, 4>::identity();
	constexpr MatrixMxN<3, 3> m7{ { { 2, 1, 3 }, { 1, 3, -3 }, { -2, 4, 4 } } };
	constexpr auto m8{ inverse(m7) };
	constexpr auto row1 = m8.get_row(1);

	std::cout << "Matrix1: rows = " << m1.rows() << ", columns = " << m1.columns() << std::endl;
	std::cout << "Matrix 1: " << m1 << std::endl;
	std::cout << "Matrix 2: " << m2 << std::endl;
	std::cout << "Matrix 3: " << m3 << std::endl;
	std::cout << "Matrix 4: " << m4 << std::endl;
	std::cout << "Matrix 4 * " << v << " = " << m4 * v << std::endl;
	std::cout << "Matrix 3 * " << n << " : " << m3 * n << std::endl;
	std::cout << "Matrix 3 / " << n << " : " << m3 * n << std::endl;
	std::cout << "Matrix 3 row 1: " << m3.get_row(1) << std::endl;
	std::cout << "Matrix 3 column 2: " << m3.get_column(2) << std::endl;
	std::cout << "Matrix 5: " << m5 << std::endl;
	std::cout << "Matrix 6 (identity): " << m6 << std::endl;
	std::cout << "Matrix 3 inverse: " << m7 << std::endl;
	std::cout << "Inversion checking: " << m8 * m7 << std::endl;
	std::cout << "Matrix 3 determinant = " << m3.det() << std::endl;
}

void test_quaternion()
{
	std::cout << "\n...Quaternion tests...\n";

	constexpr const Quaternion q1{ 2.0, 3, -1, 0 };
	constexpr const Quaternion q2{ 1, -3, 0, 0 };
	constexpr Vec3 v1{ 1,2,3 }, v2{ 4,3,2 };
	constexpr Quaternion q3 = conj(q1);
	constexpr Quaternion q4{ 1, 2, 3, 4 };
	constexpr Quaternion q5{ 0, cross(v1, v2) };
	constexpr Quaternion q6 = q1 * q2;

	std::cout << "q1: " << q1 << std::endl;
	std::cout << "q2: " << q2 << std::endl;
	std::cout << "q1*: " << q3 << std::endl;
	std::cout << "|q1| = " << q1.mod() << std::endl;
	std::cout << "q1 * q2: " << q1 * q2 << std::endl;
	std::cout << "q1^-1: " << inv(q1) << std::endl;
	std::cout << "q2 * 0.5: " << q2 * 0.5 << std::endl;

	constexpr Vec3 v{ 1, 1, 1 };
	constexpr Vec3 a{ 0, 1, 0 };
	constexpr auto t{ 45 };
	std::cout << "v: " << v << std::endl;
	std::cout << "a: " << a << std::endl;
	std::cout << v << " rotation around " << a << " by angle " << t << ": " << rotate_vector(v, a, deg_to_rad(t)) << std::endl;

	double sint{ std::sin(deg_to_rad(t)) }, cost{ std::cos(deg_to_rad(t)) };
	Matrix3x3 m{ { {cost, sint, 0 },  { -sint, cost, 0 }, { 0, 0, 1 } } };
	auto q = quaternion_from_matrix(m);
	std::cout << "initial matrix: " << m << std::endl;
	std::cout << "quaternion from matrix: " << q << std::endl;
	std::cout << "matrix from quaternion: " << matrix_from_quaternion(q) << std::endl;
}

class Signal {
	double _ampl, _phase, _freq;
public:
	Signal(
		const double ampl = 0.0, 
		const double freq = 0.0, 
		const double phase = 0.0) : _ampl{ ampl }, _freq{ freq }, _phase{ phase }{}
	double operator () (const double t) noexcept {
		return _ampl * std::sin(_freq * t + _phase);
	}

};
std::vector<double> calc_signal(const std::vector<Complex>& params, const double period) {
	auto signals = std::vector<Signal>(params.size());
	for (size_t i{}; i < params.size(); ++i) {
		signals[i] = Signal(params[i].mod(), i / period, params[i].arg());
	}
	auto result = std::vector<double>(params.size());
	for (size_t i{}; i < params.size(); ++i) {
		for (size_t k{}; k < result.size(); ++k) {
			result[i] += signals[k](k);
		}
	}
	return result;
}

void test_complex() {
	std::cout << "\n...complex numbers test...\n";

	constexpr Complex c1;
	constexpr Complex c2{ 2, -3 };
	constexpr Complex c3 = c2 * conj(c2);
	constexpr Complex c4{ -1, 2 };
	constexpr Complex c5 = c2 / c4;
	constexpr Complex c6 = c5 * 4.0;
	Complex c7; c7 = 8.8;
	const auto c8 = asin(c2);
	constexpr Complex c9 = Complex::i() * Complex::i();
	constexpr Complex c10 = zhukovskiy(c2);
	constexpr Complex c11 = -c10;
	constexpr Complex c12 = dlo(c11);
	constexpr Complex c13 = -Complex::i();
	const Complex c14 = log(c13);
	constexpr Complex c15 = 15;

	std::cout << c1 << std::endl <<
		c2 << std::endl <<
		c3 << std::endl <<
		c4 << std::endl <<
		c5 << std::endl <<
		c6 << std::endl <<
		"c2 mod = " << c2.mod() << " ; arg = " << c2.arg() << std::endl <<
		c7 << std::endl <<
		c8 << std::endl <<
		c9 << std::endl <<
		c10 << std::endl <<
		c11 << std::endl <<
		c12 << std::endl <<
		c13 << std::endl <<
		c14 << std::endl;

	std::cout << "\n...fourier transform...\n";

	const auto func = [](const double x) {
		return 3.0 + std::sin(x) + std::cos(x) - std::sin(3 * x) + std::cos(5 * x);
	};
	const auto sinus = [](const double x) {
		return std::sin(x);
	};
	const size_t count{ 256 };
	auto x = std::vector<double>(count);
	for (size_t i{}; i < count; ++i) x[i] = func(static_cast<double>(i));
	auto vals = dft(x.cbegin(), x.cend());
	for (const auto v : vals) std::cout << v << std::endl;
	std::cout << std::endl;
	vals = idft(vals);
	for (const auto v : vals) std::cout << v << std::endl;
	std::cout << std::endl;
}