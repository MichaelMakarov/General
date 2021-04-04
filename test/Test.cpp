#include "Times.h"
#include "GeneralConstants.h"
#include "Polynomial.h"
#include "Matrix3x3.h"
#include "Vector.h"
#include "Matrix.h"
#include "Geometry.h"
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

int main()
{
	test_datetime();
	test_polynomials();
	test_matrix3();
	test_vector();
	test_matrix();
	test_quaternion();
	return 0;
}

void test_datetime()
{
	std::cout << "\n...Time tests...\n";

	auto dt0{ DateTime() };
	auto dt1{ DateTime(2000, 1, 3, 23, 23, 23, 423) };
	auto dt2{ DateTime(2000, 1, 3, 0, 0, 1, 0) };
	auto dt3{ DateTime(2013, 12, 23, 1, 56, 34, 100) };
	auto dt4{ DateTime(2013, 12, 23, 11, 5, 7, 0) };
	auto dt5{ DateTime(2013, 12, 23, 20, 23, 19, 0) };

	std::cout << "dt0: " << dt0 << std::endl;
	std::cout << "dt1: " << dt1 << std::endl;
	std::cout << "dt2: " << dt2 << std::endl;
	std::cout << "dt3: " << dt3 << std::endl;
	std::cout << "dt4: " << dt4 << std::endl;
	std::cout << "dt5: " << dt5 << std::endl;

	auto jd0{ JD(JD2000) };
	auto jd1{ JD(dt1) };
	auto jd2{ JD(dt2) };
	auto jd3{ JD(dt3) };
	auto jd4{ JD(dt4) };
	auto jd5{ JD(dt5) };
	auto jd6{ jd5 + 4 / HOURS_PER_DAY };

	std::cout << "JD0 = " << jd0 << std::endl;
	std::cout << "JD1 = " << jd1 << std::endl;
	std::cout << "JD2 = " << jd2 << std::endl;
	std::cout << "JD3 = " << jd3 << std::endl;
	std::cout << "JD4 = " << jd4 << std::endl;
	std::cout << "JD5 = " << jd5 << std::endl;
	std::cout << "JD6 = " << jd6 << std::endl;

	std::cout << "JD0 to DateTime: " << jd0.to_datetime() << std::endl;
	std::cout << "JD1 to DateTime: " << jd1.to_datetime() << std::endl;
	std::cout << "JD2 to DateTime: " << jd2.to_datetime() << std::endl;
	std::cout << "JD3 to DateTime: " << jd3.to_datetime() << std::endl;
	std::cout << "JD4 to DateTime: " << jd4.to_datetime() << std::endl;
	std::cout << "JD5 to DateTime: " << jd5.to_datetime() << std::endl;
	std::cout << "JD6 to DateTime: " << jd6.to_datetime() << std::endl;

	DateTime dt6, dt7;
	if (!try_parse("2010/12/11 23:34:56.231", dt6)) {
		std::cout << "Failed to parse datetime!\n";
	}
	else std::cout << "dt6 parsed: " << dt6 << std::endl;
	if (!try_parse("13/11/2021 13:1:3.111111", dt7, "dd/MM/yyyy hh:m:s.ffffff")) {
		std::cout << "Failed to parse datetime!\n";
	}
	else std::cout << "dt7 parsed: " << dt7 << std::endl;
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

	const auto x{ VectorDyn({ 1.0, 2.0, 3.0, 4.0, 5.0}) };
	const auto y{ VectorDyn({ 0.5, -0.6, 0.52, 3.4, 8.4999 }) };
	auto poly = create_polynom(x, y, 2);
	std::cout << "Differences: ";
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

	const double a[3]{ 0.212340538, 0.590533136, 0.911412040 };
	auto m1{ Matrix3x3() };
	auto m2{ Matrix3x3(1, 2, 3, 1, 2, 4, 3, 2, 1) };
	auto m3{ Matrix3x3::inv(m2) };
	auto m4{ Matrix3x3::eye() };
	auto m5{ Matrix3x3(
		a[0], a[0] * a[0], a[0] * a[0] * a[0],
		a[1], a[1] * a[1], a[1] * a[1] * a[1],
		a[2], a[2] * a[2], a[2] * a[2] * a[2]) };
	auto m6{ Matrix3x3::inv(m5) };
	
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
}

void test_vector()
{
	std::cout << "\n...Vector tests...\n";

	auto vec1{ VectorFix<5>() };
	auto vec2{ VectorFix<10>() };
	VectorFix<7> vec3 = { 2.3, -0.1, 9.8, 12, 4.3, -1.5, 0.9 };
	std::array<double, 9> array = { 5, 5, 5, 5, 5, 1, 1, 1, 2 };
	VectorFix vec4{ array };
	auto vec5{ VectorFix<9>(array.begin(), array.end()) };
	auto vec6{ vec5 };

	std::cout << "empty vector of size 0: " << vec1 << std::endl;
	std::cout << "empty vector of size 10: " << vec2 << std::endl;
	std::cout << "vector 3: " << vec3 << std::endl;
	std::cout << "vector 4: " << vec4 << std::endl;
	std::cout << "vector 5: " << vec5 << std::endl;
	std::cout << "vector 6: " << vec6 << std::endl;

	std::cout << "vector 2 length = " << vec2.length() << std::endl;
	std::cout << "vector 3 length = " << vec3.length() << std::endl;
	std::cout << "vector 3 * vector 4 = " << vec5 * vec4 << std::endl;

	auto vec7{ VectorDyn(4) };
}

void test_matrix()
{
	std::cout << "\n...Matrix tests...\n";

	auto list = { 4.0, 5.0, -1.0, -5.6, 10.0, 2.34, -0.31, 3.33, -9.0 };
	const double array2d[2][5]{ { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 0 } };
	MatrixFix<3, 3> m1;
	auto m2{ MatrixFix<5, 4>() };
	auto m3{ MatrixFix<3, 3>(list) };
	auto m4{ MatrixFix<2, 5>(array2d) };
	auto m5{ m3 };
	auto v{ VectorFix<5>({2, -1, 0, 0, 3}) };
	const double n{ 0.5 };
	auto m6 = MatrixFix<4, 4>::identity();
	auto m7{ MatrixFix<3, 3>({ {2, 1, 3}, {1, 3, -3}, {-2, 4, 4} }) };
	auto m8{ inverse(m7) };

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

	auto q1{ Quaternion(2.0, 3, -1, 0.7) };
	auto q2{ Quaternion(0.9, -3, 0.12, -0.01) };

	std::cout << "q1: " << q1 << std::endl;
	std::cout << "q2: " << q2 << std::endl;
	std::cout << "q1*: " << Quaternion::conj(q1) << std::endl;
	std::cout << "|q1| = " << q1.mod() << std::endl;
	std::cout << "q1 * q2: " << q1 * q2 << std::endl;
	std::cout << "q1^-1: " << Quaternion::inv(q1) << std::endl;
	std::cout << "q2 * 0.5: " << q2 * 0.5 << std::endl;

	auto v{ Vec3(1, 1, 1) };
	auto a{ Vec3(0, 1, 0) };
	auto t{ 45 };
	std::cout << "v: " << v << std::endl;
	std::cout << "a: " << a << std::endl;
	std::cout << v << " rotation around " << a << " by angle " << t << ": " << rotate_vector(v, a, deg_to_rad(t)) << std::endl;

	double sint{ std::sin(deg_to_rad(t)) }, cost{ std::cos(deg_to_rad(t)) };
	auto m{ Matrix3x3(cost, sint, 0, -sint, cost, 0, 0, 0, 1) };
	auto q = quaternion_from_matrix(m);
	std::cout << "initial matrix: " << m << std::endl;
	std::cout << "quaternion from matrix: " << q << std::endl;
	std::cout << "matrix from quaternion: " << matrix_from_quaternion(q) << std::endl;
}