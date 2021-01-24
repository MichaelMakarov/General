#include "Times.h"
#include "GeneralConstants.h"
#include "Legendre.h"
#include "Matrix3x3.h"
#include "Vector.h"
#include "Matrix.h"
#include "Geometry.h"
#include <iostream>

using namespace general::geometry;
using namespace general::time;
using namespace general::math;

void test_datetime();
void test_legendre();
void test_matrix3();
void test_vector();
void test_matrix();

int main()
{
	test_datetime();
	test_legendre();
	test_matrix3();
	test_vector();
	test_matrix();
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

	Vector vec1;
	auto vec2{ Vector(10) };
	Vector vec3 = { 2.3, -0.1, 9.8, 12, 4.3, -1.5, 0.9 };
	std::vector<double> array = { 5, 5, 5, 5, 5, 1, 1, 1, 2 };
	Vector vec4{ array };
	auto vec5{ Vector(array.begin(), array.end()) };

	std::cout << "empty vector of size 0: " << vec1 << std::endl;
	std::cout << "empty vector of size 10: " << vec2 << std::endl;
	std::cout << "vector 3: " << vec3 << std::endl;
	std::cout << "vector 4: " << vec4 << std::endl;
	std::cout << "vector 5: " << vec5 << std::endl;

	std::cout << "vector 2 length = " << vec2.length() << std::endl;
	std::cout << "vector 3 length = " << vec3.length() << std::endl;
	std::cout << "vector 3 * vector 4 = " << vec3 * vec4 << std::endl;
}

void test_matrix()
{
	std::cout << "\n...Matrix tests...\n";

	auto list = { 4.0, 5.0, -1.0, -5.6, 10.0, 2.34, -0.31, 3.33, -9.0 };
	auto vector{ std::vector(list) };
	const double array2d[2][5]{ { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 0 } };
	Matrix m1;
	auto m2{ Matrix(5, 4) };
	auto m3{ Matrix(3, 3, vector) };
	auto m4{ Matrix(array2d) };
	auto v{ Vector({2, -1, 0, 0, 3}) };
	const double n{ 0.5 };

	std::cout << "Matrix 1: " << m1 << std::endl;
	std::cout << "Matrix 2: " << m2 << std::endl;
	std::cout << "Matrix 3: " << m3 << std::endl;
	std::cout << "Matrix 4: " << m4 << std::endl;
	std::cout << "Matrix 4 * " << v << " = " << m4 * v << std::endl;
	std::cout << "Matrix 3 * " << n << " : " << m3 * n << std::endl;
	std::cout << "Matrix 3 / " << n << " : " << m3 * n << std::endl;
	std::cout << "Matrix 3 row 1: " << m3.get_row(1) << std::endl;
	std::cout << "Matrix 3 column 2: " << m3.get_column(-1) << std::endl;

}
