#pragma once
#include "mathfuncs.h"
#include <ostream>
#include <istream>
#include <string>
#include <chrono>

namespace general
{
	namespace time
	{
		using llong = long long;
		using ushort = unsigned short;

		constexpr bool leap_year(const size_t year) noexcept {
			return year % 400 == 0 || (!(year % 100 == 0) && year % 4 == 0);
		}
		constexpr bool check(const llong year, const ushort month, const ushort day)
		{
			constexpr ushort ord_days[12]{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
			constexpr ushort leap_days[12]{ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
			if (day == 0) return false;
			if (leap_year(year)) return day <= leap_days[month - 1];
			else return day <= ord_days[month - 1];
		}

		class datetime;

		/// <summary>
		/// object contains year, month and day
		/// </summary>
		class date {
			friend datetime;
		private:
			llong _year{ 1 };
			ushort _month{ 1 }, _day{ 1 };

		public:
			constexpr date() noexcept = default;
			constexpr date(const llong year, const ushort month, const ushort day) {
				if (month < 1 || month > 12) throw std::invalid_argument("Invalid month!");
				if (!check(year, month, day)) throw std::invalid_argument("Invalid day!");
				_year = year;
				_month = month;
				_day = day;
			}
			date(const date& other) = default;
			date(date&& other) noexcept = default;
			~date() noexcept = default;

			date& operator = (const date& other) = default;
			date& operator = (date&& other) noexcept = default;

			constexpr llong get_year() const noexcept { return _year; }
			constexpr ushort get_month() const noexcept { return _month; }
			constexpr ushort get_day() const noexcept { return _day; }

			friend constexpr bool operator > (const date& f, const date& s) noexcept {
				if (f._year > s._year) return true;
				else if (f._year == s._year) {
					if (f._month > s._month) return true;
					else if (f._month == s._month) {
						return f._day > s._day;
					}
				}
				return false;
			}
			friend constexpr bool operator < (const date& f, const date& s) noexcept {
				if (f._year < s._year) return true;
				else if (f._year == s._year) {
					if (f._month < s._month) return true;
					else if (f._month == s._month) {
						return f._day < s._day;
					}
				}
				return false;
			}
			friend constexpr bool operator == (const date& f, const date& s) noexcept {
				return f._year == s._year && f._month == s._month && f._day == s._day;
			}
			friend constexpr bool operator >= (const date& f, const date& s) noexcept {
				return f > s || f == s;
			}
			friend constexpr bool operator <= (const date& f, const date& s) noexcept {
				return f < s || f == s;
			}

			friend std::ostream& operator << (std::ostream& o, const date& d);
			friend std::istream& operator >> (std::istream& i, date& d);
		};

		/// <summary>
		/// object contains hour, minute, second and millisecond
		/// </summary>
		class time {
			friend datetime;
		private:
			ushort _hour{}, _minute{}, _second{}, _millisec{};

		public:
			constexpr time() noexcept = default;
			constexpr time(const ushort hour, const ushort minute, const ushort second, const ushort millisecond = 0) {
				if (hour > 23) throw std::invalid_argument("Invalid hour > 23!");
				if (minute > 59) throw std::invalid_argument("Invalid minute > 59!");
				if (second > 59) throw std::invalid_argument("Invalid second > 59!");
				if (millisecond > 999) throw std::invalid_argument("Invalid millisecond > 999!");
				_hour = hour;
				_minute = minute;
				_second = second;
				_millisec = millisecond;
			}
			time(const time& t) = default;
			time(time&& t) noexcept = default;
			~time() noexcept = default;

			time& operator = (const time& t) = default;
			time& operator = (time&& t) noexcept = default;

			constexpr unsigned short get_hour() const noexcept { return _hour; }
			constexpr unsigned short get_minute() const noexcept { return _minute; }
			constexpr unsigned short get_second() const noexcept { return _second; }
			constexpr unsigned short get_millisecond() const noexcept { return _millisec; }

			friend constexpr bool operator > (const time& f, const time& s) noexcept {
				if (f._hour > s._hour) return true;
				else if (f._hour == s._hour) {
					if (f._minute > s._minute) return true;
					else if (f._minute == s._minute) {
						if (f._second > s._second) return true;
						else if (f._second == s._second) {
							return f._millisec > s._millisec;
						}
					}
				}
				return false;
			}
			friend constexpr bool operator < (const time& f, const time& s) noexcept {
				if (f._hour < s._hour) return true;
				else if (f._hour == s._hour) {
					if (f._minute < s._minute) return true;
					else if (f._minute == s._minute) {
						if (f._second < s._second) return true;
						else if (f._second == s._second) {
							return f._millisec < s._millisec;
						}
					}
				}
				return false;
			}
			friend constexpr bool operator == (const time& f, const time& s) noexcept {
				return f._hour == s._hour && f._minute == s._minute && f._second == s._second && f._millisec == s._millisec;
			}
			friend constexpr bool operator >= (const time& f, const time& s) noexcept {
				return f > s || f == s;
			}
			friend constexpr bool operator <= (const time& f, const time& s) noexcept {
				return f < s || f == s;
			}

			friend std::ostream& operator << (std::ostream& o, const time& d);
		};

		/// <summary>
		/// object implements both date and time of the day
		/// </summary>
		class datetime
		{
		private:
			date _date{};
			time _time{};

		public:
			constexpr datetime() noexcept = default;
			constexpr datetime(const date& d, const time& t) : _date{ d }, _time{ t } {}
			constexpr datetime(
				const llong year, const ushort month, const ushort day, 
				const ushort hour, const ushort minute, const ushort second, const ushort millisecond = 0) :
				_date{ date(year, month, day) }, 
				_time{ time(hour, minute, second, millisecond) }
			{}
			datetime(const datetime& dt) = default;
			datetime(datetime&& dt) noexcept = default;
			~datetime() noexcept = default;

			constexpr const date& get_date() const { return _date; }
			constexpr const time& get_time() const { return _time; }

			constexpr llong get_year() const { return _date._year; }
			constexpr ushort get_month() const { return _date._month; }
			constexpr ushort get_day() const { return _date._day; }
			constexpr ushort get_hour() const { return _time._hour; }
			constexpr ushort get_minute() const { return _time._minute; }
			constexpr ushort get_second() const { return _time._second; }
			constexpr ushort get_millisecond() const { return _time._millisec; }

			datetime& operator = (const datetime& dt) = default;
			datetime& operator = (datetime&& dt) noexcept = default;;

			friend constexpr bool operator > (const datetime& f, const datetime& s) noexcept {
				return f._date > s._date || (f._date == s._date && f._time > s._time);
			}
			friend constexpr bool operator < (const datetime& f, const datetime& s) noexcept {
				return f._date < s._date || (f._date == s._date && f._time < s._time);
			}
			friend constexpr bool operator == (const datetime& f, const datetime& s) noexcept {
				return f._date == s._date && f._time == s._time;
			}
			friend constexpr bool operator >= (const datetime& f, const datetime& s) noexcept {
				return f == s || f > s;
			}
			friend constexpr bool operator <= (const datetime& f, const datetime& s) noexcept {
				return f == s || f < s;
			}

			friend std::ostream& operator << (std::ostream& o, const datetime& d);

			static datetime now();
		};

		/// <summary>
		/// Restore datetime from string of certain format
		/// </summary>
		/// <param name="str"> - a string representation of datetime</param>
		/// <param name="format"> - a string format of datetime representation (y - year, M - month, d - day, h - hour, M - minute, s - second, f - millisecond)</param>
		/// <returns>a datetime corresponding the string (throws an exception when invalid format)</returns>
		datetime datetime_from_str(const std::string& str, const std::string& format = "y/M/d h:m:s.f");

		/// <summary>
		/// Class implements julian date refered to midnight count down
		/// </summary>
		class JD
		{
		private:
			llong _day{ 1 };
			double _time{};

			constexpr void shift_behind() noexcept { _time += 1.0; _day -= 1; }

		public:
			constexpr JD() noexcept = default;
			explicit constexpr JD(const datetime& dt) noexcept {
				llong 
					a = (14 - dt.get_month()) / 12,
					y = 4800 + dt.get_year() - a,
					m = dt.get_month() + 12 * a - 3;
				_day = dt.get_day() + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
				_time = dt.get_hour() * (1 / 24.0) + dt.get_minute() * (1 / 1440.0) + (dt.get_second() + dt.get_millisecond() * 1e-3) * (1 / 86400.0);
			}
			constexpr JD(const llong day, const double time) {
				if (time >= 0.0 && time < 1.0) {
					_day = day;
					_time = time;
				}
				else throw std::invalid_argument("invalid time of the day!");
			}
			constexpr JD(const double jd) noexcept {
				_day = static_cast<llong>(jd);
				_time = jd - _day;
			}
			JD(const JD& jd) = default;
			JD(JD&& jd) noexcept = default;
			~JD() = default;

			JD& operator = (const double jd) noexcept;
			JD& operator = (const JD& jd) = default;
			JD& operator = (JD && jd) noexcept = default;
			
			/// <summary>
			/// day number
			/// </summary>
			/// <returns>day number</returns>
			constexpr llong JDN() const noexcept { return _day; }
			/// <summary>
			/// part of day
			/// </summary>
			/// <returns>value between 0 and 1</returns>
			constexpr double T() const noexcept { return _time; }
			/// <summary>
			/// representing as a single number (an integer part is day number + fractional part is a part of the day)
			/// </summary>
			/// <returns>single value</returns>
			constexpr double to_double() const noexcept { return _day + _time; }

			/// <summary>
			/// addition period
			/// </summary>
			/// <param name="dt">is time period in seconds</param>
			/// <returns>instance of JD</returns>
			JD& operator += (const double dt) noexcept;
			/// <summary>
			/// subtracting period
			/// </summary>
			/// <param name="dt">is time period in seconds</param>
			/// <returns>instance of JD</returns>
			JD& operator -= (const double dt) noexcept;
			
			JD& add_days(const int n) noexcept;
			JD& add_hours(const int n) noexcept;
			JD& add_minutes(const int n) noexcept;
			JD& add_seconds(const int n) noexcept;
			/// <summary>
			/// addition period
			/// </summary>
			/// <param name="jd">is a juliand date refered to midnight</param>
			/// <param name="dt">is a period in seconds</param>
			/// <returns>julian date</returns>
			friend constexpr JD operator + (const JD& jd, const double dt) noexcept {
				JD result;
				double t = jd._time + dt / SEC_PER_DAY;
				result._day = static_cast<llong>(t);
				result._time = t - result._day;
				result._day += jd._day;
				return result;
			}
			/// <summary>
			/// subtraction period
			/// </summary>
			/// <param name="jd">is a juliand date refered to midnight</param>
			/// <param name="dt">is a period in seconds</param>
			/// <returns>julian date</returns>
			friend constexpr JD operator - (const JD& jd, const double dt) noexcept {
				JD result;
				double t = jd._time - dt / SEC_PER_DAY;
				result._day = static_cast<llong>(t);
				result._time = t - result._day;
				result._day += jd._day;
				if (result._time < 0.0) result.shift_behind();
				return result;
			}
			/// <summary>
			/// subtracting the julian dates
			/// </summary>
			/// <param name="f">is a left julian date</param>
			/// <param name="s">is a right julian date</param>
			/// <returns>a period between dates in seconds</returns>
			friend constexpr double operator - (const JD& f, const JD& s) noexcept {
				return ((f._time - s._time) + (f._day - s._day)) * general::time::SEC_PER_DAY;
			}

			friend constexpr bool operator < (const JD& f, const JD& s) noexcept {
				return f._day < s._day || (f._day == s._day && f._time < s._time);
			}
			friend constexpr bool operator > (const JD& f, const JD& s) noexcept {
				return f._day > s._day || (f._day == s._day && f._time > s._time);
			}
			friend constexpr bool operator == (const JD& f, const JD& s) noexcept {
				return f._day == s._day && f._time == s._time;
			}
			friend constexpr bool operator <= (const JD& f, const JD& s) noexcept {
				return f == s || f < s;
			}
			friend constexpr bool operator >= (const JD& f, const JD& s) noexcept {
				return f == s || f > s;
			}
			/// <summary>
			/// datetime to julian date conversion
			/// </summary>
			/// <param name="jd">is a julian date refered to midnight</param>
			/// <returns>corresponding datetime</returns>
			friend constexpr datetime jd_to_datetime(const JD& jd) {
				llong
					a = jd._day + 32044,
					b = (4 * a + 3) / 146097,
					c = a - (146097 * b) / 4,
					d = (4 * c + 3) / 1461,
					e = c - (1461 * d) / 4,
					m = (5 * e + 2) / 153,
					year = 100 * b + d - 4800 + (m / 10),
					t = math::roundval(jd._time * SEC_PER_DAY * 1e3);
				ushort
					day = static_cast<ushort>(e - (153 * m + 2) / 5 + 1),
					month = static_cast<ushort>(m + 3 - 12 * (m / 10));
				auto millisec = static_cast<ushort>(t % 1000);
				t /= 1000;
				auto second = static_cast<ushort>(t % 60);
				t /= 60;
				auto minute = static_cast<ushort>(t % 60);
				t /= 60;
				auto hour = static_cast<ushort>(t % 24);
				return datetime(year, month, day, hour, minute, second, millisec);
			}

			friend std::ostream& operator << (std::ostream& o, const JD& jd);
		};
		/// <summary>
		/// A period of time calculator
		/// </summary>
		class stopwatch
		{
			std::chrono::steady_clock::time_point _start, _finish;
		public:
			/// <summary>
			/// Start to calculate the period of time
			/// </summary>
			void start();
			/// <summary>
			/// Finish to calcaulate the period of time
			/// </summary>
			void finish();
			/// <summary>
			/// Calculated period of time
			/// </summary>
			/// <returns>seconds of period duration</returns>
			double duration() const;
		};
		/// <summary>
		/// conversion from decimal time fromat in hours to time structure
		/// </summary>
		/// <param name="hour">is time in hours (decimal fromat)</param>
		/// <returns>(hour, minute, second, millisecond)</returns>
		inline constexpr time decimaltime_to_hms(double hour) noexcept {
			const auto h = static_cast<ushort>(hour);
			hour -= h;
			hour *= MIN_PER_HOUR;
			const auto m = static_cast<ushort>(hour);
			hour -= m;
			hour *= SEC_PER_MIN;
			const auto s = static_cast<ushort>(hour);
			hour -= s;
			hour *= MILLISEC_PER_SEC;
			const auto ms = static_cast<ushort>(hour);
			return time{ h, m, s, ms };
		}
	}
}