#pragma once
#include <ostream>
#include <istream>
#include <string>
#include <chrono>

namespace general
{
	namespace time
	{
		class DateTime;

		/// <summary>
		/// Class implements date functionality and contains year, month and day
		/// </summary>
		class Date
		{
			friend DateTime;
		private:
			size_t _year;
			unsigned short _month;
			unsigned short _day;

			bool Check(
				const size_t year, 
				const unsigned short month, 
				const unsigned short day);

		public:
			Date() : _year{ 0 }, _month{ 0 }, _day{ 0 } {}
			Date(
				const size_t year,
				const unsigned short month,
				const unsigned short day);
			Date(const Date& d) noexcept = default;
			Date(Date&& d) noexcept;
			~Date() noexcept = default;

			Date& operator = (const Date& d) noexcept = default;
			Date& operator = (Date&& d) noexcept;

			size_t get_year() const noexcept { return _year; }
			unsigned short get_month() const noexcept { return _month; }
			unsigned short get_day() const noexcept { return _day; }

			friend bool operator > (const Date& f, const Date& s) noexcept;
			friend bool operator < (const Date& f, const Date& s) noexcept;
			friend bool operator == (const Date& f, const Date& s) noexcept;

			friend std::ostream& operator << (std::ostream& o, const Date& d);
			friend std::istream& operator >> (std::istream& i, Date& d);
		};

		/// <summary>
		/// Class implements time of the day and contains hour, minute, second and millisecond
		/// </summary>
		class Time
		{
			friend DateTime;
		private:
			unsigned short _hour, _minute, _second, _millisec;

		public:
			Time() : _hour{ 0 }, _minute{ 0 }, _second{ 0 }, _millisec{ 0 } {}
			Time(
				const unsigned short hour,
				const unsigned short minute,
				const unsigned short second,
				const unsigned short millisecond = 0);
			Time(const Time& t) noexcept = default;
			Time(Time&& t) noexcept;
			~Time() noexcept = default;

			Time& operator = (const Time& t) noexcept = default;
			Time& operator = (Time&& t) noexcept;

			unsigned short get_hour() const noexcept { return _hour; }
			unsigned short get_minute() const noexcept { return _minute; }
			unsigned short get_second() const noexcept { return _second; }
			unsigned short get_millisecond() const noexcept { return _millisec; }

			friend bool operator > (const Time& f, const Time& s) noexcept;
			friend bool operator < (const Time& f, const Time& s) noexcept;
			friend bool operator == (const Time& f, const Time& s) noexcept;

			friend std::ostream& operator << (std::ostream& o, const Time& d);
		};

		/// <summary>
		/// Class implements functionality of date and time of the day
		/// </summary>
		class DateTime
		{
		private:
			Date _date;
			Time _time;

		public:
			DateTime() = default;
			DateTime(const Date& date, const Time& time) : _date{ date }, _time{ time } {}
			DateTime(
				const size_t year,
				const unsigned short month,
				const unsigned short day,
				const unsigned short hour,
				const unsigned short minute,
				const unsigned short second,
				const unsigned short millisecond = 0) :
				_date{ Date(year, month, day) }, 
				_time{ Time(hour, minute, second, millisecond) }
			{}
			DateTime(const DateTime& dt) noexcept = default;
			DateTime(DateTime&& dt) noexcept;
			~DateTime() noexcept = default;

			const Date& get_date() const { return _date; }
			const Time& get_time() const { return _time; }

			size_t get_year() const { return _date._year; }
			unsigned short get_month() const { return _date._month; }
			unsigned short get_day() const { return _date._day; }
			unsigned short get_hour() const { return _time._hour; }
			unsigned short get_minute() const { return _time._minute; }
			unsigned short get_second() const { return _time._second; }
			unsigned short get_millisecond() const { return _time._millisec; }

			DateTime& operator = (const DateTime& dt) noexcept = default;
			DateTime& operator = (DateTime&& dt) noexcept;

			friend bool operator > (const DateTime& f, const DateTime& s);
			friend bool operator < (const DateTime& f, const DateTime& s);
			friend bool operator == (const DateTime& f, const DateTime& s);

			friend std::ostream& operator << (std::ostream& o, const DateTime& d);

			static DateTime now();
		};

		/// <summary>
		/// Restore datetime from string of certain format
		/// </summary>
		/// <param name="str"> - a string representation of datetime</param>
		/// <param name="format"> - a string format of datetime representation (y - year, M - month, d - day, h - hour, M - minute, s - second, f - millisecond)</param>
		/// <returns>a datetime corresponding the string (throws an exception when invalid format)</returns>
		DateTime datetime_from_str(const std::string& str, const std::string& format = "y/M/d h:m:s.f");

		/// <summary>
		/// Class implements julian date refered to midnight count down
		/// </summary>
		class JD
		{
			using long_t = long long;
		private:
			long_t _day;
			double _time;

			void shift_behind() noexcept { _time += 1.0; _day -= 1; }

		public:
			JD() : _day{ 1 }, _time{ 0 } {}
			explicit JD(const DateTime& datetime);
			JD(const long_t days, const double time_of_day);
			JD(const double jd);
			JD(const JD& jd) noexcept = default;
			JD(JD && jd) noexcept;
			~JD() = default;

			JD& operator = (const double jd) noexcept;
			JD& operator = (const JD& jd) noexcept = default;
			JD& operator = (JD && jd) noexcept;
			
			/// <summary>
			/// day number
			/// </summary>
			/// <returns>day number</returns>
			long_t JDN() const noexcept { return _day; }
			/// <summary>
			/// part of day
			/// </summary>
			/// <returns>value between 0 and 1</returns>
			double T() const noexcept { return _time; }
			/// <summary>
			/// conversion to corresponding datetime structure
			/// </summary>
			/// <returns>corresponding datetime</returns>
			DateTime to_datetime() const;

			/// <summary>
			/// representing as a single number (an integer part is day number + fractional part is a part of the day)
			/// </summary>
			/// <returns>single value</returns>
			double to_double() const noexcept;

			/// <summary>
			/// add time as double
			/// </summary>
			/// <param name="dt"> - time in seconds</param>
			/// <returns>current instance of JD</returns>
			JD& operator += (const double dt) noexcept;
			/// <summary>
			/// subtract time as double
			/// </summary>
			/// <param name="dt"> - time in seconds</param>
			/// <returns>current instance of JD</returns>
			JD& operator -= (const double dt) noexcept;

			JD& add_days(const int n) noexcept;
			JD& add_hours(const int n) noexcept;
			JD& add_minutes(const int n) noexcept;
			JD& add_seconds(const int n) noexcept;

			friend JD operator + (const JD& jd, const double dt) noexcept;
			friend JD operator - (const JD& jd, const double dt) noexcept;
			friend double operator - (const JD& f, const JD& s) noexcept;

			friend bool operator < (const JD& f, const JD& s) noexcept;
			friend bool operator > (const JD& f, const JD& s) noexcept;
			friend bool operator == (const JD& f, const JD& s) noexcept;
			friend bool operator <= (const JD& f, const JD& s) noexcept;
			friend bool operator >= (const JD& f, const JD& s) noexcept;

			friend std::ostream& operator << (std::ostream& o, const JD& jd);
		};
		/// <summary>
		/// A period of time calculator
		/// </summary>
		class Stopwatch
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
	}
}