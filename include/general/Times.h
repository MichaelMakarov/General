#pragma once
#include <ostream>
#include <string>

namespace general
{
	namespace time
	{
		class DateTime;

		// Class implements the functionality of the date.
		// It contains year, month and day.
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

			size_t get_year() const { return _year; }
			unsigned short get_month() const { return _month; }
			unsigned short get_day() const { return _day; }

			friend bool operator > (const Date& f, const Date& s);
			friend bool operator < (const Date& f, const Date& s);
			friend bool operator == (const Date& f, const Date& s);

			friend std::ostream& operator << (std::ostream& o, const Date& d);
		};

		// Class implements time of the day;
		// Contains hour, minute, second and millisecond.
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

			unsigned short get_hour() const { return _hour; }
			unsigned short get_minute() const { return _minute; }
			unsigned short get_second() const { return _second; }
			unsigned short get_millisecond() const { return _millisec; }

			friend bool operator > (const Time& f, const Time& s);
			friend bool operator < (const Time& f, const Time& s);
			friend bool operator == (const Time& f, const Time& s);

			friend std::ostream& operator << (std::ostream& o, const Time& d);
		};

		// Class implements both date and time of the day.
		class DateTime
		{
		private:
			Date _date;
			Time _time;

		public:
			DateTime() : _date{}, _time{} {}
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
		};

		/// <summary>
		/// Parsing the datetime from string
		/// </summary>
		/// <param name="str"> - a string representation of datetime</param>
		/// <param name="dt"> - a sample of datetime will be parsed to</param>
		/// <param name="format"> - a string format of datetime representation (y - year, M - month, d - day, h - hour, M - minute, s - second, f - millisecond)</param>
		/// <returns>true if succeeded</returns>
		bool try_parse(const std::string& str, DateTime& dt, const std::string& format = "y/M/d h:m:s.f");

		// Represents the julian date refered to midnight count down.
		class JD
		{
		private:
			size_t _day;
			double _time;

		public:
			JD() : _day{ 0 }, _time{ 0 } {}
			explicit JD(const DateTime& datetime);
			JD(const size_t day, const double time);
			JD(const double jd);
			JD(const JD& jd) noexcept = default;
			JD(JD && jd) noexcept;
			~JD() = default;

			JD& operator = (const double jd);
			JD& operator = (const JD& jd) noexcept = default;
			JD& operator = (JD && jd) noexcept;
			
			// A day number
			size_t JDN() const { return _day; }
			// Time of the day as part of the day
			double T() const { return _time; }
			// Conversion to corresponding datetime
			DateTime to_datetime() const;
			// Representing as single number
			double to_double() const;

			JD& operator += (const double dt);
			JD& operator -= (const double dt);

			void add_days(const int n);
			void add_hours(const int n);
			void add_minutes(const int n);
			void add_seconds(const int n);

			friend JD operator + (const JD& jd, const double dt);
			friend JD operator - (const JD& jd, const double dt);
			friend double operator - (const JD& f, const JD& s);

			friend bool operator < (const JD& f, const JD& s);
			friend bool operator > (const JD& f, const JD& s);
			friend bool operator == (const JD& f, const JD& s);
			friend bool operator <= (const JD& f, const JD& s);
			friend bool operator >= (const JD& f, const JD& s);

			friend std::ostream& operator << (std::ostream& o, const JD& jd);
		};
	}
}