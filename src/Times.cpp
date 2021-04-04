#include "Times.h"
#include "GeneralConstants.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <chrono>


namespace general
{
	namespace time
	{
		bool Date::Check(
			const size_t year,
			const unsigned short month,
			const unsigned short day)
		{
			const unsigned short daysOrder1[12]
			{
				31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
			};
			const unsigned short daysOrder2[12]
			{
				31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
			};
			const bool leap = year % 400 == 0 ? true :
				(year % 100 == 0 ? false :
					(year % 4 == 0 ? true : false));
			if (day > 0)
			{
				if (leap) return day <= daysOrder2[month - 1];
				else return day <= daysOrder1[month - 1];
			}
			return false;
		}

		Date::Date(
			const size_t year,
			const unsigned short month,
			const unsigned short day)
		{
			if (month < 1 || month > 12)
				throw std::invalid_argument("Invalid month!");
			if (!Check(year, month, day))
				throw std::invalid_argument("Invalid day!");
			_year = year;
			_month = month;
			_day = day;
		}

		Date::Date(Date&& d) noexcept
		{
			_day = d._day;
			_month = d._month;
			_year = d._year;
			d._year = d._day = d._month = 0;
		}

		Date& Date::operator=(Date&& d) noexcept
		{
			_day = d._day;
			_month = d._month;
			_year = d._year;
			d._year = d._day = d._month = 0;
			return *this;
		}

		bool operator > (const Date& f, const Date& s)
		{
			if (f._year > s._year) return true;
			else if (f._year == s._year) {
				if (f._month > s._month) return true;
				else if (f._month == s._month) {
					return f._day > s._day;
				}
			}
			return false;
		}
		bool operator < (const Date& f, const Date& s)
		{
			if (f._year < s._year) return true;
			else if (f._year == s._year) {
				if (f._month < s._month) return true;
				else if (f._month == s._month) {
					return f._day < s._day;
				}
			}
			return false;
		}
		bool operator == (const Date& f, const Date& s)
		{
			return f._year == s._year && f._month == s._month && f._day == s._day;
		}

		std::ostream& operator << (std::ostream& o, const Date& d)
		{
			o << d._year << '/' << d._month << '/' << d._day;
			return o;
		}
		std::istream& operator >> (std::istream& i, Date& d)
		{
			i >> d._year >> d._month >> d._day;
			return i;
		}


		Time::Time(
			const unsigned short hour,
			const unsigned short minute,
			const unsigned short second,
			const unsigned short millisecond)
		{
			if (hour > 23)
				throw std::invalid_argument("Invalid hour > 23!");
			if (minute > 59)
				throw std::invalid_argument("Invalid minute > 59!");
			if (second > 59)
				throw std::invalid_argument("Invalid second > 59!");
			if (millisecond > 999)
				throw std::invalid_argument("Invalid millisecond > 999!");
			_hour = hour;
			_minute = minute;
			_second = second;
			_millisec = millisecond;
		}

		Time::Time(Time&& t) noexcept
		{
			_hour = t._hour;
			_minute = t._minute;
			_second = t._second;
			_millisec = t._millisec;
			t._hour = t._millisec = t._minute = t._second = 0;
		}

		Time& Time::operator=(Time&& t) noexcept
		{
			_hour = t._hour;
			_minute = t._minute;
			_second = t._second;
			_millisec = t._millisec;
			t._hour = t._millisec = t._minute = t._second = 0;
			return *this;
		}

		std::ostream& operator << (std::ostream& o, const Time& t)
		{
			o << t._hour << ':' << t._minute << ':' << t._second + t._millisec * 1e-3;
			return o;
		}

		bool operator > (const Time& f, const Time& s)
		{
			if (f._hour > s._hour) return true;
			else if (f._hour == s._hour)
			{
				if (f._minute > s._minute) return true;
				else if (f._minute == s._minute)
				{
					if (f._second > s._second) return true;
					else if (f._second == s._second)
					{
						return f._millisec > s._millisec;
					}
				}
			}
			return false;
		}
		bool operator < (const Time& f, const Time& s)
		{
			if (f._hour < s._hour) return true;
			else if (f._hour == s._hour)
			{
				if (f._minute < s._minute) return true;
				else if (f._minute == s._minute)
				{
					if (f._second < s._second) return true;
					else if (f._second == s._second)
					{
						return f._millisec < s._millisec;
					}
				}
			}
			return false;
		}
		bool operator == (const Time& f, const Time& s)
		{
			return f._hour == s._hour && f._minute == s._minute && f._second == s._second && f._millisec == s._millisec;
		}


		std::ostream& operator << (std::ostream& o, const DateTime& d)
		{
			o << d._date << ' ' << d._time;
			return o;
		}
		bool operator > (const DateTime& f, const DateTime& s)
		{
			if (f._date > s._date) return true;
			if (f._date == s._date)
				return f._time > s._time;
			return false;
		}
		bool operator < (const DateTime& f, const DateTime& s)
		{
			if (f._date < s._date) return true;
			else if (f._date == s._date)
				return f._time < s._time;
			return false;
		}
		bool operator == (const DateTime& f, const DateTime& s)
		{
			return f._date == s._date && f._time == s._time;
		}

		bool try_parse(const std::string& str, DateTime& dt, const std::string& format)
		{
			char separators[6]{ 0 };
			auto get_index = [](const char symbol) {
				switch (symbol) {
					case 'y': return 0;
					case 'M': return 1;
					case 'd': return 2;
					case 'h': return 3;
					case 'm': return 4;
					case 's': return 5;
					case 'f': return 6;
				}
				return 7;
			};
			size_t data[7]{ 1, 1, 1, 0, 0, 0 }, indices[7]{ 0 };
			size_t index{ 0 };
			size_t value;
			for (const char& c : format) {
				value = get_index(c);
				if (value == 7) separators[index++] = c;
				else indices[index] = value;
			}
			if (index < 3) return false;
			std::string buf;
			index = 0;
			for (const char& c : str) {
				if (c != separators[index]) buf.push_back(c);
				else {
					if (indices[index] == 6) value = buf.size();
					data[indices[index++]] = std::atoi(buf.c_str());
					buf.clear();
				}
			}
			if (buf.size() != 0) data[indices[index]] = std::atoi(buf.c_str()); 
			if (indices[index] == 6) value = buf.size();
			try {
				dt = DateTime(
					data[0],
					static_cast<unsigned short>(data[1]),
					static_cast<unsigned short>(data[2]),
					static_cast<unsigned short>(data[3]),
					static_cast<unsigned short>(data[4]),
					static_cast<unsigned short>(data[5]),
					static_cast<unsigned short>(data[6] * std::pow(10, 3 - static_cast<double>(value)))
				);
			} catch (std::exception) { return false; }
			return true;
		}

		DateTime DateTime::now()
		{
			auto time{ std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) };
			auto t = std::localtime(&time);
			return DateTime(
				size_t(1900) + static_cast<size_t>(t->tm_year), 
				static_cast<unsigned short>(t->tm_mon + 1), 
				static_cast<unsigned short>(t->tm_mday), 
				static_cast<unsigned short>(t->tm_hour), 
				static_cast<unsigned short>(t->tm_min), 
				static_cast<unsigned short>(t->tm_sec));
		}

		JD::JD(const double jd)
		{
			if (jd < 0) throw std::invalid_argument("Invalid julian date < 0!");
			double jdn;
			_time = std::modf(jd, &jdn);
			_day = static_cast<long_t>(jdn);
		}
		JD::JD(const DateTime& dt)
		{
			Date d = dt.get_date();
			Time t = dt.get_time();
			size_t a = (14 - d.get_month()) / 12,
				y = 4800 + d.get_year() - a,
				m = d.get_month() + 12 * a - 3;
			_day = d.get_day() + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
			_time = (t.get_hour() + (t.get_minute() + (t.get_second() + t.get_millisecond() * 1e-3) / 60) / 60) / 24;
		}

		JD::JD(const long_t day, const double time)
		{
			if (time > 0 && time < 1.0) {
				_day = static_cast<long_t>(day);
				_time = time;
			}
			else throw std::invalid_argument("Time is invalid (should be 0 < t < 1)!");
		}

		DateTime JD::to_datetime() const
		{
			size_t
				a = static_cast<size_t>(_day) + 32044,
				b = (4 * a + 3) / 146097,
				c = a - (146097 * b) / 4,
				d = (4 * c + 3) / 1461,
				e = c - (1461 * d) / 4,
				m = (5 * e + 2) / 153,
				year = 100 * b + d - 4800 + (m / 10),
				t = static_cast<size_t>(_time * SEC_PER_DAY * 1e3);
			unsigned short 
				day = static_cast<unsigned short>(e - (153 * m + 2) / 5 + 1),
				month = static_cast<unsigned short>(m + 3 - 12 * (m / 10));
			auto millisec = static_cast<unsigned short>(t % 1000);
			t /= 1000;
			auto second = static_cast<unsigned short>(t % 60);
			t /= 60;
			auto minute = static_cast<unsigned short>(t % 60);
			t /= 60;
			auto hour = static_cast<unsigned short>(t % 24);
			return DateTime(year, month, day, hour, minute, second, millisec);
		}

		double JD::to_double() const
		{
			return _time + _day;
		}

		JD& JD::operator = (const double jd)
		{
			double jdn;
			_time = std::modf(jd, &jdn);
			_day = static_cast<long_t>(jdn);
			return *this;
		}
		JD::JD(JD&& jd) noexcept : _day{ jd._day }, _time{ jd._time }
		{
			jd._time = 0; jd._day = 1;
		}

		JD& JD::operator = (JD&& jd) noexcept
		{
			_day = jd._day;
			_time = jd._time;
			jd._time = 0; jd._day = 1;
			return *this;
		}

		JD& JD::operator += (const double dt)
		{
			double day;
			_time = std::modf(_time + dt, &day);
			_day += static_cast<long_t>(day);
			return *this;
		}

		JD& JD::operator -= (const double dt)
		{
			double day;
			_time = std::modf(_time - dt, &day);
			_day += static_cast<long_t>(day);
			return *this;
		}

		JD& JD::add_days(const int n)
		{
			_day += n;
			return *this;
		}
		JD& JD::add_hours(const int n)
		{
			double day;
			_time = std::modf(_time + n / HOURS_PER_DAY, &day);
			_day += static_cast<long_t>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}
		JD& JD::add_minutes(const int n)
		{
			double day;
			_time = std::modf(_time + n / MIN_PER_DAY, &day);
			_day += static_cast<long_t>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}
		JD& JD::add_seconds(const int n)
		{
			double day;
			_time = std::modf(_time + n / SEC_PER_DAY, &day);
			_day += static_cast<long_t>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}

		JD operator + (const JD& jd, const double dt)
		{
			auto result{ jd };
			double day;
			result._time = std::modf(jd._time + dt, &day);
			result._day += static_cast<long long>(day);
			return result;
		}
		JD operator - (const JD& jd, const double dt)
		{
			auto result{ jd };
			double day;
			result._time = std::modf(jd._time - dt, &day);
			result._day += static_cast<long long>(day);
			return result;
		}

		double operator - (const JD& f, const JD& s)
		{
			using long_t = long long;
			return (f._time - s._time) + (f._day - s._day);
		}

		bool operator < (const JD& f, const JD& s)
		{
			return f._day < s._day || (f._day == s._day && f._time < s._time);
		}
		bool operator > (const JD& f, const JD& s)
		{
			return f._day > s._day || (f._day == s._day && f._time > s._time);
		}
		bool operator == (const JD& f, const JD& s)
		{
			return f._day == s._day && std::abs(f._time - s._time) < (1e-3 / SEC_PER_DAY);
		}
		bool operator <= (const JD& f, const JD& s)
		{
			return f._day < s._day || (f._day == s._day && f._time <= s._time);
		}
		bool operator >= (const JD& f, const JD& s)
		{
			return f._day > s._day || (f._day == s._day && f._time >= s._time);
		}

		std::ostream& operator << (std::ostream& o, const JD& jd)
		{
			o << jd._day;
			if (jd._time > 1e-10) {
				std::ostringstream str;
				str << std::setprecision(10) << jd._time;
				auto val{ str.str() };
				o << ".";
				for (size_t i = 2; i < val.size(); ++i)
					o << val[i];
			}
			return o;
		}
		DateTime::DateTime(DateTime&& dt) noexcept
		{
			_date = std::move(dt._date);
			_time = std::move(dt._time);
		}
		DateTime& DateTime::operator=(DateTime&& dt) noexcept
		{
			_date = std::move(dt._date);
			_time = std::move(dt._time);
			return *this;
		}
}
}