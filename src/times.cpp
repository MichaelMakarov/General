#include "times.h"
#include <iomanip>
#include <cmath>

namespace general
{
	namespace time
	{
		const char zero = '0';

		std::ostream& operator << (std::ostream& o, const date& d)
		{
			return o << 
				std::setfill(zero) << std::setw(4) << d._year << '/' << 
				std::setfill(zero) << std::setw(2) << d._month << '/' << 
				std::setfill(zero) << std::setw(2) << d._day;
		}
		std::istream& operator >> (std::istream& i, date& d)
		{
			return i >> d._year >> d._month >> d._day;
		}


		std::ostream& operator << (std::ostream& o, const time& t)
		{
			double integ, decim;
			decim = std::modf(std::round(t._second * 1e3 + t._millisec) * 1e-3, &integ);
			return o << 
				std::setfill(zero) << std::setw(2) << t._hour << ':' <<
				std::setfill(zero) << std::setw(2) << t._minute << ':' <<
				std::setfill(zero) << std::setw(2) << integ << '.' <<
				std::setfill(zero) << std::setw(3) << decim * 1e3;
		}

		std::ostream& operator << (std::ostream& o, const datetime& d)
		{
			return o << d._date << ' ' << d._time;
		}

		datetime datetime_from_str(const std::string& str, const std::string& format)
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
			if (index < 3) throw std::invalid_argument(str + " does not correspond format " + format);
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
			return datetime(
				data[0],
				static_cast<ushort>(data[1]),
				static_cast<ushort>(data[2]),
				static_cast<ushort>(data[3]),
				static_cast<ushort>(data[4]),
				static_cast<ushort>(data[5]),
				static_cast<ushort>(data[6] * std::pow(10, 3 - static_cast<double>(value)))
			);
		}

		datetime datetime::now()
		{
			auto time{ std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) };
			auto t = std::localtime(&time);
			return datetime(
				size_t(1900) + static_cast<size_t>(t->tm_year), 
				static_cast<unsigned short>(t->tm_mon + 1), 
				static_cast<unsigned short>(t->tm_mday), 
				static_cast<unsigned short>(t->tm_hour), 
				static_cast<unsigned short>(t->tm_min), 
				static_cast<unsigned short>(t->tm_sec));
		}

		JD& JD::operator = (const double jd) noexcept
		{
			double jdn;
			_time = std::modf(jd, &jdn);
			_day = static_cast<llong>(jdn);
			return *this;
		}

		JD& JD::operator += (const double dt) noexcept
		{
			double day;
			_time = std::modf(_time + dt / SEC_PER_DAY, &day);
			_day += static_cast<llong>(day);
			return *this;
		}

		JD& JD::operator -= (const double dt) noexcept
		{
			double day;
			_time = std::modf(_time - dt / SEC_PER_DAY, &day);
			_day += static_cast<llong>(day);
			return *this;
		}

		JD& JD::add_days(const int n) noexcept
		{
			_day += n;
			return *this;
		}
		JD& JD::add_hours(const int n) noexcept
		{
			double day;
			_time = std::modf(_time + n / HOURS_PER_DAY, &day);
			_day += static_cast<llong>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}
		JD& JD::add_minutes(const int n) noexcept
		{
			double day;
			_time = std::modf(_time + n / MIN_PER_DAY, &day);
			_day += static_cast<llong>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}
		JD& JD::add_seconds(const int n) noexcept
		{
			double day;
			_time = std::modf(_time + n / SEC_PER_DAY, &day);
			_day += static_cast<llong>(day);
			if (_time < 0.0) shift_behind();
			return *this;
		}


		std::ostream& operator << (std::ostream& ostr, const JD& jd)
		{
			ostr << jd._day;
			if (jd._time > 1e-9) {
				ostr << std::setprecision(9) << '.' << std::round(jd._time * 1e9);
			}
			return ostr;
		}

		void stopwatch::start()
		{
			_finish = _start = std::chrono::high_resolution_clock::now();
		}
		void stopwatch::finish()
		{
			_finish = std::chrono::high_resolution_clock::now();
		}
		double stopwatch::duration() const
		{
			return std::chrono::duration<double, std::nano>(_finish - _start).count() * 1e-9;
		}
	}
}