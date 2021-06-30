#pragma once

namespace general
{
    namespace time
    {
        inline constexpr double SEC_PER_MIN{ 60 };
        inline constexpr double SEC_PER_HOUR{ 3600 };
        inline constexpr double SEC_PER_DAY{ 86400 };

        inline constexpr double MILLISEC_PER_HOUR{ 3600000 };
        inline constexpr double MILLISEC_PER_MIN{ 60000 };
        inline constexpr double MILLISEC_PER_SEC{ 1000 };

        inline constexpr double HOURS_PER_DAY{ 24 };

        inline constexpr double MIN_PER_HOUR{ 60 };
        inline constexpr double MIN_PER_DAY{ 1440 };

        inline constexpr double DAYS_PER_YEAR{ 365 };

        // Julian date for 1st January 1900, 12 am
        inline constexpr double JD1900{ 2415021.5 };
        // Julian date for 1st January 2000, 12 am
        inline constexpr double JD2000{ 2451545.5 };
        // Julian date for 30th December 1899, 12 am
        inline constexpr double JD1899{ 2415019.5 };
    }

    namespace math
    {
        inline constexpr double PI{ 3.1415926535897932 };
	    inline constexpr double PI1_2{ PI * 0.5 };
	    inline constexpr double PI3_2{ PI * 1.5 };
	    inline constexpr double PI2{ PI * 2 };
        // seconds per 360 degrees
        inline constexpr double SEC_PER_ROUND{ 1296000 };
        // radians per angle second
        inline constexpr double RAD_PER_SEC{ PI2 / SEC_PER_ROUND };
    }
}