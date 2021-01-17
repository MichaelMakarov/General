#include "GeneralConstants.h"
#include <vector>
#include <ostream>
#include <istream>

namespace general
{
	namespace math
	{
        inline double deg_to_rad(const double degrees)
        {
            return degrees * PI / 180.0;
        }
        inline double rad_to_deg(const double radians)
        {
            return radians * 180.0 / PI;
        }

        long double factorial(const size_t x);
	}
}