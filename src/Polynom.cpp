#include "Polynom.h"
#include "Matrix.h"

namespace general
{
	namespace math
	{
		Polynom create_polynom(const VectorDyn& X, const VectorDyn& Y, const size_t degree)
		{
			if (X.size() != Y.size())
				throw std::invalid_argument("X and Y have different dimensions!");
			auto A{ MatrixDyn(X.size(), degree + 1) };
			auto polynom{ Polynom(degree) };
			for (size_t i = 0; i < X.size(); ++i) {
				A(i, 0) = 1.0;
				for (size_t k = 1; k < degree + 1; ++k)
					A(i, k) = A(i, k - 1) * X[i];
			}
			auto AT = transpose(A);
			auto M{ AT * A };
			auto D{ MatrixDyn(M.rows(), M.columns()) };
			for (size_t i = 0; i < M.rows(); ++i) D(i, i) = 1 / std::sqrt(M(i, i));
			auto C = DxA(D, AxD(inverse(DxA(D, AxD(M, D))), D)) * AT * Y;
			for (size_t i = 0; i < degree + 1; ++i)
				polynom[i] = C[i];
			return polynom;
		}
	}
}