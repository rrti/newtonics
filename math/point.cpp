#include "./point.hpp"

namespace newtonics {
	namespace lib_math {
		static const t_point4f null_point_4f = t_point4f(  0.0f,   0.0f,   0.0f,   0.0f); // w=0!
		static const t_point4f zero_point_4f = t_point4f(  0.0f,   0.0f,   0.0f,   1.0f);
		static const t_point4f ones_point_4f = t_point4f(  1.0f,   1.0f,   1.0f,   1.0f);
		static const t_point4f  eps_point_4f = t_point4f(M_FEPS, M_FEPS, M_FEPS, M_FEPS);

		template<> const t_point4f& t_point4f::null_point() { return null_point_4f; }
		template<> const t_point4f& t_point4f::zero_point() { return zero_point_4f; }
		template<> const t_point4f& t_point4f::ones_point() { return ones_point_4f; }
		template<> const t_point4f& t_point4f:: eps_point() { return  eps_point_4f; }
	};
};

