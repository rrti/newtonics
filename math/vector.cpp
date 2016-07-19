#include "./vector.hpp"

namespace newtonics {
	namespace lib_math {
		static const t_vector4f      eps_vector_4f = t_vector4f(M_FEPS, M_FEPS, M_FEPS,  0.0f);
		static const t_vector4f      inf_vector_4f = t_vector4f(M_FINF, M_FINF, M_FINF,  0.0f);
		static const t_vector4f     zero_vector_4f = t_vector4f(  0.0f,   0.0f,   0.0f,  0.0f);
		static const t_vector4f     ones_vector_4f = t_vector4f(  1.0f,   1.0f,   1.0f,  0.0f);
		static const t_vector4f    error_vector_4f = t_vector4f( -1.0f,  -1.0f,  -1.0f,  0.0f);
		static const t_vector4f   x_axis_vector_4f = t_vector4f(  1.0f,   0.0f,   0.0f,  0.0f);
		static const t_vector4f   y_axis_vector_4f = t_vector4f(  0.0f,   1.0f,   0.0f,  0.0f);
		static const t_vector4f   z_axis_vector_4f = t_vector4f(  0.0f,   0.0f,   1.0f,  0.0f);
		static const t_vector4f   w_axis_vector_4f = t_vector4f(  0.0f,   0.0f,   0.0f,  1.0f);
		static const t_vector4f  xz_axis_vector_4f = x_axis_vector_4f + z_axis_vector_4f;
		static const t_vector4f  xy_axis_vector_4f = x_axis_vector_4f + y_axis_vector_4f;
		static const t_vector4f  yz_axis_vector_4f = y_axis_vector_4f + z_axis_vector_4f;
		static const t_vector4f xyz_axis_vector_4f = x_axis_vector_4f + y_axis_vector_4f + z_axis_vector_4f;

		template<> const t_vector4f& t_vector4f::     eps_vector() { return      eps_vector_4f; }
		template<> const t_vector4f& t_vector4f::     inf_vector() { return      inf_vector_4f; }
		template<> const t_vector4f& t_vector4f::    zero_vector() { return     zero_vector_4f; }
		template<> const t_vector4f& t_vector4f::    ones_vector() { return     ones_vector_4f; }
		template<> const t_vector4f& t_vector4f::   error_vector() { return    error_vector_4f; }
		template<> const t_vector4f& t_vector4f::  x_axis_vector() { return   x_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f::  y_axis_vector() { return   y_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f::  z_axis_vector() { return   z_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f::  w_axis_vector() { return   w_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f:: xz_axis_vector() { return  xz_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f:: xy_axis_vector() { return  xy_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f:: yz_axis_vector() { return  yz_axis_vector_4f; }
		template<> const t_vector4f& t_vector4f::xyz_axis_vector() { return xyz_axis_vector_4f; }
	};
};

