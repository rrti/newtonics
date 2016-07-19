#ifndef MATH_VECTOR_HDR
#define MATH_VECTOR_HDR

#include <cmath>
#include <cstdint>

#include <algorithm>
#include <functional>

#include "./math_defs.hpp"
#include "./template_funcs.hpp"

namespace newtonics {
	namespace lib_math {
		template<typename type> struct t_tuple;
		template<typename type> struct t_point;
		template<typename type> struct t_matrix;

		// represents a free vector (magnitude & direction)
		template<typename type> struct t_vector {
		public:
			t_vector<type>(const type _x = 0, const type _y = 0, const type _z = 0, const type _w = 0) {
				x() = _x; y() = _y;
				z() = _z; w() = _w;
			}
			t_vector<type>(const type v[MATH_VECTOR_SIZE]) {
				x() = v[0]; y() = v[1];
				z() = v[2]; w() = v[3];
			}
			t_vector<type>(const t_vector<type>& v) {
				*this = v;
			}

			bool operator == (const t_vector<type>& v) const { return ( equals(v)); }
			bool operator != (const t_vector<type>& v) const { return (!equals(v)); }

			t_vector<type> operator + (const t_vector<type>& v) const { return (t_vector<type>( x() + v.x(),  y() + v.y(),  z() + v.z(),  w() + v.w())); }
			t_vector<type> operator - (const t_vector<type>& v) const { return (t_vector<type>( x() - v.x(),  y() - v.y(),  z() - v.z(),  w() - v.w())); }
			t_vector<type> operator * (const t_vector<type>& v) const { return (t_vector<type>( x() * v.x(),  y() * v.y(),  z() * v.z(),  w() * v.w())); } // UNOFFICIAL
			t_vector<type> operator / (const t_vector<type>& v) const { return (t_vector<type>( x() / v.x(),  y() / v.y(),  z() / v.z(),  w() / v.w())); } // UNOFFICIAL
			t_vector<type> operator - (                       ) const { return (t_vector<type>(-x(),         -y(),         -z()        , -w()        )); }

			t_vector<type>& operator  = (const t_vector<type>& v) { x()  = v.x(); y()  = v.y(); z()  = v.z(); w()  = v.w(); return *this; }
			t_vector<type>& operator += (const t_vector<type>& v) { x() += v.x(); y() += v.y(); z() += v.z(); w() += v.w(); return *this; }
			t_vector<type>& operator -= (const t_vector<type>& v) { x() -= v.x(); y() -= v.y(); z() -= v.z(); w() -= v.w(); return *this; }
			t_vector<type>& operator *= (const t_vector<type>& v) { x() *= v.x(); y() *= v.y(); z() *= v.z(); w() *= v.w(); return *this; } // UNOFFICIAL
			t_vector<type>& operator /= (const t_vector<type>& v) { x() /= v.x(); y() /= v.y(); z() /= v.z(); w() /= v.w(); return *this; } // UNOFFICIAL

			// t_vector<type> operator + (const type s) const { return (t_vector<type>(x() + s, y() + s, z() + s, w() + s)); } // UNOFFICIAL
			// t_vector<type> operator - (const type s) const { return (t_vector<type>(x() - s, y() - s, z() - s, w() - s)); } // UNOFFICIAL
			t_vector<type> operator * (const type s) const { return (t_vector<type>(x() * s, y() * s, z() * s, w() * s)); }
			t_vector<type> operator / (const type s) const { return (t_vector<type>(x() / s, y() / s, z() / s, w() / s)); }

			// t_vector<type>& operator += (const type s) { x() += s; y() += s; z() += s; w() += s; return *this; } // UNOFFICIAL
			// t_vector<type>& operator -= (const type s) { x() -= s; y() -= s; z() -= s; w() -= s; return *this; } // UNOFFICIAL
			t_vector<type>& operator *= (const type s) { x() *= s; y() *= s; z() *= s; w() *= s; return *this; }
			t_vector<type>& operator /= (const type s) { x() /= s; y() /= s; z() /= s; w() /= s; return *this; }

			// calculates v * M = v' (with v a row-vector)
			// same as v' = M * v (with v a column-vector)
			// UNOFFICIAL! (matrices operate on points)
			t_vector<type> operator * (const t_matrix<type>& m) const {
				const type _x = x() * m[0] + y() * m[1] + z() * m[ 2];
				const type _y = x() * m[4] + y() * m[5] + z() * m[ 6];
				const type _z = x() * m[8] + y() * m[9] + z() * m[10];
				return (t_vector<type>(_x, _y, _z, 0));
			}
			// calculates the matrix (this * v^T) where
			// <this> is interpreted as a column-vector
			// UNOFFICIAL! (matrices operate on points)
			t_matrix<type> operator ^ (const t_vector<type>& v) const {
				const t_vector<type> xv = t_vector<type>(x() * v.x(), y() * v.x(), z() * v.x(), 0);
				const t_vector<type> yv = t_vector<type>(x() * v.y(), y() * v.y(), z() * v.y(), 0);
				const t_vector<type> zv = t_vector<type>(x() * v.z(), y() * v.z(), z() * v.z(), 0);
				return (t_matrix<type>(xv, yv, zv, t_vector<type>::w_axis_vector()));
			}

			// cf. point::to_vector
			t_tuple<type> to_tuple() const {                   return (t_tuple<type>(                 *this)); }
			t_point<type> to_point() const { assert(w() == 0); return (t_point<type>::zero_point() + (*this)); }

			// inner product is only defined between vectors and vectors
			// w-components are normally zero so they do not contribute!
			// (so dedicated "*_xyz" convenience wrappers can be omitted)
			type inner_product(const t_vector<type>& v, const t_vector<type>& mask = xyz_axis_vector()) const {
				return ((x() * v.x() * mask.x()) + (y() * v.y() * mask.y()) + (z() * v.z() * mask.z()) + (w() * v.w() * mask.w()));
			}

			type sq_magnit(const t_vector<type>& mask = xyz_axis_vector()) const { return (sq_len(mask)); }
			type    magnit(const t_vector<type>& mask = xyz_axis_vector()) const { return (   len(mask)); }
			type sq_length(const t_vector<type>& mask = xyz_axis_vector()) const { return (sq_len(mask)); }
			type    length(const t_vector<type>& mask = xyz_axis_vector()) const { return (   len(mask)); }
			type sq_len   (const t_vector<type>& mask = xyz_axis_vector()) const { return (inner_product(*this, mask)); }
			type    len   (const t_vector<type>& mask = xyz_axis_vector()) const { return (lib_math::sqrt(sq_len(mask))); }

			#if 0
			type dot_xy(const t_vector<type>& v) const { return inner_product(v, xy_axis_vector()); }
			type dot_yz(const t_vector<type>& v) const { return inner_product(v, yz_axis_vector()); }
			type dot_xz(const t_vector<type>& v) const { return inner_product(v, xz_axis_vector()); }

			type sq_magnit_xy() const { return (inner_product(*this, xy_axis_vector())); }
			type sq_magnit_yz() const { return (inner_product(*this, yz_axis_vector())); }
			type sq_magnit_xz() const { return (inner_product(*this, xz_axis_vector())); }

			type magnit_xy() const { return (lib_math::sqrt(sq_magnit_xy())); }
			type magnit_yz() const { return (lib_math::sqrt(sq_magnit_yz())); }
			type magnit_xz() const { return (lib_math::sqrt(sq_magnit_xz())); }
			#endif


			// returns length, normalizes at same time
			// no point adding a const-version of this
			// ("length" == "magnitude")
			type magnit_normalize_ref(const t_vector<type>& mask = xyz_axis_vector(), const type eps = t_tuple<type>::eps_scalar()) {
				type sq_mag = sq_magnit(mask);
				type    mag = sq_mag;

				if (normalizable(mask, eps * eps)) {
					(*this) *= mask;
					(*this) /= (mag = lib_math::sqrt(sq_mag));
				}

				return mag;
			}


			// note: these are EXTRINSIC rotation operators
			// they rotate <this> around the WORLD'S (fixed)
			// coordinate axes (ie. rot_y_ in xz-plane etc.)
			t_vector<type> rotate_x_ext(type angle_radians) const {
				const type ca = std::cos(angle_radians);
				const type sa = std::sin(angle_radians);

				// counter-clockwise rotation (+y becomes +z) for positive angles
				return (t_vector<type>(x(),  ca * y() - sa * z(), sa * y() + ca * z(),  0));
			}
			t_vector<type> rotate_y_ext(type angle_radians) const {
				const type ca = std::cos(angle_radians);
				const type sa = std::sin(angle_radians);

				// counter-clockwise rotation (+z becomes +x) for positive angles
				// +x is (1, 0, 0) and points RIGHT (seen from above)
				// +z is (0, 0, 1) and points DOWN (seen from above)
				return (t_vector<type>(ca * x() + sa * z(),  y(),  -sa * x() + ca * z(),  0));
			}
			t_vector<type> rotate_z_ext(type angle_radians) const {
				const type ca = std::cos(angle_radians);
				const type sa = std::sin(angle_radians);

				// counter-clockwise rotation (+x becomes +y) for positive angles
				return (t_vector<type>(ca * x() - sa * y(), sa * x() + ca * y(),  z(),  0));
			}


			// reflect this vector about the vector <N>
			t_vector<type> reflect(const t_vector<type>& n) const {
				const t_vector<type>& v = *this;
				const t_vector<type>  r = n * (n.inner_product(v) * 2);
				return (v - r);
			}

			// refract this vector about the vector <N>
			//   n1 is the refr. index of the material being exited
			//   n2 is the refr. index of the material being entered
			t_vector<type> refract(const t_vector<type>& n, type r1, type r2) const {
				const t_vector<type>& v = *this;

				const type r  = r1 / r2;
				const type ct = -(n.inner_product(v));
				const type st = 1 - (r * r * (1 - ct * ct));

				if (st > 0)
					return ((v * r) + (n * (r * ct - lib_math::sqrt(st))));

				// total internal reflection
				return n;
			}


			// NOTE: incorrect use of mask here
			t_vector<type> outer_product(const t_vector<type>& vec, const t_vector<type>& mask = xyz_axis_vector()) const {
				const type _x = ((y() * vec.z()) - (z() * vec.y())) * mask.x();
				const type _y = ((z() * vec.x()) - (x() * vec.z())) * mask.y();
				const type _z = ((x() * vec.y()) - (y() * vec.x())) * mask.z();
				return (t_vector<type>(_x, _y, _z, 0));
			}

			// note: these divide instead of multiplying by
			// the reciprocal of length (<type> can be int)
			t_vector<type> normalize(const t_vector<type>& mask = xyz_axis_vector(), const type eps = t_tuple<type>::eps_scalar()) const {
				t_vector<type> r = *this;

				if (normalizable(mask, eps * eps)) {
					return ((r * mask) / lib_math::sqrt(sq_magnit(mask)));
				}

				return r;
			}
			/*
			t_vector<type>& normalize_ref(const t_vector<type>& mask = xyz_axis_vector(), const type eps = t_tuple<type>::eps_scalar()) {
				t_vector<type>& r = *this;

				if (normalizable(mask, eps * eps)) {
					return ((r *= mask) /= lib_math::sqrt(sq_magnit(mask)));
				}

				return r;
			}
			*/

			t_vector<type>& normalize_ref() { return ((*this) = normalize()); }
			t_vector<type>& outer_product_ref(const t_vector<type>& v) { return ((*this) = outer_product(v)); }



			#if 0
			t_vector<type> project(const t_vector<type>& m) const {
				return (t_vector<type>(x() * m.x(), y() * m.y(), z() * m.z(), w() * m.w()));
			}

			// these are redundant with multiplicative-masks
			t_vector<type> xz() const { return (t_vector<type>(x(), 0,   z(), w())); }
			t_vector<type> xy() const { return (t_vector<type>(x(), y(), 0,   w())); }
			t_vector<type> yz() const { return (t_vector<type>(0,   y(), z(), w())); }
			#endif


			#if 0
			t_vector<type> project_xz() const { return (project(xz_axis_vector())); }
			t_vector<type> project_xy() const { return (project(xy_axis_vector())); }
			t_vector<type> project_yz() const { return (project(yz_axis_vector())); }

			t_vector<type>& project_xz_ref() { return ((*this) = project_xz()); }
			t_vector<type>& project_xy_ref() { return ((*this) = project_xy()); }
			t_vector<type>& project_yz_ref() { return ((*this) = project_yz()); }
			#endif

			t_vector<type> abs(                       ) const { return (t_vector<type>(std::abs(      x()       ), std::abs(      y()       ), std::abs(      z()       ), w())); }
			t_vector<type> max(const t_vector<type>& v) const { return (t_vector<type>(std::max<type>(x(), v.x()), std::max<type>(y(), v.y()), std::max<type>(z(), v.z()), w())); }
			t_vector<type> min(const t_vector<type>& v) const { return (t_vector<type>(std::min<type>(x(), v.x()), std::min<type>(y(), v.y()), std::min<type>(z(), v.z()), w())); }

			t_vector<type> signs() const { return (t_vector<type>(signum(x()), signum(y()), signum(z()), signum(w()))); }
			t_vector<type> clamp(const t_vector<type>& vmin, const t_vector<type>& vmax) const {
				t_vector<type> ret;
				ret.x() = clamp(x(), vmin.x(), vmax.x());
				ret.y() = clamp(y(), vmin.y(), vmax.y());
				ret.z() = clamp(z(), vmin.z(), vmax.z());
				ret.w() = clamp(w(), vmin.w(), vmax.w());
				return ret;
			}


			t_vector<type>      randomize(std::function<type(void)>& rng) const { return (t_vector<type>((rng() * 2) - 1, (rng() * 2) - 1, (rng() * 2) - 1)); }
			t_vector<type> unit_randomize(std::function<type(void)>& rng, const type eps = t_tuple<type>::eps_scalar()) const {
				t_vector<type> r = xyz_axis_vector();

				while (std::abs(r.sq_magnit() - 1) > eps) {
					r.x() = (rng() * 2) - 1;
					r.y() = (rng() * 2) - 1;
					r.z() = (rng() * 2) - 1;
				}

				return r;
			}

			t_vector<type>&      randomize_ref(std::function<type(void)>& rng                                              ) { return ((*this) =      randomize(   )); }
			t_vector<type>& unit_randomize_ref(std::function<type(void)>& rng, const type eps = t_tuple<type>::eps_scalar()) { return ((*this) = unit_randomize(eps)); }


			unsigned int hash(const t_vector<type>& mask = ones_vector()) const {
				return ((t_tuple<type>(*this)).hash(t_tuple<type>(mask)));
			}

			bool equals(const t_vector<type>& vec, const t_vector<type>& eps = eps_vector()) const {
				return ((t_tuple<type>(*this)).equals(t_tuple<type>(vec), t_tuple<type>(eps)));
			}
			// vector is normalizable if its (masked) magnitude != {0 or 1}
			bool normalizable(const t_vector<type>& mask, const type eps = t_tuple<type>::eps_scalar()) const {
				return (sq_magnit(mask) > eps && std::abs(1 - sq_magnit(mask)) > eps);
			}

			#if 0
			bool normalizable_xz(const type eps = t_tuple<type>::eps_scalar()) const { return (normalizable(xz_axis_vector, eps * eps)); }
			bool normalizable_xy(const type eps = t_tuple<type>::eps_scalar()) const { return (normalizable(xy_axis_vector, eps * eps)); }
			bool normalizable_yz(const type eps = t_tuple<type>::eps_scalar()) const { return (normalizable(yz_axis_vector, eps * eps)); }
			#endif

			void print() const { (t_tuple<type>(*this)).print(); }
			void sanity_assert() const { (t_tuple<type>(*this)).sanity_assert(); }

			const type* xyzw() const { return &m_values[0]; }
			      type* xyzw()       { return &m_values[0]; }

			type  operator [] (unsigned int n) const { return m_values[n]; }
			type& operator [] (unsigned int n)       { return m_values[n]; }

			type  x() const { return m_values[0]; }
			type  y() const { return m_values[1]; }
			type  z() const { return m_values[2]; }
			type  w() const { return m_values[3]; }
			type& x()       { return m_values[0]; }
			type& y()       { return m_values[1]; }
			type& z()       { return m_values[2]; }
			type& w()       { return m_values[3]; }


			static const t_vector<type>&      eps_vector();
			static const t_vector<type>&      inf_vector();
			static const t_vector<type>&     zero_vector();
			static const t_vector<type>&     ones_vector();
			static const t_vector<type>&    error_vector();
			static const t_vector<type>&   x_axis_vector();
			static const t_vector<type>&   y_axis_vector();
			static const t_vector<type>&   z_axis_vector();
			static const t_vector<type>&   w_axis_vector();
			static const t_vector<type>&  xz_axis_vector();
			static const t_vector<type>&  xy_axis_vector();
			static const t_vector<type>&  yz_axis_vector();
			static const t_vector<type>& xyz_axis_vector();

		private:
			type m_values[MATH_VECTOR_SIZE];
		};


		typedef t_vector< int32_t> t_vector4i;
		typedef t_vector<uint32_t> t_vector4ui;

		typedef t_vector< int64_t> t_vector4li;
		typedef t_vector<uint64_t> t_vector4lui;

		typedef t_vector< float> t_vector4f;
		typedef t_vector<double> t_vector4d;

		// aliases for less typing
		typedef t_vector4i  t_vec4i;
		typedef t_vector4ui t_vec4ui;

		typedef t_vector4li  t_vec4li;
		typedef t_vector4lui t_vec4lui;

		typedef t_vector4f t_vec4f;
		typedef t_vector4d t_vec4d;
	};
};

#endif

