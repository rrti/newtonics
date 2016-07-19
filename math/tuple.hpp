#ifndef MATH_TUPLE_HDR
#define MATH_TUPLE_HDR

#include <cassert>
#include <cmath>
#include <cstdint>

#include "./math_defs.hpp"
#include "./template_funcs.hpp"

namespace newtonics {
	namespace lib_math {
		template<typename type> struct t_point;
		template<typename type> struct t_vector;

		template<typename type> struct t_tuple {
		public:
			t_tuple<type>(const type _x = 0, const type _y = 0, const type _z = 0, const type _w = 0) {
				x() = _x; y() = _y;
				z() = _z; w() = _w;
			}
			// safer than const type*
			t_tuple<type>(const type v[MATH_TUPLE_SIZE]) {
				x() = v[0]; y() = v[1];
				z() = v[2]; w() = v[3];
			}
			t_tuple<type>(const t_point<type>& p) {
				x() = p.x(); y() = p.y();
				z() = p.z(); w() = p.w();
			}
			t_tuple<type>(const t_vector<type>& v) {
				x() = v.x(); y() = v.y();
				z() = v.z(); w() = v.w();
			}
			t_tuple<type>(const t_tuple<type>& t) {
				*this = t;
			}

			bool operator == (const t_tuple<type>& t) const { return ( equals(t)); }
			bool operator != (const t_tuple<type>& t) const { return (!equals(t)); }

			t_tuple<type>& operator = (const t_tuple<type>& t) {
				x() = t.x(); y() = t.y();
				z() = t.z(); w() = t.w();
				return *this;
			}

			t_tuple<type> operator + (const t_tuple<type>& t) const { return (t_tuple<type>(x() + t.x(), y() + t.y(), z() + t.z(), w() + t.w())); }
			t_tuple<type> operator - (const t_tuple<type>& t) const { return (t_tuple<type>(x() - t.x(), y() - t.y(), z() - t.z(), w() - t.w())); }

			t_tuple<type> operator * (const type s) const { return (t_tuple<type>(x() * s, y() * s, z() * s, w() * s)); }
			t_tuple<type> operator / (const type s) const { return (t_tuple<type>(x() / s, y() / s, z() / s, w() / s)); }


			unsigned int hash(const t_tuple<type>& mask = ones_tuple()) const {
				const unsigned int hx = (*reinterpret_cast<unsigned int*>(const_cast<type*>(&m_values[0]))) * mask.x();
				const unsigned int hy = (*reinterpret_cast<unsigned int*>(const_cast<type*>(&m_values[1]))) * mask.y();
				const unsigned int hz = (*reinterpret_cast<unsigned int*>(const_cast<type*>(&m_values[2]))) * mask.z();
				const unsigned int hw = (*reinterpret_cast<unsigned int*>(const_cast<type*>(&m_values[3]))) * mask.w();
				return (hx ^ hy ^ hz ^ hw);
			}
			static unsigned int hash(const t_point<type>& pnt, const t_point<type>& mask) { return 0; }
			static unsigned int hash(const t_vector<type>& vec, const t_vector<type>& mask) { return 0; }

			bool equals(const t_tuple<type>& t, const t_tuple<type>& eps = eps_tuple()) const {
				unsigned int mask = 0;
				mask += lib_math::fp_eq<type>(x(), t.x(), eps.x());
				mask += lib_math::fp_eq<type>(y(), t.y(), eps.y());
				mask += lib_math::fp_eq<type>(z(), t.z(), eps.z());
				mask += lib_math::fp_eq<type>(w(), t.w(), eps.w());
				return (mask == 4);
			}

			bool is_point() const { return (w() == 1); }
			bool is_vector() const { return (w() == 0); }

			static bool equals(const t_point<type>& pnt_a, const t_point<type>& pnt_b, const t_point<type>& eps) { return false; }
			static bool equals(const t_vector<type>& vec_a, const t_vector<type>& vec_b, const t_vector<type>& eps) { return false; }

			void print() const;
			void sanity_assert() const {
				assert(!std::isnan(x()) && !std::isinf(x()));
				assert(!std::isnan(y()) && !std::isinf(y()));
				assert(!std::isnan(z()) && !std::isinf(z()));
				assert(!std::isnan(w()) && !std::isinf(w()));
			}

			const type* xyzw() const { return &m_values[0]; }
			      type* xyzw()       { return &m_values[0]; }

			// no need to make functors from tuples
			// operator const type* () const { return xyzw(); }
			// operator       type* ()       { return xyzw(); }

			type  operator [] (unsigned int n) const { return m_values[n]; }
			type& operator [] (unsigned int n)       { return m_values[n]; }

			// same as operator []
			// type  i(std::size_t n) const { return m_values[n]; }
			// type& i(std::size_t n)       { return m_values[n]; }
			type  x() const { return m_values[0]; }
			type  y() const { return m_values[1]; }
			type  z() const { return m_values[2]; }
			type  w() const { return m_values[3]; }
			type& x()       { return m_values[0]; }
			type& y()       { return m_values[1]; }
			type& z()       { return m_values[2]; }
			type& w()       { return m_values[3]; }

			static const t_tuple<type>& zero_tuple();
			static const t_tuple<type>& ones_tuple();
			static const t_tuple<type>&  eps_tuple();

			static type eps_scalar() { return ((eps_tuple()).x()); }

		private:
			type m_values[MATH_TUPLE_SIZE];
		};


		typedef t_tuple< int32_t> t_tuple4i;
		typedef t_tuple<uint32_t> t_tuple4ui;

		typedef t_tuple< int64_t> t_tuple4li;
		typedef t_tuple<uint64_t> t_tuple4lui;

		typedef t_tuple< float> t_tuple4f;
		typedef t_tuple<double> t_tuple4d;

		// aliases for less typing
		typedef t_tuple4i  t_tup4i;
		typedef t_tuple4ui t_tup4ui;

		typedef t_tuple4li  t_tup4li;
		typedef t_tuple4lui t_tup4lui;

		typedef t_tuple4f t_tup4f;
		typedef t_tuple4d t_tup4d;
	};
};

#endif

