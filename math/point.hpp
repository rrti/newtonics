#ifndef MATH_POINT_HDR
#define MATH_POINT_HDR

#include <cmath>
#include <cstdint>

#include "./math_defs.hpp"

namespace newtonics {
	namespace lib_math {
		// we do *not* want these to inherit from t_tuple
		// for type-safety reasons: (point - point) should
		// always be a vector, etc.
		template<typename type> struct t_tuple;
		template<typename type> struct t_vector;

		template<typename type> struct t_point {
		public:
			t_point<type>(const type _x = 0, const type _y = 0, const type _z = 0, const type _w = 1) {
				x() = _x; y() = _y;
				z() = _z; w() = _w;
			}
			t_point<type>(const type v[MATH_POINT_SIZE]) {
				x() = v[0]; y() = v[1];
				z() = v[2]; w() = v[3];
			}
			t_point<type>(const t_point<type>& p) {
				*this = p;
			}

			bool operator == (const t_point<type>& p) const { return ( equals(p)); }
			bool operator != (const t_point<type>& p) const { return (!equals(p)); }

			// rules:
			//   type(point + vector) := type(point)
			//   type(point - vector) := type(point + -vector)
			//   type(point - point) := type(vector)
			//   type(point + point) := UNDEFINED
			//   points: w=1, vectors: w=0
			t_point<type> operator + (const t_vector<type>& v) const { return (t_point<type>(x() + v.x(), y() + v.y(), z() + v.z(), w() + v.w())); }
			t_point<type> operator - (const t_vector<type>& v) const { return (t_point<type>(x() - v.x(), y() - v.y(), z() - v.z(), w() - v.w())); }
			t_vector<type> operator - (const t_point<type>& p) const { return (t_vector<type>(x() - p.x(), y() - p.y(), z() - p.z(), w() - p.w())); }

			t_point<type>& operator  = (const t_point<type>& p)  { x()  = p.x(); y()  = p.y(); z()  = p.z(); w()  = p.w(); return *this; }
			t_point<type>& operator += (const t_vector<type>& v) { x() += v.x(); y() += v.y(); z() += v.z(); w() += v.w(); return *this; }
			t_point<type>& operator -= (const t_vector<type>& v) { x() -= v.x(); y() -= v.y(); z() -= v.z(); w() -= v.w(); return *this; }

			// t_point<type> operator - () const { return (t_point<type>(-x(), -y(), -z(), -w())); }

			// cf. vector::to_point
			t_tuple<type> to_tuple() const { return (t_tuple<type>(*this)); }
			t_vector<type> to_vector() const { assert(w() == 1); return ((*this) - zero_point()); }


			unsigned int hash(const t_point<type>& mask = ones_point()) const {
				return ((t_tuple<type>(*this)).hash(t_tuple<type>(mask)));
			}

			bool equals(const t_point<type>& pnt, const t_point<type>& eps = eps_point()) const {
				return ((t_tuple<type>(*this)).equals(t_tuple<type>(pnt), t_tuple<type>(eps)));
			}

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


			static const t_point<type>&  eps_point();
			static const t_point<type>& null_point();
			static const t_point<type>& zero_point();
			static const t_point<type>& ones_point();

		private:
			type m_values[MATH_POINT_SIZE];
		};


		typedef t_point< int32_t> t_point4i;
		typedef t_point<uint32_t> t_point4ui;

		typedef t_point< int64_t> t_point4li;
		typedef t_point<uint64_t> t_point4lui;

		typedef t_point< float> t_point4f;
		typedef t_point<double> t_point4d;

		// aliases for less typing
		typedef t_point4i  t_pos4i;
		typedef t_point4ui t_pos4ui;

		typedef t_point4li  t_pos4li;
		typedef t_point4lui t_pos4lui;

		typedef t_point4f t_pos4f;
		typedef t_point4d t_pos4d;
	};
};

#endif

