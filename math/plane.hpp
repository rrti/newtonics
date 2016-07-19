#ifndef MATH_PLANE_HDR
#define MATH_PLANE_HDR

#include "./tuple.hpp"
#include "./point.hpp"
#include "./vector.hpp"

namespace newtonics {
	namespace lib_math {
		struct t_raw_triangle;
		struct t_idx_triangle;

		struct t_plane {
		typedef float m_coor_type;

		typedef t_tuple<m_coor_type> m_tuple_type;
		typedef t_point<m_coor_type> m_point_type;
		typedef t_vector<m_coor_type> m_vector_type;

		public:
			t_plane(const m_vector_type& n = m_vector_type::y_axis_vector()) {
				m_normal = n;
			}

			const m_point_type get_point() const { return (m_point_type::zero_point() + m_normal * m_vector_type::xyz_axis_vector() * m_normal.w()); }
			const m_vector_type& get_normal() const { return m_normal; }

			// note: assumes |m_normal.xyz| = 1
			m_coor_type origin_distance() const { return m_normal.w(); }
			m_coor_type point_distance(const m_point_type& p) const {
				return ((p.to_vector()).inner_product(m_normal) - m_normal.w());
			}

			// test the interval [pos, pos + vel*t] with t in [0, 1]
			// t < 0 implies no collision at all or not this time-step
			m_coor_type get_collision_time(
				const m_point_type& pos,
				const m_vector_type& vel,
				const m_coor_type tmin = m_coor_type(0),
				const m_coor_type tmax = m_coor_type(1),
				const m_coor_type eps = m_tuple_type::eps_scalar()
			) const;


			void clip_triangle_aux(
				const t_raw_triangle& triangle,
				const t_vector4i& indices,
				      t_raw_triangle* clipped_tris
			) const;


			// clip <triangle> against arbitrarily oriented (infinite) plane <this>
			unsigned int clip_triangle(
				const t_raw_triangle& triangle,
				      t_raw_triangle* clipped_tris,
				const m_coor_type eps = m_tuple_type::eps_scalar()
			) const;
			unsigned int clip_triangle(
				const t_idx_triangle& idx_triangle,
				const m_point_type* idx_vertices,
				      t_raw_triangle* clipped_tris,
				const m_coor_type eps = m_tuple_type::eps_scalar()
			) const;

		private:
			// represents a plane ax + by + cz + d = 0 (when N = <a,b,c>
			// is a unit-vector, d is the raw orthogonal distance along
			// N from the plane to origin and plane is in Hessian normal
			// form)
			//
			//   nx = a / SQRT(a*a + b*b + c*c)
			//   ny = b / SQRT(a*a + b*b + c*c)
			//   nz = c / SQRT(a*a + b*b + c*c)
			//   nd = d / SQRT(a*a + b*b + c*c)
			//
			// unit normal vector; w-component stores distance to origin!
			// (direction defines which side of the plane is its "front")
			m_vector_type m_normal;
		};
	};
};

#endif

