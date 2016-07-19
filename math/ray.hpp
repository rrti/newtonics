#ifndef MATH_RAY_HDR
#define MATH_RAY_HDR

#include "./tuple.hpp"
#include "./point.hpp"
#include "./vector.hpp"

namespace newtonics {
	namespace lib_math {
		template<typename type> struct t_matrix;

		struct t_plane;
		struct t_raw_triangle;
		struct t_idx_triangle;

		// represents a finite-length line segment
		struct t_ray {
		typedef float m_coor_type;

		typedef t_tuple<m_coor_type> m_tuple_type;
		typedef t_point<m_coor_type> m_point_type;
		typedef t_vector<m_coor_type> m_vector_type;
		typedef t_matrix<m_coor_type> m_matrix_type;

		public:
			t_ray();
			t_ray(const m_point_type& p, const m_vector_type& d);
			t_ray(const m_point_type& p0, const m_point_type& p1);
			t_ray(const t_ray& r) { *this = r; }

			t_ray& operator = (const t_ray& r);

			// flip direction, keep start-position
			t_ray invert_dir() const;
			// flip direction, flip start-position
			// undefined if ray is infinitely long
			t_ray invert_pos() const;
			t_ray transform(const m_matrix_type& matrix) const;

			const m_point_type& pos() const { return m_pos; }
			const m_vector_type& dir() const { return m_dir; }

			m_coor_type len() const { return m_len; }

			// distance along ray to point in <plane>
			//
			// distance is negative whenever ray can not logically intersect
			// (e.g. if it starts above plane and points away from normal or
			// if it starts anywhere below plane regardless of direction)
			m_coor_type plane_intersect_distance(const t_plane& plane) const;
			m_coor_type plane_intersect_distance(
				const m_vector_type& plane_normal,
				const m_coor_type plane_origin_distance,
				const m_coor_type eps = m_tuple_type::eps_scalar()
			) const;


			m_coor_type plane_distance(const t_plane& plane);
			m_coor_type point_distance(const m_point_type& point) const;
			m_coor_type ray_distance(const t_ray& ray) const;


			// point on ray at which we intersect <plane>
			m_point_type plane_intersection_point(const t_plane& plane) const;

			// distance along ray to point in plane of <triangle>
			m_coor_type triangle_intersect_distance(const t_raw_triangle& triangle) const;
			m_coor_type triangle_intersect_distance(const t_idx_triangle& idx_triangle, const m_point_type* idx_verts) const;

			// point on ray at which we intersect <triangle>
			m_point_type triangle_intersection_point(const t_raw_triangle& triangle) const;
			m_point_type triangle_intersection_point(const t_idx_triangle& idx_triangle, const m_point_type* idx_verts) const;

			// point where <this> intersects <ray>
			bool ray_intersect(const t_ray& ray, m_point_type* point, const t_ray::m_coor_type eps) const;

		private:
			m_point_type m_pos; // origin
			m_vector_type m_dir; // direction

			// must always be >= 0 so we represent a finite line segment
			// (but setting it to +infinity makes us a semi-infinite ray)
			// note: m_dir is ALWAYS assumed to be normalized!
			m_coor_type m_len;
		};
	};
};

#endif

