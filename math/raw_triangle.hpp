#ifndef MATH_TRIANGLE_HDR
#define MATH_TRIANGLE_HDR

#include "./tuple.hpp"
#include "./point.hpp"
#include "./vector.hpp"
#include "./plane.hpp"

namespace newtonics {
	namespace lib_math {
		template<typename type> struct t_matrix;

		struct t_plane;
		struct t_ray;
		struct t_idx_triangle;

		struct t_raw_triangle {
		typedef float m_coor_type;

		typedef t_tuple<m_coor_type> m_tuple_type;
		typedef t_point<m_coor_type> m_point_type;
		typedef t_vector<m_coor_type> m_vector_type;
		typedef t_matrix<m_coor_type> m_matrix_type;

		public:
			t_raw_triangle();
			t_raw_triangle(const m_point_type& vert_a, const m_point_type& vert_b, const m_point_type& vert_c);
			t_raw_triangle(const m_point_type& vert_a, const m_point_type& vert_b, const m_point_type& vert_c, const m_vector_type& normal);
			t_raw_triangle(const t_raw_triangle& t) { *this = t; }

			t_raw_triangle& operator = (const t_raw_triangle& triangle);

			t_raw_triangle& invert_winding();
			t_raw_triangle transform(const m_matrix_type& matrix) const;

			t_plane to_plane() const {
				t_vec4f n;
				n.x() = m_normal.x();
				n.y() = m_normal.y();
				n.z() = m_normal.z();
				n.w() = m_normal.inner_product(m_verts[0] - t_pos4f::zero_point());
				return (t_plane(n));
			}

			m_point_type calc_midpoint() const {
				m_point_type mid_point;
				mid_point += (m_verts[0].to_vector() / m_coor_type(3));
				mid_point += (m_verts[1].to_vector() / m_coor_type(3));
				mid_point += (m_verts[2].to_vector() / m_coor_type(3));
				return mid_point;
			}

			const m_point_type& get_vertex(unsigned int idx) const { return m_verts[idx]; }
			      m_point_type& get_vertex(unsigned int idx)       { return m_verts[idx]; }
			const m_vector_type& get_normal() const { return m_normal; }
			      m_vector_type& get_normal()       { return m_normal; }

			m_coor_type  get_radius() const { return m_radius; }
			m_coor_type& get_radius()       { return m_radius; }

			m_coor_type calc_radius() const;
			m_coor_type origin_distance() const;


			bool point_in_triangle(const m_point_type& point, const m_coor_type eps = m_tuple_type::eps_scalar()) const;

			// test if <this> intersects <plane>
			bool intersect_plane(const t_plane& plane, const m_coor_type eps = m_tuple_type::eps_scalar()) const;

			// test if <this> intersects <triangle>; both
			// triangles must be in same coordinate-space
			//
			// if so, exactly one or two edges of <this> will
			// cross the plane in which <triangle> is embedded
			// and their intersection point will lie inside it
			//
			// note that is not a symmetric test!
			//
			unsigned int intersect_triangle_cop(const t_raw_triangle& triangle, m_point_type* points, const m_coor_type eps = m_tuple_type::eps_scalar()) const;
			unsigned int intersect_triangle(const t_raw_triangle& triangle, m_point_type* points, const m_coor_type eps = m_tuple_type::eps_scalar()) const;

			unsigned int intersect_triangle(
				const t_idx_triangle& idx_triangle,
				const m_point_type* idx_verts,
				      m_point_type* points,
				const m_coor_type eps = m_tuple_type::eps_scalar()
			) const;

			// maximum number of intersection points this triangle can have with another
			// for two coplanar triangles A and B up to six edge-intersection points can
			// exist, and up to three vertices of triangle A can lie inside B (or v.v.)
			//
			static constexpr unsigned int get_max_intersect_points() { return 6; }

		private:
			// the usual indexed representation is (more) complex
			// since we need access to the geometry for intersect
			// testing, etc
			m_point_type m_verts[3];
			m_vector_type m_normal;

			// either calculated or set via normal.w
			m_coor_type m_radius;
		};

		typedef t_raw_triangle t_raw_tri;
	};
};

#endif

