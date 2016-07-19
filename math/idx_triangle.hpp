#ifndef MATH_IDX_TRIANGLE_HDR
#define MATH_IDX_TRIANGLE_HDR

#include <algorithm>
#include <cassert>
#include <vector>

#include "./tuple.hpp"
#include "./point.hpp"
#include "./vector.hpp"
#include "./raw_triangle.hpp"

namespace newtonics {
	namespace lib_math {
		struct t_idx_triangle {
		typedef float m_coor_type;
		typedef t_point<m_coor_type> m_point_type;
		typedef t_vector<m_coor_type> m_vector_type;
		typedef t_tuple<unsigned int> m_tuple_type;
		public:
			t_idx_triangle(const std::vector<m_point_type>& vertices, const m_tuple_type& indices) {
				assert(indices.x() < vertices.size());
				assert(indices.y() < vertices.size());
				assert(indices.z() < vertices.size());
				assert(vertices.size() >= 3);

				*this = t_idx_triangle(&vertices[0], indices);
			}
			t_idx_triangle(const m_point_type* vertices, const m_tuple_type& indices) {
				const m_vector_type v10 = (vertices[indices.y()] - vertices[indices.x()]).normalize();
				const m_vector_type v20 = (vertices[indices.z()] - vertices[indices.x()]).normalize();

				// BIG FAT NOTE:
				//   the cross-product of two normalized vectors is *NOT*
				//   always normalized itself, so we *MUST* re-normalize!
				//   e.g.: <1, 0, 0> x <0.71, 0, 0.71> = <0, 0.71, 0>
				//
				m_indices = indices;
				m_normal = (v20.outer_product((v10))).normalize();
				m_radius = calc_radius(vertices, indices);
			}
			t_idx_triangle(const m_point_type* vertices, const m_tuple_type& indices, const m_vector_type& normal) {
				m_indices = indices;
				m_normal = normal;
				m_radius = calc_radius(vertices, indices);
			}

			t_raw_triangle to_raw_triangle(const m_point_type* vertices) const {
				const m_point_type& va = vertices[m_indices.x()];
				const m_point_type& vb = vertices[m_indices.y()];
				const m_point_type& vc = vertices[m_indices.z()];
				const m_vector_type  n = m_normal + m_vector_type::w_axis_vector() * m_radius;
				return (t_raw_triangle(va, vb, vc, n));
			}


			const m_tuple_type& get_indices() const { return m_indices; }
			const m_vector_type& get_normal() const { return m_normal; }

			m_coor_type get_radius() const { return m_radius; }
			m_coor_type calc_radius(const m_point_type* vertices, const m_tuple_type& indices) const {
				m_point_type mid_point = calc_midpoint(vertices, indices);
				m_vector_type mid_vector = vertices[indices.x()] - mid_point;

				if ((vertices[indices.y()] - mid_point).sq_magnit() > mid_vector.sq_magnit())
					mid_vector = vertices[indices.y()] - mid_point;
				if ((vertices[indices.z()] - mid_point).sq_magnit() > mid_vector.sq_magnit())
					mid_vector = vertices[indices.z()] - mid_point;

				return (mid_vector.magnit());
			}

			m_point_type calc_midpoint(const m_point_type* vertices, const m_tuple_type& indices) const {
				m_point_type mid_point;
				mid_point += (vertices[indices.x()].to_vector() / m_coor_type(3));
				mid_point += (vertices[indices.y()].to_vector() / m_coor_type(3));
				mid_point += (vertices[indices.z()].to_vector() / m_coor_type(3));
				return mid_point;
			}

		private:
			m_tuple_type m_indices;
			m_vector_type m_normal;

			m_coor_type m_radius;
		};

		typedef t_idx_triangle t_idx_tri;
	};
};

#endif

