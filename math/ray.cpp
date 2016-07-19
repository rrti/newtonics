#include <cmath>
#include <limits>

#include "./ray.hpp"
#include "./matrix.hpp"
#include "./plane.hpp"
#include "./raw_triangle.hpp"
#include "./idx_triangle.hpp"
#include "./template_funcs.hpp"

namespace newtonics {
	namespace lib_math {
		t_ray::t_ray() {
			m_dir.z() = m_coor_type(1);
			m_len     = m_coor_type(1);
		}
		t_ray::t_ray(const m_point_type& p, const m_vector_type& d) {
			m_pos = p;
			m_dir = d;
			m_len = std::numeric_limits<m_coor_type>::max();
			// m_len = std::numeric_limits<m_coor_type>::infinity();
		}
		t_ray::t_ray(const m_point_type& p0, const m_point_type& p1) {
			assert(p1 != p0);

			m_pos  = p0;
			m_dir  = p1 - p0;
			m_len  = m_dir.len();
			m_dir /= m_len;
		}

		t_ray& t_ray::operator = (const t_ray& r) {
			m_pos = r.m_pos;
			m_dir = r.m_dir;
			m_len = r.m_len;
			return *this;
		}

		// flip direction, keep start-position
		t_ray t_ray::invert_dir() const {
			t_ray r;
			r.m_pos =  m_pos;
			r.m_dir = -m_dir;
			r.m_len =  m_len;
			return r;
		}
		// flip direction, flip start-position
		// undefined if ray is infinitely long
		t_ray t_ray::invert_pos() const {
			t_ray r;
			r.m_pos =  m_pos + m_dir * m_len;
			r.m_dir = -m_dir;
			r.m_len =  m_len;
			return r;
		}
		t_ray t_ray::transform(const m_matrix_type& matrix) const {
			assert(m_pos.w() == m_coor_type(1));

			t_ray r;
			r.m_pos = matrix * (m_pos        );
			r.m_dir = matrix * (m_pos + m_dir) - r.m_pos;
			r.m_len = m_len;
			return r;
		}



		t_ray::m_coor_type t_ray::plane_intersect_distance(const t_plane& plane) const {
			return (plane_intersect_distance(plane.get_normal(), plane.origin_distance()));
		}

		t_ray::m_coor_type t_ray::plane_intersect_distance(
			const m_vector_type& plane_normal,
			const m_coor_type plane_origin_distance,
			const t_ray::m_coor_type eps
		) const {
			// get distance from ray-origin to world-origin along normal
			// then offset this distance by how far plane itself is from
			// world-origin
			const m_coor_type p = plane_normal.inner_product(m_pos.to_vector());
			const m_coor_type d = plane_normal.inner_product(m_dir);
			const m_coor_type n = plane_origin_distance - p;

			m_coor_type r = m_coor_type(-1);

			if (std::abs(d) > eps) {
				r = ((lib_math::sign(p) != m_coor_type(1))? -std::abs(n / d): (n / d));
			} else {
				// we are parallel to plane and can never intersect
				// sign depends on whether pos is above or below it
				//
				// don't use infinity, it allows NaN's to percolate
				// r = std::numeric_limits<m_coor_type>::infinity() * s;
				r = std::numeric_limits<m_coor_type>::max() * lib_math::sign(p);
			}

			return r;
		}


		/*
		m_coor_type t_ray::plane_distance(const t_plane& plane) {
			const m_coor_type dp = m_dir.inner_product(plane.normal());

			// parallel to plane
			if (dp == 0)
				return type(-1);

			// "n / d"
			return ((plane.origin_distance() - (m_pos - zero_point).inner_product(plane.normal())) / dp);
		}
		*/
		/*
		t_ray::m_coor_type t_ray::point_distance(const m_point_type& point) const {
			const m_vector_type  v = point - m_pos;
			const m_vector_type& d = m_dir;

			// orthogonally project <v> onto <d> (which is normalized)
			const m_point_type p = m_pos + (d * v.inner_product(d));
			const m_vector_type o = point - p;

			// NOTE:
			//   the returned distance is always non-negative, so we
			//   do not know if <point> lies above or below the line
			//
			//   the projected point might not fall on the segment
			//   (for half-infinite segments this happens whenever
			//   v.dot(d) is negative), what to return then?
			//
			return (o.len());
		}
		*/
		/*
		t_ray::m_coor_type t_ray::ray_distance(const t_ray& ray) const {
			const m_vector_type c = ray.pos() - m_pos;
			const m_vector_type v = m_dir.outer_product(ray.dir());

			// let
			//   line A = x1 + a*s
			//   line B = x3 + b*t
			// where
			//   a = x2 - x1
			//   b = x4 - x3
			//   c = x3 - x1
			//   v = outer_product(a, b)
			// then
			//   D = |inner_product(c, v)| / |v|
			//
			// NOTE:
			//   this only returns the minimum distance, not the points
			//   on the lines *at which* it is minimized and ONLY works
			//   for NON-intersecting NON-parallel lines (in 3D)
			return (std::abs(c.inner_product(v)) / v.len());
		}
		*/


		// point on ray at which we intersect <plane>
		t_ray::m_point_type t_ray::plane_intersection_point(const t_plane& plane) const {
			const m_coor_type  d = plane_intersect_distance(plane);
			const m_point_type p = m_pos + m_dir * d;
			const bool         b = (d >= m_coor_type(0) && d <= m_len);

			// w-coordinate of 0 means points is NOT valid (caller must check)
			return (m_point_type(p.x(), p.y(), p.z(), p.w() * b));
		}


		// distance along ray to point in plane of <triangle>
		t_ray::m_coor_type t_ray::triangle_intersect_distance(const t_raw_triangle& triangle) const {
			return (plane_intersect_distance(triangle.get_normal(), triangle.origin_distance()));
		}

		t_ray::m_coor_type t_ray::triangle_intersect_distance(const t_idx_triangle& idx_triangle, const m_point_type* idx_verts) const {
			return (triangle_intersect_distance(idx_triangle.to_raw_triangle(idx_verts)));
		}


		// point on ray at which we intersect <triangle>
		//
		// intersection point is always within the PLANE
		// of <triangle> but not! necessarily within its
		// edges
		t_ray::m_point_type t_ray::triangle_intersection_point(const t_raw_triangle& triangle) const {
			const m_coor_type  d = triangle_intersect_distance(triangle);
			const m_point_type p = m_pos + m_dir * d;
			const bool         b = (d >= m_coor_type(0) && d <= m_len && triangle.point_in_triangle(p));

			// w-coordinate of 0 means points is NOT valid (caller must check)
			return (m_point_type(p.x(), p.y(), p.z(), p.w() * b));
		}

		t_ray::m_point_type t_ray::triangle_intersection_point(const t_idx_triangle& idx_triangle, const m_point_type* idx_verts) const {
			return (triangle_intersection_point(idx_triangle.to_raw_triangle(idx_verts)));
		}


		bool t_ray::ray_intersect(const t_ray& ray, t_ray::m_point_type* point, const t_ray::m_coor_type eps) const {
			const m_point_type&  p0 =     pos();
			const m_point_type&  p1 = ray.pos();
			const m_vector_type& d0 =     dir();
			const m_vector_type& d1 = ray.dir();

			bool ret = false;

			// find t0 (on r0) and t1 (on r1) by equating r0 to r1
			// (we have three equ's and two vars, so can be solved)
			//   (0) p0   + d0   * t0  ==  p1   + d1   * t1
			//   (1) p0.x + d0.x * t0  ==  p1.x + d1.x * t1
			//   (2) p0.y + d0.y * t0  ==  p1.y + d1.y * t1
			//   (3) p0.z + d0.z * t0  ==  p1.z + d1.z * t1
			//
			// move p1.x to LHS in (1), then divide by d1.x
			//   (1) (p0.x - p1.x)      +  d0.x       * t0  ==  d1.x * t1
			//   (1) (p0.x - p1.x)/d1.x + (d0.x/d1.x) * t0  ==         t1
			//
			// substitute expression for t1 into either (2) or (3)
			//   (2) p0.y        + d0.y * t0  ==  p1.y        + d1.y *  ((p0.x - p1.x)/d1.x  +        (d0.x/d1.x) * t0)
			//   (2) p0.y        + d0.y * t0  ==  p1.y        + d1.y *  ((p0.x - p1.x)/d1.x) + d1.y * (d0.x/d1.x) * t0
			//   (2) p0.y        + d0.y * t0  ==  p1.y        + B                            + C                  * t0
			//   (2)               d0.y * t0  ==  p1.y - p0.y + B                            + C                  * t0
			//   (2)               d0.y * t0  ==  A           + B                            + C                  * t0
			//   (2)  -C * t0    + d0.y * t0  ==  A           + B
			//   (2) (-C + d0.y)        * t0  ==  A           + B
			//   (2) D                  * t0  ==  A           + B
			//   (2)                      t0  == (A           + B) / D
			//
			for (unsigned int i = 0; (i < 3 && !ret); i++) {
				if (std::fabs(d1[i]) < eps)
					continue;

				for (unsigned int j = 0; (j < 2 && !ret); j++) {
					// try each of the two substitution targets
					const unsigned n = (i + j) % 3;

					const float An = p1[n] - p0[n];
					const float Bn = d1[n] * ((p0[i] - p1[i]) / d1[i]);
					const float Cn = d1[n] * (d0[i] / d1[i]);
					const float Dn = -Cn + d0[n];

					if (std::fabs(Dn) < eps)
						continue;

					const float t0 = (An + Bn) / Dn;
					const float t1 = (p0[i] - p1[i]) / d1[i] + (d0[i] / d1[i]) * t0;

					// check if t0 and t1 are valid (rays represent finite line-segments)
					ret |= (t0 >= 0.0f && t0 <=     len());
					ret |= (t1 >= 0.0f && t1 <= ray.len());

					if (ret) {
						const m_point_type x = p0 + d0 * t0;
						const m_point_type y = p1 + d1 * t1;

						*point = ((x.to_vector() + y.to_vector()) * 0.5f).to_point();
						break;
					}
				}
			}

			return ret;
		}
	};
};

