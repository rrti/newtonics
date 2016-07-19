#include <algorithm>
#include <limits>

#include "./plane.hpp"
#include "./ray.hpp"
#include "./raw_triangle.hpp"
#include "./idx_triangle.hpp"
#include "./template_funcs.hpp"

namespace newtonics {
	namespace lib_math {
		// test the interval [pos, pos + vel*t] with t in [0, 1]
		// t < 0 implies no collision at all or not this time-step
		t_plane::m_coor_type t_plane::get_collision_time(
			const m_point_type& pos,
			const m_vector_type& vel,
			const m_coor_type tmin,
			const m_coor_type tmax,
			const m_coor_type eps
		) const {
			const m_coor_type d0 = point_distance(pos + vel * tmin);
			const m_coor_type d1 = point_distance(pos + vel * tmax);

			// segment interval starts behind plane, no collision possible
			if (d0 < m_coor_type(0))
				return (-std::numeric_limits<m_coor_type>::max());
			// segment interval ends in front of plane, no collision possible
			if (d1 > m_coor_type(0))
				return (-std::numeric_limits<m_coor_type>::max());
			// segment interval starts in front of and ends behind plane
			// and has been narrowed to a point, stop further subdivision
			if (std::abs(tmax - tmin) < eps)
				return ((tmin + tmax) * m_coor_type(0.5f));

			// binary-search approach, recurse (note that for
			// point-vs-plane we could solve for t analytically)
			const m_coor_type t0 = get_collision_time(pos, vel, tmin, (tmin + tmax) * m_coor_type(0.5f));
			const m_coor_type t1 = get_collision_time(pos, vel, (tmin + tmax) * m_coor_type(0.5f), tmax);

			return (std::max<m_coor_type>(t0, t1));
		}



		void t_plane::clip_triangle_aux(
			const t_raw_triangle& triangle,
			const t_vector4i& indices,
			      t_raw_triangle* clipped_tris
		) const {
			const m_vector_type& norm  = triangle.get_normal();
			const m_point_type& vert_a = triangle.get_vertex(indices.x()); // 1st vertex "above" <this>
			const m_point_type& vert_b = triangle.get_vertex(indices.y()); // 2nd vertex "above" <this>
			const m_point_type& vert_c = triangle.get_vertex(indices.z()); //     vertex "below" <this>

			const t_ray edge_ca = t_ray(vert_c, vert_a), edge_ac = edge_ca.invert_pos();
			const t_ray edge_cb = t_ray(vert_c, vert_b), edge_bc = edge_cb.invert_pos();

			// planes are one-sided!
			const t_raw_triangle::m_point_type vert_ca = ((edge_ca.dir()).inner_product(norm) >= t_raw_triangle::m_coor_type(0))?
				edge_ac.plane_intersection_point(*this):
				edge_ca.plane_intersection_point(*this);
			const t_raw_triangle::m_point_type vert_cb = ((edge_cb.dir()).inner_product(norm) >= t_raw_triangle::m_coor_type(0))?
				edge_bc.plane_intersection_point(*this):
				edge_cb.plane_intersection_point(*this);

			// TODO: directly create triangles with same ordering as <triangle>
			clipped_tris[0] = t_raw_triangle(vert_c,  vert_ca, vert_cb/*, norm*/);
			clipped_tris[1] = t_raw_triangle(vert_a,  vert_cb, vert_ca/*, norm*/);
			clipped_tris[2] = t_raw_triangle(vert_b,  vert_ca, vert_cb/*, norm*/);
		}

		unsigned int t_plane::clip_triangle(
			const t_raw_triangle& triangle,
			      t_raw_triangle* clipped_tris,
			const m_coor_type eps
		) const {
			const m_coor_type vertex_dists[3] = {
				point_distance(triangle.get_vertex(0)),
				point_distance(triangle.get_vertex(1)),
				point_distance(triangle.get_vertex(2)),
			};
			const m_coor_type vertex_signs[3] = {
				lib_math::sign(vertex_dists[0]),
				lib_math::sign(vertex_dists[1]),
				lib_math::sign(vertex_dists[2]),
			};

			const unsigned int num_clipped_verts =
				(std::abs(vertex_dists[0]) <= eps) +
				(std::abs(vertex_dists[1]) <= eps) +
				(std::abs(vertex_dists[2]) <= eps);

			unsigned int num_clipped_tris = 0;

			// plane does not cut through any edge or vertex
			// (ie. all vertices located on one side of plane)
			// so no need to clip
			//
			// when we _do_ need to clip, there are _always_
			// two vertices with equal sign which determine
			// the edges to clip along
			if (std::abs(vertex_signs[0] + vertex_signs[1] + vertex_signs[2]) == 3)
				return num_clipped_tris;

			switch (num_clipped_verts) {
				case 0: {
					// case #1; plane is more than epsilon-distance
					// from each vertex so it cuts strictly through
					// triangle edges to create a pyramid (t0) and
					// a trapezoid; latter is a quadrilateral which
					// can be split into two triangles t1 and t2
					//
					if (vertex_signs[0] == vertex_signs[1]) {
						// A and B are on opposite side of plane as C
						// and will become part of trapezoid quadrilat
						//
						clip_triangle_aux(triangle, t_vector4i(0, 1, 2, -1), clipped_tris);
					}

					if (vertex_signs[0] == vertex_signs[2]) {
						// A and C are on opposite side of plane as B
						// and will become part of trapezoid quadrilat
						//
						clip_triangle_aux(triangle, t_vector4i(0, 2, 1, -1), clipped_tris);
					}

					if (vertex_signs[1] == vertex_signs[2]) {
						// B and C are on opposite side of plane as A
						// and will become part of trapezoid quadrilat
						//
						clip_triangle_aux(triangle, t_vector4i(1, 2, 0, -1), clipped_tris);
					}

					num_clipped_tris = 3;
				} break;


				case 1: {
					// case #2; plane cuts through ONE (OR TWO) edges and ONE
					// vertex of triangle --> creates only two new triangles
					// (third would be degenerate)
					//
					// look up the epsilon-vertex ('a')
					for (unsigned int n = 0; n < 3; n++) {
						if (vertex_dists[n] < -eps || vertex_dists[n] > eps)
							continue;

						// now need to determine if we are dealing with
						// a trapezoidal degeneracy or legitimate split
						//
						// from <n> we can determine opposite edge, and
						// if that edge intersects <plane> then split is
						// "legitimate" (ie. if vertices n+1 and n+2 are
						// on opposite sides of <plane>)
						//
						if (vertex_signs[(n + 1) % 3] == vertex_signs[(n + 2) % 3]) {
							clipped_tris[0] = triangle;
							num_clipped_tris = 1;
						} else {
							// find intersection point with opposite edge ('bc')
							const t_ray opp_edge_bc = t_ray(triangle.get_vertex((n + 1) % 3), triangle.get_vertex((n + 2) % 3));
							const t_ray opp_edge_cb = opp_edge_bc.invert_pos();

							const t_raw_triangle::m_point_type opp_vert = ((opp_edge_bc.dir()).inner_product(m_normal) >= m_coor_type(0))?
								opp_edge_cb.plane_intersection_point(*this):
								opp_edge_bc.plane_intersection_point(*this);

							assert(opp_vert.w() == 1.0f);

							clipped_tris[0] = t_raw_triangle(triangle.get_vertex(n), triangle.get_vertex((n + 1) % 3), opp_vert/*, triangle.normal()*/);
							clipped_tris[1] = t_raw_triangle(triangle.get_vertex(n), opp_vert, triangle.get_vertex((n + 2) % 3)/*, triangle.normal()*/);
							num_clipped_tris = 2;
						}

						break;
					}
				} break;


				case 2: {
					// case #3; plane cuts strictly through TWO vertices of
					// the triangle --> no new triangles need to be created
					// so just make a copy
					clipped_tris[0] = triangle;
					num_clipped_tris = 1;
				} break;


				case 3: {
					// case #4; plane is nearly entirely parallel to triangle
					// do nothing
				} break;
			}

			return num_clipped_tris;
		}

		unsigned int t_plane::clip_triangle(
			const t_idx_triangle& idx_triangle,
			const m_point_type* idx_vertices,
			      t_raw_triangle* clipped_tris,
			const m_coor_type eps
		) const {
			return (clip_triangle(idx_triangle.to_raw_triangle(idx_vertices), clipped_tris, eps));
		}
	};
};

