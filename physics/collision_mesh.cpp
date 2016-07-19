#include <set>

#include "../global_defs.hpp"
#include "collision_mesh.hpp"
#include "../math/matrix.hpp"
#include "../math/ray.hpp"
#include "../math/template_funcs.hpp"
#include "../misc/scoped_timer.hpp"

namespace newtonics {
	namespace lib_physics {
		float t_collision_mesh::calc_bounding_box() {
			m_bounding_box.clear();

			for (size_t n = 0; n < m_vertices.size(); n++) {
				m_bounding_box.add_point(m_vertices[n]);
			}

			// calculate and assign the BB's radius
			return (m_bounding_box.finalize());
		}

		void t_collision_mesh::set_vertex_normals() {
			// NOTE:
			//   this is unreliable because vertices are shared unevenly
			//   for coldet we do not need them, only face-normals matter
			assert(!m_vnormals.empty());

			for (size_t n = 0; n < m_polygons.size(); n++) {
				const t_tup4ui& indices = m_polygons[n].get_indices();
				const t_vec4f& normal = m_polygons[n].get_normal();

				// vertices xyz all share triangle n's normal
				m_vnormals[indices.x()] += normal;
				m_vnormals[indices.y()] += normal;
				m_vnormals[indices.z()] += normal;
			}

			for (size_t n = 0; n < m_vertices.size(); n++) {
				m_vnormals[n] = m_vnormals[n].normalize();
			}
		}

		bool t_collision_mesh::point_inside_mesh(const t_pos4f& point_mspace) const {
			bool ret = true;

			for (size_t n = 0; n < m_polygons.size(); n++) {
				const t_tup4ui& indices = m_polygons[n].get_indices();
				const t_vec4f& normal = m_polygons[n].get_normal();
				const t_vec4f& vector = point_mspace - m_vertices[indices.y()];

				// LEQ, triangle surfaces count as inclusive
				ret &= (normal.inner_product(vector) <= 0.0f);
			}

			return ret;
		}


		t_pos4f t_collision_mesh::calc_ray_intersection(const t_ray& wspace_ray, const t_mat44f& wspace_to_mspace_mat) const {
			// transform ray to model-space
			const t_ray mspace_ray = wspace_ray.transform(wspace_to_mspace_mat);

			// use null-point s.t. w=0 if no valid intersection is found
			t_pos4f int_point = t_pos4f::null_point();
			t_vec4f int_vector = m_bounding_box.get_mpos() - mspace_ray.pos();

			// early-out test; collision meshes are usually simple so we do
			// not need fancy acceleration structures (for meshes which are
			// more geometrically complex than cubes we should however also
			// use bounding-box intersection)
			if (((mspace_ray.dir() * int_vector.inner_product(mspace_ray.dir())) - int_vector).sq_len() > lib_math::square(m_bounding_box.get_radius()))
				return int_point;

			// magnitude of this will correspond to closest IP
			int_vector = t_vec4f::inf_vector();

			// find closest intersection point
			for (size_t n = 0; n < m_polygons.size(); n++) {
				const t_raw_tri& raw_tri_poly   = m_polygons[n].to_raw_triangle(&m_vertices[0]);
				const t_pos4f&   tri_int_point  = mspace_ray.triangle_intersection_point(raw_tri_poly);
				const t_vec4f    tri_int_vector = tri_int_point - mspace_ray.pos();

				if (tri_int_point.w() == 0.0f)
					continue;

				if (tri_int_vector.sq_len() >= int_vector.sq_len())
					continue;

				int_point = tri_int_point;
				int_vector = tri_int_vector;
			}

			return (wspace_to_mspace_mat.invert_affine() * int_point);
		}

		t_pos4f t_collision_mesh::calc_center_of_mass(float mass) const {
			assert(!m_vertices.empty());
			assert(m_vertices.size() == m_vmasses.size());

			// defaults to (0, 0, 0)
			t_pos4f center_of_mass;

			if (mass > 0.0f) {
				for (size_t n = 0; n < m_vertices.size(); n++) {
					center_of_mass += (m_vertices[n].to_vector() * m_vmasses[n]);
				}

				center_of_mass = (center_of_mass.to_vector() / mass).to_point();
			}

			assert(point_inside_mesh(center_of_mass));
			return center_of_mass;
		}


		t_mat44f t_collision_mesh::calc_inertia_tensor() const {
			const t_mat44f I; // identity

			t_mat44f M; // 1st moment
			t_mat44f C; // 2nd moment of the mass distribution (mass-weighted covariance)

			for (size_t n = 0; n < m_vertices.size(); n++) {
				const t_vec4f  v   = m_vertices[n].to_vector();
				const t_mat44f vTv = v ^ v; // v * v^T

				M += (((I * v.inner_product(v)) - vTv) * m_vmasses[n]);
				C += (vTv * m_vmasses[n]);
			}

			return M;
		}


		float t_collision_mesh::calc_moi_scalar_com(const t_vec4f& axis) const {
			float moment_of_inertia = 0.0f;

			for (size_t n = 0; n < m_vertices.size(); n++) {
				const t_vec4f vert = m_vertices[n].to_vector();

				// project vector from origin O to point mass
				// onto axis of rotation which runs through O
				//
				// rv . rv = (v . v) - (v . a)^2
				const float vv = vert.inner_product(vert);
				const float va = vert.inner_product(axis);

				moment_of_inertia += (m_vmasses[n] * (vv - lib_math::square(va)));
			}

			return moment_of_inertia;
		}

		float t_collision_mesh::calc_moi_scalar(const t_pos4f& point) const {
			float moment_of_inertia = 0.0f;

			for (size_t n = 0; n < m_vertices.size(); n++) {
				moment_of_inertia += ((m_vertices[n] - point).sq_len() * m_vmasses[n]);
			}

			return moment_of_inertia;
		}

		float t_collision_mesh::calc_moi_scalar(const t_pos4f& point, const t_vec4f& axis) const {
			float moment_of_inertia = 0.0f;

			for (size_t n = 0; n < m_vertices.size(); n++) {
				// project each vertex onto axis running through <point>
				// use the squared (orthogonal) distance to this axis as
				// weighing factor for the MOI; <axis> is any normalized
				// (unit-length) axis of rotation
				const t_vec4f v = (m_vertices[n] - point);
				const t_vec4f o = v - (axis * axis.inner_product(v));

				moment_of_inertia += (o.sq_len() * m_vmasses[n]);
			}

			return moment_of_inertia;
		}


		float t_collision_mesh::calc_mass() const {
			float vert_mass_sum = 0.0f;

			for (size_t n = 0; n < m_vertices.size(); n++) {
				vert_mass_sum += m_vmasses[n];
			}

			return vert_mass_sum;
		}



		void t_collision_mesh::create_plane_mesh(t_collision_mesh& mesh, const t_vec4f& scales) {
			assert(scales.x() > 0.0f);
			assert(scales.z() > 0.0f);

			// first dispose of any old data
			mesh.clear();

			// corners [indices 0-3], clockwise order
			mesh.add_vertex((t_vec4f(-1.0f,  1.0f, -1.0f) * scales).to_point()); // BL(4)
			mesh.add_vertex((t_vec4f( 1.0f,  1.0f, -1.0f) * scales).to_point()); // BR(5)
			mesh.add_vertex((t_vec4f( 1.0f,  1.0f,  1.0f) * scales).to_point()); // TR(6)
			mesh.add_vertex((t_vec4f(-1.0f,  1.0f,  1.0f) * scales).to_point()); // TL(7)

			mesh.add_normal(t_vec4f::y_axis_vector());
			mesh.add_normal(t_vec4f::y_axis_vector());
			mesh.add_normal(t_vec4f::y_axis_vector());
			mesh.add_normal(t_vec4f::y_axis_vector());

			mesh.add_mass(M_FINF);
			mesh.add_mass(M_FINF);
			mesh.add_mass(M_FINF);
			mesh.add_mass(M_FINF);

			// triangles [indices 0-1]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(0, 1, 2)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(2, 3, 0)));

			mesh.calc_bounding_box();
		}

		// NOTE:
		//   we operate directly on the vertices because non-uniform
		//   scaling would make our object matrices non-orthonormal
		//   (using a scale matrix would require re-calculating all
		//   normals, etc)
		//
		//   furthermore object-space mesh calculations would need
		//   to know the scale factors but mesh_t structure itself
		//   carries no information about them (due to separation
		//   of concerns) so this is the simpler option
		//
		void t_collision_mesh::create_cube_mesh(t_collision_mesh& mesh, const t_vec4f& scales) {
			assert(scales.x() > 0.0f);
			assert(scales.y() > 0.0f);
			assert(scales.z() > 0.0f);
			assert(scales.w() > 0.0f); // mass

			// first dispose of any old data
			mesh.clear();

			// bottom corners [indices 0-3], clockwise order
			mesh.add_vertex((t_vec4f(-1.0f, -1.0f, -1.0f) * scales).to_point()); // BL(0)
			mesh.add_vertex((t_vec4f( 1.0f, -1.0f, -1.0f) * scales).to_point()); // BR(1)
			mesh.add_vertex((t_vec4f( 1.0f, -1.0f,  1.0f) * scales).to_point()); // TR(2)
			mesh.add_vertex((t_vec4f(-1.0f, -1.0f,  1.0f) * scales).to_point()); // TL(3)
			// top corners [indices 4-7], clockwise order
			mesh.add_vertex((t_vec4f(-1.0f,  1.0f, -1.0f) * scales).to_point()); // BL(4)
			mesh.add_vertex((t_vec4f( 1.0f,  1.0f, -1.0f) * scales).to_point()); // BR(5)
			mesh.add_vertex((t_vec4f( 1.0f,  1.0f,  1.0f) * scales).to_point()); // TR(6)
			mesh.add_vertex((t_vec4f(-1.0f,  1.0f,  1.0f) * scales).to_point()); // TL(7)

			for (unsigned int n = 0; n < mesh.get_num_verts(); n++) {
				mesh.add_mass(scales.w() / mesh.get_num_verts());
				mesh.add_normal((mesh.get_vertex(n)).to_vector() / scales);
			}

			// note: all triangles are specified such that
			// vertices have clockwise winding when viewed
			// from outside mesh looking in
			//
			// bottom triangles [indices 0-1]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(2, 1, 0)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(0, 3, 2)));
			// top triangles [indices 2-3]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 5, 6)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(6, 7, 4)));
			// mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 6, 7)));

			// left triangles [indices 4-5]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(7, 3, 0)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(0, 4, 7)));
			// right triangles [indices 6-7]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(1, 2, 6)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(6, 5, 1)));
			// mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(1, 6, 5)));

			// back triangles [indices 8-9]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(0, 1, 5)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(5, 4, 0)));
			// front triangles [indices 10-11]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(6, 2, 3)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(3, 7, 6)));

			// normals are added manually
			// mesh.set_vertex_normals();
			mesh.calc_bounding_box();
		}

		void t_collision_mesh::create_pyramid_mesh(t_collision_mesh& mesh, const t_vec4f& scales) {
			assert(scales.x() > 0.0f);
			assert(scales.y() > 0.0f);
			assert(scales.z() > 0.0f);
			assert(scales.w() > 0.0f); // mass

			// first dispose of any old data
			mesh.clear();

			// bottom corners [indices 0-3], clockwise order
			mesh.add_vertex((t_vec4f(-1.0f, -1.0f, -1.0f) * scales).to_point()); // BL(0)
			mesh.add_vertex((t_vec4f( 1.0f, -1.0f, -1.0f) * scales).to_point()); // BR(1)
			mesh.add_vertex((t_vec4f( 1.0f, -1.0f,  1.0f) * scales).to_point()); // TR(2)
			mesh.add_vertex((t_vec4f(-1.0f, -1.0f,  1.0f) * scales).to_point()); // TL(3)
			// top corner [index 4]
			mesh.add_vertex((t_vec4f( 0.0f,  1.0f,  0.0f) * scales).to_point());

			for (unsigned int n = 0; n < mesh.get_num_verts(); n++) {
				mesh.add_mass(scales.w() / mesh.get_num_verts());
			}

			// bottom triangles [indices 0-1]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(2, 1, 0)));
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(0, 3, 2)));

			// left triangle [index 2]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 3, 0)));
			// right triangle [index 3]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 1, 2)));
			// front triangle [index 4]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 2, 3)));
			// back triangle [index 5]
			mesh.add_poly(t_idx_tri(mesh.get_verts(), t_tup4ui(4, 0, 1)));

			mesh.add_normal(((mesh.get_tri(0)).get_normal() + (mesh.get_tri(2)).get_normal() + (mesh.get_tri(5)).get_normal()) / 3.0f);
			mesh.add_normal(((mesh.get_tri(0)).get_normal() + (mesh.get_tri(3)).get_normal() + (mesh.get_tri(5)).get_normal()) / 3.0f);
			mesh.add_normal(((mesh.get_tri(0)).get_normal() + (mesh.get_tri(3)).get_normal() + (mesh.get_tri(4)).get_normal()) / 3.0f);
			mesh.add_normal(((mesh.get_tri(0)).get_normal() + (mesh.get_tri(2)).get_normal() + (mesh.get_tri(4)).get_normal()) / 3.0f);
			mesh.add_normal(
				((mesh.get_tri(2)).get_normal() + (mesh.get_tri(3)).get_normal() +
				 (mesh.get_tri(4)).get_normal() + (mesh.get_tri(5)).get_normal()) / 4.0f
			);

			// normals are added manually
			// mesh.set_vertex_normals();
			mesh.calc_bounding_box();
		}



		bool t_collision_mesh::check_collision(t_coltest_params* tp, unsigned int t_type, unsigned int p_indx, unsigned int q_indx) {
			const t_collision_mesh* p_mesh = tp->m_meshes[p_indx];
			const t_collision_mesh* q_mesh = tp->m_meshes[q_indx];
			const t_mat44f* p_matr = tp->m_matrices[p_indx];
			const t_mat44f* q_matr = tp->m_matrices[q_indx];

			// output variables
			t_tup4ui& p_sep_inds = tp->m_sep_inds[0];
			t_tup4ui& q_sep_inds = tp->m_sep_inds[1];
			t_vec4f&  p_sep_vect = tp->m_sep_vecs[0];
			t_vec4f&  q_sep_vect = tp->m_sep_vecs[1];
			t_vec4f&  p_sep_dsts = tp->m_sep_dsts[0];
			t_vec4f&  q_sep_dsts = tp->m_sep_dsts[1];

			// type.y() is only relevant for SAT_INNER
			// for the (p=0,q=1) test, p_indx=0 and q_indx=1
			// for the (q=0,p=1) test, p_indx=1 and q_indx=0
			tp->m_test_type.y() = p_indx * (t_type == COL_TEST_TYPE_SAT_INNER);
			tp->m_test_type.x() = t_type;

			assert(p_mesh != nullptr);
			assert(q_mesh != nullptr);
			assert(p_matr != nullptr);
			assert(p_matr != nullptr);

			switch (t_type) {
				case COL_TEST_TYPE_SAT_INNER: {
					return (check_collision_inner(*p_mesh, *q_mesh, *p_matr, *q_matr,  p_sep_inds, q_sep_inds, p_sep_vect, q_sep_vect, p_sep_dsts, q_sep_dsts));
				} break;
				case COL_TEST_TYPE_SAT_OUTER: {
					return (check_collision_outer(*p_mesh, *q_mesh, *p_matr, *q_matr,  p_sep_inds, q_sep_inds, p_sep_vect, q_sep_vect, p_sep_dsts, q_sep_dsts));
				} break;
				case COL_TEST_TYPE_BOUND_BOX: {
					const t_bounding_box& p_bbox = p_mesh->get_bounding_box();
					const t_bounding_box& q_bbox = q_mesh->get_bounding_box();
					const t_vec4f& pq_vec = q_matr->get_t_vector() - p_matr->get_t_vector();

					return (pq_vec.sq_len() <= lib_math::square(p_bbox.get_radius() + q_bbox.get_radius()));
				} break;
				default: {
				} break;
			}

			return false;
		}

		bool t_collision_mesh::check_collision_inner(
			const t_collision_mesh& p_mesh,
			const t_collision_mesh& q_mesh,
			const t_mat44f& p_matr,
			const t_mat44f& q_matr,
			t_tup4ui& p_sep_inds,
			t_tup4ui& q_sep_inds,
			t_vec4f&  p_sep_vect,
			t_vec4f&  q_sep_vect,
			t_vec4f&  p_sep_dsts,
			t_vec4f&  q_sep_dsts
		) {
			if (p_sep_vect != t_vec4f::zero_vector()) {
				assert(p_sep_inds.x() != -1u);

				const t_idx_tri& mesh_tri_p = p_mesh.get_tri(p_sep_inds.x());
				const t_tup4ui&  tri_inds_p = mesh_tri_p.get_indices();
				const t_pos4f&   tri_vert_p = p_matr * p_mesh.get_vertex(tri_inds_p.x());
				const t_vec4f&   tri_norm_p = p_sep_vect; // p_matr * mesh_tri_p.get_normal();

				bool have_sep_axis = true;

				p_sep_dsts.x() = M_FINF;
				p_sep_dsts.y() = M_FINF;

				// project q's vertices only on the given axis (a previous sep_axis)
				// if the projected distances are all positive then p and q are still not
				// intersecting (and the minimum positive distance gives their separation)
				// if at least one is negative, they have since collided (and the maximum
				// negative distance gives their interpenetration)
				for (unsigned int k = 0; k < q_mesh.get_num_verts() && have_sep_axis; k++) {
					const t_pos4f ws_vert_q = q_matr * q_mesh.get_vertex(k);
					const t_vec4f ws_vect_q = ws_vert_q - tri_vert_p;

					p_sep_dsts.y() = ws_vect_q.inner_product(tri_norm_p);
					q_sep_dsts.y() = p_sep_dsts.y();
					p_sep_dsts.x() = std::min(p_sep_dsts.y(), p_sep_dsts.x());
					q_sep_dsts.x() = p_sep_dsts.x();

					// note: GT instead of GEQ such that object faces are inclusive
					// using a positive epsilon makes object stacks much more stable
					have_sep_axis &= (p_sep_dsts.y() > (0.0f * t_tup4f::eps_scalar()));
				}

				// clamps
				p_sep_dsts.x() = std::max(p_sep_dsts.x(), 0.0f);
				q_sep_dsts.x() = std::max(q_sep_dsts.x(), 0.0f);

				return (!have_sep_axis);
			}


			// test all faces of m1 (p) as candidate planes; always
			// use vertex 0 (inds.x()) as representative of each tri
			for (unsigned int n = 0; n < p_mesh.get_num_polys(); n++) {
				const t_idx_tri& mesh_tri_p = p_mesh.get_tri(n);
				const t_tup4ui&  tri_inds_p = mesh_tri_p.get_indices();
				const t_pos4f&   tri_vert_p = p_matr * p_mesh.get_vertex(tri_inds_p.x());
				const t_vec4f&   tri_norm_p = p_matr * mesh_tri_p.get_normal();

				bool have_sep_axis = true;

				// {.x := minimum, .y := current} distance of any of q's vertices to p's face
				p_sep_dsts.x() = M_FINF;
				p_sep_dsts.y() = M_FINF;
				q_sep_dsts.x() = M_FINF;
				q_sep_dsts.y() = M_FINF;

				// project each vertex from m2 (q) onto m1's face-normals
				//
				// no need to normalize vectors, we care only about sign
				// as soon as we find a separating axis (all vertices of
				// m1 in front of face --> positive distances), the test
				// can be stopped
				//
				for (unsigned int k = 0; k < q_mesh.get_num_verts() && have_sep_axis; k++) {
					const t_pos4f ws_vert_q = q_matr * q_mesh.get_vertex(k);
					const t_vec4f ws_vect_q = ws_vert_q - tri_vert_p;

					p_sep_dsts.y() = ws_vect_q.inner_product(tri_norm_p);
					q_sep_dsts.y() = p_sep_dsts.y();
					p_sep_dsts.x() = std::min(p_sep_dsts.y(), p_sep_dsts.x());
					q_sep_dsts.x() = p_sep_dsts.x();

					// note: GT instead of GEQ such that object faces are inclusive
					// using a positive epsilon makes object stacks much more stable
					have_sep_axis &= (p_sep_dsts.y() > (0.0f * t_tup4f::eps_scalar()));
				}

				// return early if we concluded there is *no* collision
				// (i.e. if we found a valid axis separating the objects)
				if (!have_sep_axis)
					continue;

				// clamps
				p_sep_dsts.x() = std::max(p_sep_dsts.x(), 0.0f);
				q_sep_dsts.x() = std::max(q_sep_dsts.x(), 0.0f);

				p_sep_inds = t_tup4ui(n,  -1u, -1u,  0);
				q_sep_inds = t_tup4ui(n,  -1u, -1u,  0);

				p_sep_vect =  t_vec4f(tri_norm_p.x(), tri_norm_p.y(), tri_norm_p.z());
				q_sep_vect = -t_vec4f(tri_norm_p.x(), tri_norm_p.y(), tri_norm_p.z());
				return false;
			}

			// collision, no valid axis
			p_sep_inds = t_tup4ui(-1u, -1u, -1u, -1u);
			q_sep_inds = t_tup4ui(-1u, -1u, -1u, -1u);

			p_sep_vect = t_vec4f::zero_vector();
			q_sep_vect = t_vec4f::zero_vector();
			return true;
		}

		bool t_collision_mesh::check_collision_outer(
			const t_collision_mesh& p_mesh,
			const t_collision_mesh& q_mesh,
			const t_mat44f& p_matr,
			const t_mat44f& q_matr,
			t_tup4ui& p_sep_inds,
			t_tup4ui& q_sep_inds,
			t_vec4f&  p_sep_vect,
			t_vec4f&  q_sep_vect,
			t_vec4f&  p_sep_dsts,
			t_vec4f&  q_sep_dsts
		) {
			// test the cross-product of each face-normal of m1 with
			// each face-normal of m2 as candidate axis and remember
			// the non-colliding combination (if any)
			for (unsigned int n = 0; n < p_mesh.get_num_polys(); n++) {
				const t_idx_tri& mesh_tri_p = p_mesh.get_tri(n);
				const t_vec4f&   tri_norm_p = p_matr * mesh_tri_p.get_normal();
				const t_pos4f&   tri_vert_p = p_matr * p_mesh.get_vertex((mesh_tri_p.get_indices()).x());

				for (unsigned int k = 0; k < q_mesh.get_num_polys(); k++) {
					const t_idx_tri& mesh_tri_q = q_mesh.get_tri(k);
					const   t_vec4f& tri_norm_q = q_matr * mesh_tri_q.get_normal();
					const   t_pos4f& tri_vert_q = q_matr * q_mesh.get_vertex((mesh_tri_q.get_indices()).x());

					const t_vec4f sep_axis = (tri_norm_p.outer_product(tri_norm_q)).normalize();

					int sign_sum_p = 0;
					int sign_sum_q = 0;

					bool have_sep_axis = false;

					// {.x := minimum, .y := current} distance of any of {p,q}'s vertices to {q,p}'s face
					p_sep_dsts.x() = M_FINF;
					p_sep_dsts.y() = M_FINF;
					q_sep_dsts.x() = M_FINF;
					q_sep_dsts.y() = M_FINF;

					// test if the axis is non-degenerate (orthogonal faces)
					// NOTE:
					//   with two axis-aligned cubes this is always true, however we
					//   should not even get to this test-stage if they are separated
					if (sep_axis.sq_len() < 1.0f)
						continue;

					// project a's vertices onto candidate axis
					// if candidate passes, abs(sign_sum_p) will
					// equal |p_mesh_verts|
					for (unsigned int v = 0; v < p_mesh.get_num_verts(); v++) {
						const t_pos4f ws_vert_p = p_matr * p_mesh.get_vertex(v);
						const t_vec4f ws_vect_p = ws_vert_p - tri_vert_q;

						p_sep_dsts.y() = ws_vect_p.inner_product(sep_axis); // or tri_norm_p?
						p_sep_dsts.x() = std::min(p_sep_dsts.y(), p_sep_dsts.x());

						sign_sum_p += sign(p_sep_dsts.y());
					}

					// project b's vertices onto candidate axis
					// if candidate passes, abs(sign_sum_q) will
					// equal |q_mesh_verts|
					for (unsigned int v = 0; v < q_mesh.get_num_verts(); v++) {
						const t_pos4f ws_vert_q = q_matr * q_mesh.get_vertex(v);
						const t_vec4f ws_vect_q = ws_vert_q - tri_vert_p;

						q_sep_dsts.y() = ws_vect_q.inner_product(sep_axis); // or tri_norm_q?
						q_sep_dsts.x() = std::min(q_sep_dsts.y(), q_sep_dsts.x());

						sign_sum_q += sign(q_sep_dsts.y());
					}

					// all of a's vertices must lie on one side, all of b's on opposite side
					// signs must be non-equal or axis does not represent a separating plane
					// it can not be assumed that both objects have meshes with equal vertex
					// counts
					have_sep_axis |= (static_cast<unsigned int>(std::abs(sign_sum_p)) == p_mesh.get_num_verts());
					have_sep_axis |= (static_cast<unsigned int>(std::abs(sign_sum_q)) == q_mesh.get_num_verts());
					have_sep_axis |= (sign(sign_sum_p) != sign(sign_sum_q));

					if (!have_sep_axis)
						continue;

					// clamps
					p_sep_dsts.x() = std::max(p_sep_dsts.x(), 0.0f);
					q_sep_dsts.x() = std::max(q_sep_dsts.x(), 0.0f);

					p_sep_inds = t_tup4ui(n,  -1u, -1u,  1);
					q_sep_inds = t_tup4ui(k,  -1u, -1u,  1);

					p_sep_vect =  sep_axis; // or tri_norm_p?
					q_sep_vect = -sep_axis; // or tri_norm_q?
					return false;
				}
			}

			// collision, no valid axis
			p_sep_inds = t_tup4ui(-1u, -1u, -1u, -1u);
			q_sep_inds = t_tup4ui(-1u, -1u, -1u, -1u);

			p_sep_vect = t_vec4f::zero_vector();
			q_sep_vect = t_vec4f::zero_vector();
			return true;
		}



		template<typename contact_vertex_type>
		static bool add_contact(std::vector<contact_vertex_type>& cvs, const contact_vertex_type& cv, bool forced) {
			// NOTE:
			//   set::find and set::insert might not agree about finding duplicates
			//   depends on definition of t_contact_vertex::operator< (very tricky)
			//
			//   any contact might be a duplicate for the collider but not for the
			//   collidee (or vice versa) which means we should always insert, can
			//   only affect the force-distribution (but this will create an insane
			//   amount of contacts)
			//
			//   when two cube-meshes are touching face-to-face, there are twelve
			//   contacts: three unique <point, normal> pairs per corner per mesh
			//   pairs whose normals are not anti-parallel will be filtered out so
			//   they are never added to the final collection
			//
			//   duplicate-testing costs O(1) + O(2) + ... + O(N) == O(N^2)
			//   complex meshes would need spatial binning to speed this up
			if (forced || std::find(cvs.begin(), cvs.end(), cv) == cvs.end()) {
				cvs.push_back(cv);
				return true;
			}

			return false;
		}


		unsigned int t_collision_mesh::calc_contacts_wspace(
			const t_collision_mesh& collider_mesh,
			const t_collision_mesh& collidee_mesh,
			const t_mat44f& collider_mspace_to_wspace_mat,
			const t_mat44f& collidee_mspace_to_wspace_mat,
			std::vector<t_contact_vertex>& collider_contacts,
			std::vector<t_contact_vertex>& collidee_contacts
		) {
			t_scoped_timer timer(__PRETTY_FUNCTION__);

			const std::vector<t_pos4f>& collider_verts = collider_mesh.get_verts();
			const std::vector<t_pos4f>& collidee_verts = collidee_mesh.get_verts();

			std::vector<t_pos4f> collider_tri_int_points;
			std::vector<t_pos4f> collidee_tri_int_points;

			#if 0
			const t_mat44f collider_wspace_to_mspace_mat = collider_mspace_to_wspace_mat.invert_affine();
			const t_mat44f collidee_wspace_to_mspace_mat = collidee_mspace_to_wspace_mat.invert_affine();
			#endif

			collider_contacts.clear();
			collidee_contacts.clear();
			collider_contacts.reserve(collider_mesh.get_num_polys() * 3);
			collidee_contacts.reserve(collidee_mesh.get_num_polys() * 3);

			collider_tri_int_points.resize(t_raw_tri::get_max_intersect_points());
			collidee_tri_int_points.resize(t_raw_tri::get_max_intersect_points());

			// collision meshes are all convex so they can only ever be in
			// actual contact 1) in a point, 2) along an edge, or 3) along
			// a planar surface (assuming interpenetrations are prevented)
			//
			// note: contact-points are NOT always shared (e.g. consider the
			// tip of a pyramid vs. the face of a cube --> pyramid will have
			// point of contact but cube will not), so we calculate the face
			// intersections *bi*-directionally which should result in an
			// an equal number for both
			//
			// matrix_p_to_q * X means X transformed from p's local space to q's space
			// matrix_q_to_p * X means X transformed from q's local space to p's space
			//
			// note: transforming both members of every pair to world-space
			// is less efficient than testing all triangles of B in A's space
			// (or v.v.) but we want the contact data in WS
			//
			// for each collider-mesh triangle
			for (unsigned int n = 0; n < collider_mesh.get_num_polys(); n++) {
				const t_idx_tri& collider_idx_tri = collider_mesh.get_tri(n);
				const t_raw_tri& collider_raw_tri = collider_idx_tri.to_raw_triangle(&collider_verts[0]);
				const t_raw_tri  collider_rws_tri = collider_raw_tri.transform(collider_mspace_to_wspace_mat);
				const t_vec4f&   collider_tri_nrm = collider_rws_tri.get_normal();

				// for each collidee-mesh triangle
				for (unsigned int k = 0; k < collidee_mesh.get_num_polys(); k++) {
					const t_idx_tri& collidee_idx_tri = collidee_mesh.get_tri(k);
					const t_raw_tri& collidee_raw_tri = collidee_idx_tri.to_raw_triangle(&collidee_verts[0]);
					const t_raw_tri  collidee_rws_tri = collidee_raw_tri.transform(collidee_mspace_to_wspace_mat);
					const t_vec4f&   collidee_tri_nrm = collidee_rws_tri.get_normal();

					unsigned int num_intersect_points = 0;

					#if 1
					// triangle normals should be anti-parallel for "valid" contacts
					// contact-gathering is purely geometrical, unaffected by linear
					// or angular projected velocities
					if (collider_tri_nrm.inner_product(collidee_tri_nrm) >= (0.0f + (t_tup4f::eps_scalar() * 0.0f)))
						continue;
					// always ignore ~exactly parallel triangles
					if (collider_tri_nrm.inner_product(collidee_tri_nrm) >= (1.0f - (t_tup4f::eps_scalar() * 1.0f)))
						continue;
					#endif

					// test triangle-pair (n,k) for intersection
					if ((num_intersect_points = collider_rws_tri.intersect_triangle(collidee_rws_tri, &collider_tri_int_points[0])) != 0) {
						for (unsigned int t = 0; t < num_intersect_points; t++) {
							assert(collider_tri_int_points[t].w() != 0.0f);
							#if 0
							// a valid contact must be inside BOTH objects
							if (!collider_mesh.point_inside_mesh(collider_wspace_to_mspace_mat * collider_tri_int_points[t]))
								continue;
							if (!collidee_mesh.point_inside_mesh(collidee_wspace_to_mspace_mat * collider_tri_int_points[t]))
								continue;
							#endif

							const t_contact_vertex collider_vert = t_contact_vertex(collider_tri_int_points[t], collider_rws_tri.get_normal());
							const t_contact_vertex collidee_vert = t_contact_vertex(collider_tri_int_points[t], collidee_rws_tri.get_normal());

							const bool b0 = add_contact(collider_contacts, collider_vert, false);
							const bool b1 = add_contact(collidee_contacts, collidee_vert, false);

							(b0 && !b1 && add_contact(collidee_contacts, collidee_vert, true));
							(b1 && !b0 && add_contact(collider_contacts, collider_vert, true));
						}
					}

					// test triangle-pair (n,k) for reverse-intersection
					if ((num_intersect_points = collidee_rws_tri.intersect_triangle(collider_rws_tri, &collidee_tri_int_points[0])) != 0) {
						for (unsigned int t = 0; t < num_intersect_points; t++) {
							assert(collidee_tri_int_points[t].w() != 0.0f);
							#if 0
							if (!collider_mesh.point_inside_mesh(collider_wspace_to_mspace_mat * collidee_tri_int_points[t]))
								continue;
							if (!collidee_mesh.point_inside_mesh(collidee_wspace_to_mspace_mat * collidee_tri_int_points[t]))
								continue;
							#endif

							const t_contact_vertex collider_vert = t_contact_vertex(collidee_tri_int_points[t], collider_rws_tri.get_normal());
							const t_contact_vertex collidee_vert = t_contact_vertex(collidee_tri_int_points[t], collidee_rws_tri.get_normal());

							const bool b0 = add_contact(collider_contacts, collider_vert, false);
							const bool b1 = add_contact(collidee_contacts, collidee_vert, false);

							(b0 && !b1 && add_contact(collidee_contacts, collidee_vert, true));
							(b1 && !b0 && add_contact(collider_contacts, collider_vert, true));
						}
					}
				}
			}

			// NOTE:
			//   the order and number of {collider, collidee} contacts is
			//   identical by construction, which most code also relies on
			assert(collider_contacts.empty() == collidee_contacts.empty());
			assert(collider_contacts.size() == collidee_contacts.size());

			return (collider_contacts.size() + collidee_contacts.size());
		}
	};
};


