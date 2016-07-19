#include "../global_defs.hpp"
#include "planar_body.hpp"
#include "rigid_body.hpp"
#include "contact_data.hpp"
#include "../math/template_funcs.hpp"
#include "../misc/scoped_timer.hpp"

namespace newtonics {
	namespace lib_physics {
		bool t_rigid_body::check_collision(const t_rigid_body& collider, const t_planar_body& collidee) {
			return (t_planar_body::check_collision(collidee, collider));
		}

		// implements the separating-axes theorem
		//   given two convex objects (collision meshes) m1 and m2, m1 does
		//   not collide with m2 (and v.v.) if and *only* if all vertices of
		//   m1 lie on one side of some separating plane p and all vertices
		//   of m2 lie on the opposite side
		bool t_rigid_body::check_collision(const t_rigid_body& p, const t_rigid_body& q, t_coltest_params* tp) {
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);
			t_coltest_params test_params;

			const t_collision_mesh& p_mesh = p.get_collision_mesh();
			const t_collision_mesh& q_mesh = q.get_collision_mesh();
			const t_mat44f& p_mat = p.get_transform_mat();
			const t_mat44f& q_mat = q.get_transform_mat();

			if (p_mesh.get_num_polys() == 0)
				return false;
			if (q_mesh.get_num_polys() == 0)
				return false;

			// locally point somewhere non-null
			if (tp == nullptr)
				tp = &test_params;

			// initialize parameters
			tp->m_meshes[0] = &p_mesh;
			tp->m_meshes[1] = &q_mesh;
			tp->m_matrices[0] = &p_mat;
			tp->m_matrices[1] = &q_mat;

			// <*_sep_inds> stores index of face(s) of object(s) found to be separating-plane(s)
			//   if SAT_INNER(p,q) call returns false, p_sep_indices.x() holds the index of a's separating face
			//   if SAT_INNER(q,p) call returns false, q_sep_indices.x() holds the index of b's separating face
			//   if SAT_OUTER(p,q) call returns false, p_sep_indices.x() and q_sep_indices.x() hold the indices
			//   the w-component of *_sep_inds holds the helper-test type (outer_products = false || true)
			//
			// a sub-test returning true means it could not find a separating
			// axis, but that might be a false positive and we need to run at
			// least two non-crossproduct cases since objects can be oriented
			// in such a way that test(A,B) fails but test(B,A) succeeds
			//
			// nothing to do if bounding-box spheres do not overlap
			// note: in this case we do *not* get a separation axis
			// (caller should set bits.x() to zero to calculate one)
			// TODO: test the *boxes* with SAT, 6 quads vs 12 tris
			//
			if (tp->m_test_bits.x() && !t_collision_mesh::check_collision(tp, COL_TEST_TYPE_BOUND_BOX, 0, 1)) return false;
			if (tp->m_test_bits.y() && !t_collision_mesh::check_collision(tp, COL_TEST_TYPE_SAT_INNER, 0, 1)) return false; // (p,q)
			if (tp->m_test_bits.z() && !t_collision_mesh::check_collision(tp, COL_TEST_TYPE_SAT_INNER, 1, 0)) return false; // (q,p)
			if (tp->m_test_bits.w() && !t_collision_mesh::check_collision(tp, COL_TEST_TYPE_SAT_OUTER, 0, 1)) return false;

			// objects are colliding
			return true;
		}


		// this should only be called immediately after check_collision() returns false
		t_vec4f t_rigid_body::get_separating_axis(const t_rigid_body& p, const t_rigid_body& q, const t_coltest_params& tp) {
			#if 0
			const t_collision_mesh* meshes[2] = {&p.get_collision_mesh(), &q.get_collision_mesh()};
			const t_mat44f*          matrs[2] = {&p.get_transform_mat(), &q.get_transform_mat()};
			#else
			(void) p;
			(void) q;
			#endif

			switch (tp.m_test_type.x()) {
				case COL_TEST_TYPE_SAT_INNER: {
					const unsigned int obj_indx = tp.m_test_type.y();
					const unsigned int tri_indx = tp.m_sep_inds[obj_indx].x();

					assert(obj_indx != -1u);
					assert(tri_indx != -1u);

					#if 0
					const t_collision_mesh& mesh = *meshes[obj_indx];
					const t_mat44f&         matr =  *matrs[obj_indx];

					// sat_inner(p_idx=0, q_idx=1) or sat_inner(p_idx=1, q_idx=0) succeeded
					const t_idx_tri& tri = mesh.get_tri(tri_indx);
					const t_vec4f& nrm = matr * tri.get_normal();
					return nrm;
					#else
					return (tp.m_sep_vecs[obj_indx]);
					#endif
				} break;
				case COL_TEST_TYPE_SAT_OUTER: {
					// sat_outer(p_idx=0, q_idx=1) succeeded
					assert(false);
				} break;
				case COL_TEST_TYPE_BOUND_BOX: {
					// box-check succeeded, no axis (!)
				} break;
				default: {
				} break;
			}

			return (t_vec4f::zero_vector());
		}


		std::pair<float, float> t_rigid_body::calc_collision_time_rev(t_rigid_body& p, t_rigid_body& q) {
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

			if (!check_collision(p, q))
				return (std::make_pair(-1.0f, -1.0f));

			const t_pos4f p_cur_pos = p.get_position();
			const t_pos4f q_cur_pos = q.get_position();

			// determine how far both objects moved during the last frame
			// (due to collision-handling we can not use velocities here)
			//
			// note: we also want to (be able to) rewind orientation which
			// requires interpolating between matrices (or the translation
			// and rotation components separately, needs quaternions)
			const t_vec4f p_dif_pos = p.get_dif_position();
			const t_vec4f q_dif_pos = q.get_dif_position();

			// v.x() := lower bound, v.y() := upper bound, v.z() := current
			t_vec4f t_value = t_vec4f(0.0f, 1.0f, 0.0f, 0.0f);
			// m.x() := move object p, m.y() := move object q
			t_vec4f t_mults = t_vec4f(1.0f, 1.0f, 0.0f, 0.0f);


			// objects are colliding, but neither actually moved
			if ((p_dif_pos - t_vec4f::zero_vector()).sq_len() < M_FEPS && (q_dif_pos - t_vec4f::zero_vector()).sq_len() < M_FEPS)
				return (std::make_pair(0.0f, 0.0f));


			// objects are in contact but move as a single body
			// in this case only one object should be moved, or
			// the interpenetration will not be resolved
			//
			// first figure out which object can be moved back
			// this is not guaranteed to be correct (at higher
			// object velocities)
			//
			// if we already have contacts leave objects barely
			// non-touching, otherwise leave them barely touching
			// (within numerical accuracy limits) or calculate
			// contacts between iterations
			if ((p_dif_pos - q_dif_pos).sq_len() < M_FEPS) {
				p.set_position(p_cur_pos - p_dif_pos * 1.0f);

				if (check_collision(p, q)) {
					// still colliding --> only q must be moved
					t_mults.x() = 0.0f;
				} else {
					// no longer colliding --> only p must be moved
					t_mults.y() = 0.0f;
				}

				p.set_position(p_cur_pos);

				#if (ENABLE_RIGIDBODY_DEBUG_PRINTS == 1)
				printf("[%s][p=%u][q=%d][comoving] |p_vel|=%f |q_vel|=%f\n", __func__, p.get_id(), q.get_id(), p_dif_pos.sq_len(), q_dif_pos.sq_len());
				#endif
			}


			while (t_value.z() != (t_value.x() + t_value.y()) * 0.5f) {
				t_value.z() = (t_value.x() + t_value.y()) * 0.5f;

				p.set_position(p_cur_pos - (p_dif_pos * t_value.z() * t_mults.x()));
				q.set_position(q_cur_pos - (q_dif_pos * t_value.z() * t_mults.y()));

				if (check_collision(p, q)) {
					// not far enough forward or backward in time, adjust t-min
					t_value.x() = t_value.z();
				} else {
					// too far forward or backward in time, adjust t-max
					t_value.y() = t_value.z();
				}
			}

			// this fails numerically, but is maybe not a requirement
			// assert(!check_collision(p, q));

			p.set_position(p_cur_pos);
			q.set_position(q_cur_pos);

			return (std::make_pair(t_value.z() * t_mults.x(), t_value.z() * t_mults.y()));
		}

		bool t_rigid_body::handle_clipping_objects(
			t_rigid_body& p,
			t_rigid_body& q,
			t_vec4f& p_sep_axis,
			t_vec4f& q_sep_axis,
			std::vector<t_contact_vertex>& pq_contacts_ws,
			std::vector<t_contact_vertex>& qp_contacts_ws
		) {
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

			// we might have processed this pair already during
			// a previous iteration this frame, *or* p or q are
			// also part of other pairs that were handled
			if (!check_collision(p, q))
				return false;

			const t_collision_mesh& p_mesh = p.get_collision_mesh();
			const t_collision_mesh& q_mesh = q.get_collision_mesh();

			const t_pos4f p_pos = p.get_position();
			const t_pos4f q_pos = q.get_position();

			t_coltest_params tp;

			t_vec4f  p_mov_vec = -p_sep_axis;
			t_vec4f  q_mov_vec = -q_sep_axis;
			t_vec4f pq_mov_mul = {0.0f, 0.0f, 0.0f, 0.0f};

			if (!pq_contacts_ws.empty()) {
				p_mov_vec *= 0.0f;
				q_mov_vec *= 0.0f;

				for (size_t n = 0; n < pq_contacts_ws.size(); n++) {
					p_mov_vec -= pq_contacts_ws[n].get_normal();
					q_mov_vec -= qp_contacts_ws[n].get_normal();
				}

				p_mov_vec.normalize_ref();
				q_mov_vec.normalize_ref();

				// contacts generated during interpenetration are less usable, clear
				// this assumes p and q are already interpenetrating which might not
				// be the case (and no contacts also means no reaction forces, etc)
				pq_contacts_ws.clear();
				qp_contacts_ws.clear();
			}

			{
				assert(!p.ignore_lin_motion() || q.ignore_lin_motion());
				assert(p_mov_vec != t_vec4f::zero_vector());
				assert(q_mov_vec != t_vec4f::zero_vector());

				unsigned int num_iters = 0;

				// first push objects apart along their separation-axes (which
				// must exist, assuming the pair did not interpenetrate during
				// the previous frame) by doubling the step-size each iteration
				// then pull them back together so calc_contacts() has a better
				// input state
				// this a minor hack in lieu of more robust contact-calculation
				tp.m_test_bits.x() = 0;
				tp.m_test_bits.w() = 0;

				do {
					pq_mov_mul.z() = pq_mov_mul.x();
					pq_mov_mul.w() = pq_mov_mul.y();
					pq_mov_mul.x() = (1 << num_iters) * 0.125f;
					pq_mov_mul.y() = (1 << num_iters) * 0.125f;

					p.set_position(p_pos + p_mov_vec * pq_mov_mul.x());
					q.set_position(q_pos + q_mov_vec * pq_mov_mul.y());

					num_iters += 1;
				} while (check_collision(p, q, &tp));
			}

			// the axes are currently mirrored (p=-q) so either both are valid
			// or neither is; in the latter case we need to calculate them OTF
			// during/after separation
			// note: this executes "*axis = tp.axis", not "axis = &tp.axis";
			if (p_sep_axis == t_vec4f::zero_vector()) {
				p_sep_axis =  get_separating_axis(p, q, tp);
				q_sep_axis = -get_separating_axis(p, q, tp);

				assert(tp.m_test_bits.x() == 1 || p_sep_axis != t_vec4f::zero_vector());
				assert(tp.m_test_bits.x() == 1 || q_sep_axis != t_vec4f::zero_vector());
			}

			{
				// enable the early-out again, use axes to speed up next calls
				tp.m_test_bits.x() = 1;

				// search the interval between *_mov_mul_max and *_mov_mul_min
				// (final pq_mov_mul could also be returned, currently unused)
				//
				// while (std::fabs(pq_mov_mul.x() - pq_mov_mul.z()) > (1.0f / 65536.0f)) {
				for (unsigned int n = 0; n < 32; n++) {
					const float p_mov_mul_mid = (pq_mov_mul.x() + pq_mov_mul.z()) * 0.5f;
					const float q_mov_mul_mid = (pq_mov_mul.y() + pq_mov_mul.w()) * 0.5f;

					p.set_position(p_pos + p_mov_vec * p_mov_mul_mid);
					q.set_position(q_pos + q_mov_vec * q_mov_mul_mid);

					if (check_collision(p, q, &tp)) {
						// back too far
						pq_mov_mul.z() = p_mov_mul_mid;
						pq_mov_mul.w() = q_mov_mul_mid;
					} else {
						// not far enough
						pq_mov_mul.x() = p_mov_mul_mid;
						pq_mov_mul.y() = q_mov_mul_mid;
					}
				}
			}

			#if 0
			// final epsilon-pull; should *always* move objects closer
			// instead of this, find distance along p_sep_axis and use
			// that for a single final step
			for (unsigned int n = 0; n < 32 && !check_collision(p, q, &tp); n++) {
				p.set_position(p.get_position() - p_mov_vec * (1.0f / 8192.0f));
				q.set_position(q.get_position() - q_mov_vec * (1.0f / 8192.0f));
			}

			#else

			check_collision(p, q, &tp);

			p.set_position(p.get_position() - p_sep_axis * tp.m_sep_dsts[ (tp.m_test_type.y()    ) & 1 ].x());
			q.set_position(q.get_position() - q_sep_axis * tp.m_sep_dsts[ (tp.m_test_type.y() + 1) & 1 ].x());
			#endif

			{
				const t_mat44f& p_matrix = p.get_transform_mat();
				const t_mat44f& q_matrix = q.get_transform_mat();

				return (t_collision_mesh::calc_contacts_wspace(p_mesh, q_mesh, p_matrix, q_matrix, pq_contacts_ws, qp_contacts_ws) == 0);
			}
		}




		bool t_rigid_body::apply_reaction_forces(
			t_rigid_body& p,
			t_rigid_body& q,
			const std::vector<t_contact_vertex>& pq_contacts_ws,
			const std::vector<t_contact_vertex>& qp_contacts_ws
		) {
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

			// get the current linear (TODO: angular) forces acting on p and q
			const t_vec4f p_lf = p.get_lin_force();
			const t_vec4f q_lf = q.get_lin_force();
			const t_vec4f   ev = t_vec4f::eps_vector() * 10.0f;

			// forces and distribs [0] and [2] correspond to actions
			//   (a) p exerts a force p_lf on q  (can be zero)
			//   (b) q exerts a force q_lf on p  (can be zero)
			// forces and distribs [1] and [3] correspond to reactions
			//   (c) q reacts to (a) by exerting force -p_lf on p
			//   (d) p reacts to (b) by exerting force -q_lf on q
			// vertical stacks should remain stable just by exchanging these
			const std::array<t_vec4f, 2 + 2> force_vectors = {q_lf, -p_lf,  p_lf, -q_lf};

			p.calc_contact_force_distrib(force_vectors[0], pq_contacts_ws, FORCE_TYPE_REAC); // (b)
			p.calc_contact_force_distrib(force_vectors[1], pq_contacts_ws, FORCE_TYPE_REAC); // (c)
			q.calc_contact_force_distrib(force_vectors[2], qp_contacts_ws, FORCE_TYPE_REAC); // (a)
			q.calc_contact_force_distrib(force_vectors[3], qp_contacts_ws, FORCE_TYPE_REAC); // (d)

			// these are epsilon-tolerance comparisons: true does not imply zero
			if ((p.get_lin_force()).equals(t_vec4f::zero_vector(), ev)) {
				q.add_lin_force(p.get_lin_force());
				p.set_lin_force(t_vec4f::zero_vector());
			}

			if ((q.get_lin_force()).equals(t_vec4f::zero_vector(), ev)) {
				p.add_lin_force(q.get_lin_force());
				q.set_lin_force(t_vec4f::zero_vector());
			}

			return ((!p_lf.equals(p.get_lin_force(), ev)) || (!q_lf.equals(q.get_lin_force(), ev)));
		}


		bool t_rigid_body::apply_collision_forces(
			t_rigid_body& p,
			t_rigid_body& q,
			const std::vector<t_contact_vertex>& pq_contacts_ws,
			const std::vector<t_contact_vertex>& qp_contacts_ws
		) {
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

			// for each matching pair of contacts, apply collision response
			// at least one of the objects needs to be moving, otherwise no
			// momentum will be exchanged
			//
			// applying forces through corner contacts will cause both parties
			// to receive smaller acceleration / decceleration than they would
			// with COM-aligned forces, so multiple collisions will occur over
			// time until velocities have equalized --> iterate until momentum
			// exchange stops
			std::vector<t_vec4f> pq_contact_weights = std::move(p.calc_contact_force_weights(pq_contacts_ws));
			std::vector<t_vec4f> qp_contact_weights = std::move(q.calc_contact_force_weights(qp_contacts_ws));

			std::vector< std::pair<t_vec4f, t_vec4f> > collision_forces(pq_contacts_ws.size());

			const t_mat44f& p_raw_mat = p.get_transform_mat();
			const t_mat44f& q_raw_mat = q.get_transform_mat();
			const t_mat44f  p_inv_mat = p_raw_mat.invert_affine();
			const t_mat44f  q_inv_mat = q_raw_mat.invert_affine();

			bool ret = false;

			for (size_t k = 0; k < pq_contacts_ws.size(); k++) {
				collision_forces[k] = calc_collision_forces(p, q, p_inv_mat, q_inv_mat, pq_contacts_ws[k], qp_contacts_ws[k]);

				ret |= (!(collision_forces[k].first ).equals(t_vec4f::zero_vector(), t_vec4f::eps_vector() * 10.0f));
				ret |= (!(collision_forces[k].second).equals(t_vec4f::zero_vector(), t_vec4f::eps_vector() * 10.0f));
			}

			for (size_t k = 0; k < pq_contacts_ws.size(); k++) {
				const t_contact_response& p_contact_resp = p.calc_contact_response(collision_forces[k].first  * pq_contact_weights[k]/*.w()*/, pq_contacts_ws[k], FORCE_TYPE_COLL);
				const t_contact_response& q_contact_resp = q.calc_contact_response(collision_forces[k].second * qp_contact_weights[k]/*.w()*/, qp_contacts_ws[k], FORCE_TYPE_COLL);

				p.add_contact_response(p_contact_resp);
				q.add_contact_response(q_contact_resp);
			}

			// simulate accelerations p and q are going to experience from forces *this* iteration
			p.set_rmoi((p.get_collision_mesh()).calc_moi_scalar(p.get_center_of_mass(), p.calc_rot_axis()));
			q.set_rmoi((q.get_collision_mesh()).calc_moi_scalar(q.get_center_of_mass(), q.calc_rot_axis()));
			p.set_lin_velocity(p.calc_lin_velocity(1.0f), CURR_ATTRIB_IDX);
			p.set_ang_velocity(p.calc_ang_velocity(1.0f), CURR_ATTRIB_IDX);
			q.set_lin_velocity(q.calc_lin_velocity(1.0f), CURR_ATTRIB_IDX);
			q.set_ang_velocity(q.calc_ang_velocity(1.0f), CURR_ATTRIB_IDX);

			// add current forces to temporaries
			p.add_lin_force(p.get_lin_force(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
			p.add_ang_force(p.get_ang_force(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
			q.add_lin_force(q.get_lin_force(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
			q.add_ang_force(q.get_ang_force(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
			// clear current forces
			p.clear_forces(CURR_ATTRIB_IDX);
			q.clear_forces(CURR_ATTRIB_IDX);

			return ret;
		}


		void t_rigid_body::apply_friction_forces(
			t_rigid_body& p,
			t_rigid_body& q,
			const std::vector<t_contact_vertex>& pq_contacts_ws,
			const std::vector<t_contact_vertex>& qp_contacts_ws
		) {
			(void) p;
			(void) q;
			(void) pq_contacts_ws;
			(void) qp_contacts_ws;

			#if 0
			t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

			for (size_t k = 0; k < pq_contacts_ws.size(); k++) {
			}
			for (size_t k = 0; k < qp_contacts_ws.size(); k++) {
			}
			#endif
		}


		std::pair<t_vec4f, t_vec4f> t_rigid_body::calc_collision_forces(
			t_rigid_body& p,
			t_rigid_body& q,
			const t_mat44f& p_inv_mat,
			const t_mat44f& q_inv_mat,
			const t_contact_vertex& p_contact_ws,
			const t_contact_vertex& q_contact_ws
		) {
			std::pair<t_vec4f, t_vec4f> ret;

			// (WS) normal of face *of p* which q collided with
			const t_vec4f& p_contact_normal_ws = p_contact_ws.get_normal();
			// (WS) normal of face *of q* which p collided with
			const t_vec4f& q_contact_normal_ws = q_contact_ws.get_normal();

			// raw absolute linear velocities at this contact-point
			const t_vec4f p_v_raw = p.get_lin_velocity() + p.calc_point_lin_velocity(p_inv_mat * p_contact_ws.get_point());
			const t_vec4f q_v_raw = q.get_lin_velocity() + q.calc_point_lin_velocity(q_inv_mat * q_contact_ws.get_point());

			// project current raw velocity-vectors onto contact-surface normals
			//
			// for 2D and 3D collisions all velocities must be the components
			// perpendicular to the tangent line (or plane) at contact points
			//
			// const t_vec4f p_v = p_v_raw.normalize() * -std::min(0.0f, p_v_raw.inner_product(q_contact_normal_ws));
			// const t_vec4f q_v = q_v_raw.normalize() * -std::min(0.0f, q_v_raw.inner_product(p_contact_normal_ws));
			//
			const t_vec4f p_v = p_v_raw.normalize() * std::fabs(p_v_raw.inner_product(q_contact_normal_ws));
			const t_vec4f q_v = q_v_raw.normalize() * std::fabs(q_v_raw.inner_product(p_contact_normal_ws));

			// turn projected velocities into relative velocities
			//   p_v = v1
			//   q_v = v2
			//   qp_v = v2 - v1 = -v'
			//   pq_v = v1 - v2 =  v'
			const t_vec4f qp_v = q_v - p_v;
			const t_vec4f pq_v = p_v - q_v;

			const float p_mass = p.get_mass();
			const float q_mass = q.get_mass();
			const float s_mass = p_mass + q_mass;

			if (qp_v.inner_product(p_contact_normal_ws) >= 0.0f)
				return ret;
			if (pq_v.inner_product(q_contact_normal_ws) >= 0.0f)
				return ret;

			// bail if total relative momentum is zero
			if ((p_v + q_v) == t_vec4f::zero_vector())
				return ret;


			// EXTREME ELASTIC (OBJECTS BOUNCE AWAY, c_r=1) CASE, RELATIVE VELOCITIES:
			//   v1' = v1-v2 = pq_v
			//   v2' = v2-v2 = 0
			//   u1' = ((m1-m2) / (m1+m2))*v1' [SIC!]
			//   u2' = ((m1+m1) / (m1+m2))*v1' [SIC!]
			// const t_vec4f p_u = pq_v * (p_mass - q_mass) / (s_mass);
			// const t_vec4f q_u = pq_v * (p_mass + p_mass) / (s_mass);
			//
			// EXTREME INELASTIC (OBJECTS STICK TOGETHER, c_r=0) CASE, RELATIVE VELOCITIES:
			//   v1' = v1-v2 = pq_v
			//   v2' = v2-v2 = 0
			//   u'  = (m1*v1')/(m1+m2) [FROM m1*v1' = (m1+m2)*u']
			// const t_vec4f p_u = (pq_v * p_mass) / (s_mass);
			// const t_vec4f q_u = (pq_v * p_mass) / (s_mass);
			//
			// GENERAL {IN}ELASTIC CASE, RELATIVE VELOCITIES
			//   v1' = v1-v2 = pq_v
			//  -v1' = v2-v1 = qp_v
			//   v2' = v2-v2 = 0
			//   u1' = ((m1*v1') + m2*c_r*-v1') / (m1+m2)
			//   u2' = ((m1*v1') + m1*c_r* v1') / (m1+m2)
			//   u1  = ?
			//   u2  = ?
			// const t_vec4f p_u = (pq_v * p_mass + qp_v * q_mass * c_r) / (s_mass);
			// const t_vec4f q_u = (pq_v * p_mass + pq_v * p_mass * c_r) / (s_mass);
			//
			// GENERAL {IN}ELASTIC CASE, ABSOLUTE VELOCITIES
			//   v1' = v1-v2 = pq_v
			//  -v1' = v2-v1 = qp_v
			//   Cr  = (u2-u1) / (v1')
			//   u1  = ((m1*v1 + m2*v2) + m2*Cr*(-v1')) / (m1+m2)
			//   u2  = ((m1*v1 + m2*v2) + m1*Cr*( v1')) / (m1+m2)
			//
			// 'x_v' := pre-collision speed for object x
			// 'x_u' := post-collision speed for object x
			//
			const t_vec4f p_u = ((p_v * p_mass + q_v * q_mass) + (qp_v * q_mass * (q.get_force_coeffs()).x())) / s_mass;
			const t_vec4f q_u = ((p_v * p_mass + q_v * q_mass) + (pq_v * p_mass * (p.get_force_coeffs()).x())) / s_mass;

			const t_vec4f p_post_mom = p_u * p_mass;
			const t_vec4f q_post_mom = q_u * q_mass;

			// m1*v1 + m2*v2 (projected total momentum; must equal m1*u1 + m2*u2)
			//   [m1=1 * v1=0.10, m2=1 * v2=0.00] --> [m1=1 * v1= 0.00, m2=1 * v2=0.10] is a valid solution when c_r=1
			//   [m1=1 * v1=0.10, m2=1 * v2=0.00] --> [m1=1 * v1= 0.05, m2=1 * v2=0.05] is a valid solution when c_r=0
			//   [m1=1 * v1=0.10, m2=1 * v2=0.00] --> [m1=1 * v1=-0.05, m2=1 * v2=0.05] is never a valid solution!
			const t_vec4f pq_pre_mom = p_v * p_mass + q_v * q_mass;

			// figure out forces from the momentum changes
			// a change in momentum wrt. time equals force:
			//   F = m*a = m*dv (a=dv/dt)
			//   P = m*v
			//   d[P] = d[m*v] = m*dv = m*a = F
			// (assuming object mass stays the same)
			const t_vec4f p_force = p_post_mom - (p_v * p_mass);
			const t_vec4f q_force = q_post_mom - (q_v * q_mass);

			// check for momentum conservation
			// assert(std::fabs(pq_pre_mom.sq_len() - (p_post_mom + q_post_mom).sq_len()) < t_tup4f::eps_scalar());
			//
			// this is numerically more stable for large mass-values
			// assert((pq_pre_mom.sq_len() / (p_post_mom + q_post_mom).sq_len()) > (1.0f - t_tup4f::eps_scalar()));
			//
			assert(lib_math::fp_eq(pq_pre_mom.sq_len(), (p_post_mom + q_post_mom).sq_len(), t_tup4f::eps_scalar()));

			ret.first  = p_force;
			ret.second = q_force;
			return ret;
		}


		std::pair<t_vec4f, t_vec4f> t_rigid_body::calc_friction_forces(
			t_rigid_body& p,
			t_rigid_body& q,
			const t_contact_vertex& p_contact_ws,
			const t_contact_vertex& q_contact_ws
		) {
			(void) p;
			(void) q;
			(void) p_contact_ws;
			(void) q_contact_ws;

			// TODO: like planar_body, but for mutual contacts
			std::pair<t_vec4f, t_vec4f> ret;
			return ret;
		}


		// force applied through an arbitrary contact point, assumed
		// to be parallel or orthogonal to normal at the given point
		//
		// this creates an angular acceleration (moment=<lever, scale>)
		// *and* possibly a linear acceleration, which are represented
		// in the t_contact_response structure
		//
		t_contact_response  t_rigid_body::calc_contact_response_raw(
			const t_vec4f& contact_force_ws,
			const t_pos4f& contact_point_ws,
			unsigned int force_type
		) const {
			t_contact_response force_response;

			// vector from point at which force applies *to* object COM
			const t_vec4f  force_moment_arm = (get_transform_mat(CURR_ATTRIB_IDX) * m_center_of_mass) - contact_point_ws;
			const t_vec4f& force_moment_dir = force_moment_arm.normalize();

			(void) force_type;

			// find out if this force is going to cause a rotation
			//
			// |force_moment_arm| is length of the moment arm (arm dot dir)
			// the outer product of force_moment_arm and angular_force will
			// be the rotation axis, scaled by the torque's magnitude
			const float lin_force_scale    = contact_force_ws.inner_product(force_moment_dir);
			const float ang_force_scale    = force_moment_arm.inner_product(force_moment_dir);
			const float raw_force_scale_sq = contact_force_ws.sq_len() - t_tup4f::eps_scalar();

			if (square(lin_force_scale) < raw_force_scale_sq) {
				// linear-component force is parallel to moment-arm
				// angular-component force is orthogonal to moment-arm
				(void) lin_force_scale;
				(void) ang_force_scale;

				#if (ENABLE_ANGULAR_FORCES == 1)
				// this is the proper formulation, but we do not support
				// friction among objects yet so e.g. skewed stacks are
				// not stable (ground contact-responses cause movements
				// along the z-dimension because contact-points are not
				// aligned with COM)
				const t_vec4f lin_force_vector = force_moment_dir * lin_force_scale; // component of F parallel to moment-arm
				const t_vec4f ang_force_vector = contact_force_ws - lin_force_vector; // component of F orthogonal to moment-arm
				#else
				const t_vec4f lin_force_vector = contact_force_ws.normalize() * lin_force_scale;
				const t_vec4f ang_force_vector = contact_force_ws - lin_force_vector;
				#endif

				// assert((lin_force_vector + ang_force_vector) == contact_force_ws);
				// redundant: lin_force_vector is simply force_moment_dir scaled
				// assert(std::fabs(lin_force_vector.inner_product(force_moment_dir)) > (1.0f - t_tup4f::eps_scalar()));
				// this fails numerically for larger object masses
				// assert(std::fabs(ang_force_vector.inner_product(force_moment_dir)) < (0.0f + t_tup4f::eps_scalar()));

				force_response.set_lin_force(lin_force_vector);
				force_response.set_ang_force(ang_force_vector.outer_product(-force_moment_arm));

				// force_response.set_ang_force(ang_force_vector);
				// force_response.set_ang_lever(force_moment_arm);
			} else {
				// contact point coincides with center of mass
				// so force only induces a linear acceleration
				force_response.set_lin_force(contact_force_ws);
			}

			return force_response;
		}

		t_contact_response  t_rigid_body::calc_contact_response(
			const t_vec4f& contact_force_ws,
			const t_contact_vertex& contact_vertex_ws,
			unsigned int force_type
		) const {
			// find magnitude of force projected onto contact-normal
			const float raw_contact_force_scale = contact_force_ws.inner_product(contact_vertex_ws.get_normal());
			// const float abs_contact_force_scale =  std::fabs(raw_contact_force_scale);
			// const float min_contact_force_scale = -std::min(raw_contact_force_scale, 0.0f);

			// use contact-normal to split force into its parallel and orthogonal components
			// (these components themselves can each induce linear and angular accelerations)
			//
			// contact-normals always face outward, forces can be either parallel (dot(N, F) >= 0)
			// or anti-parallel (dot(N, F) < 0) to normal and cause object to be pulled or pushed
			// respectively
			//
			const t_vec4f para_contact_force = contact_vertex_ws.get_normal() * raw_contact_force_scale;
			const t_vec4f orth_contact_force = contact_force_ws - para_contact_force;

			const t_contact_response para_force_resp = calc_contact_response_raw(para_contact_force, contact_vertex_ws.get_point(), force_type);
			const t_contact_response orth_force_resp = calc_contact_response_raw(orth_contact_force, contact_vertex_ws.get_point(), force_type);

			t_contact_response force_response;

			switch (force_type) {
				case FORCE_TYPE_REAC: {
					// no pulling forces for now (objects are not physically linked)
					// can only happen inside a reaction-chain with rest. coeff. of 0
					//
					// NOTE: friction will pull (at some contacts) rather than push
					if (raw_contact_force_scale >= 0.0f)
						return force_response;

					// NOTE: orthogonal responses omitted (horizontal stack, two
					// objects side by side, gravity acting on both would induce
					// rotation if included)
					force_response.set_lin_force(para_force_resp.get_lin_force());
					force_response.set_ang_force(para_force_resp.get_ang_force());
				} break;

				default: {
					// sum up the linear and angular contributions from each component
					// for a single contact this is allowed, since angular forces over
					// the same arm obey  M x (F1 + F2 + ...) = (M x F1) + (M x F2) + ...
					force_response.set_lin_force(para_force_resp.get_lin_force() + orth_force_resp.get_lin_force());
					force_response.set_ang_force(para_force_resp.get_ang_force() + orth_force_resp.get_ang_force());
				} break;
			}

			return force_response;
		}


		std::vector<t_vec4f> t_rigid_body::calc_contact_force_weights(
			const std::vector<t_contact_vertex>& contacts_ws
		) const {
			const size_t k = contacts_ws.size();

			const t_pos4f com_ws_pos = get_transform_mat(CURR_ATTRIB_IDX) * m_center_of_mass;
			      t_vec4f com_ws_vec;

			// last element stores the sum; initially all zero-vectors
			//
			// note: len(sum(x1+x2+...,y1+y2+...,z1+z2+...)) is *not* equal to
			// sum(len(x1,y1,z1)+len(x2,y2,z2)+...) which means vector-summing
			// will produce a different result than scalar-summing
			//
			std::vector<t_vec4f> contact_weights(k + 1);

			for (size_t n = 0; n < k; n++) {
				com_ws_vec = (contacts_ws[n].get_point() - com_ws_pos).abs();

				// xyz-distances are unsigned, use abs
				contact_weights[k]     += (contact_weights[n]     = com_ws_vec.max(t_vec4f::eps_vector()));
				contact_weights[k].w() += (contact_weights[n].w() = com_ws_vec.len(                     ));
			}

			// empty contact-set can occur if ENABLE_FRICTION_FORCES is 0
			assert(contacts_ws.empty() || contact_weights[k].w() != 0.0f);

			// normalize to find each contact's ratio; w-components store |com_vec[n]|/sum_dist
			//
			// a contact represents an (infinitesimal) surface area; any given
			// linear force will act through all of these but not (necessarily)
			// uniformly depending on each contact's position relative to the
			// COM (assume objects have homogeneous uniform mass distribution)
			for (size_t n = 0; n < k; n++) {
				contact_weights[n] = contact_weights[n] / contact_weights.back();
			}

			return contact_weights;
		}

		void t_rigid_body::calc_contact_force_distrib(
			const t_vec4f& contact_force_ws,
			const std::vector<t_contact_vertex>& contacts_ws,
			unsigned int force_type
		) {
			std::vector<t_vec4f> contact_weights = std::move(calc_contact_force_weights(contacts_ws));

			for (size_t n = 0; n < contacts_ws.size(); n++) {
				const t_vec4f nrm_contact_force = contact_force_ws * contact_weights[n]/*.w()*/;
				const t_contact_response nrm_contact_resp = calc_contact_response(nrm_contact_force, contacts_ws[n], force_type);

				add_contact_response(nrm_contact_resp);
			}
		}


		static t_contact_response clamp_force_response(t_contact_response force_response) {
			#if (CLAMP_FORCE_RESPONSES == 1)
			t_vec4f lf = force_response.get_lin_force();
			t_vec4f af = force_response.get_ang_force();

			lf.x() *= (std::fabs(lf.x()) >= (t_tup4f::eps_scalar()));
			lf.y() *= (std::fabs(lf.y()) >= (t_tup4f::eps_scalar()));
			lf.z() *= (std::fabs(lf.z()) >= (t_tup4f::eps_scalar()));
			af.x() *= (std::fabs(af.x()) >= (t_tup4f::eps_scalar()));
			af.y() *= (std::fabs(af.y()) >= (t_tup4f::eps_scalar()));
			af.z() *= (std::fabs(af.z()) >= (t_tup4f::eps_scalar()));

			force_response.set_lin_force(lf);
			force_response.set_ang_force(af);
			#endif
			return force_response;
		}


		void t_rigid_body::add_clamped_contact_force(
			const t_vec4f& contact_force_ws,
			const t_contact_vertex& contact_vertex_ws,
			const unsigned int force_type
		) {
			add_contact_response(clamp_force_response(calc_contact_response(contact_force_ws, contact_vertex_ws, force_type)));
		}

		void t_rigid_body::add_contact_response(const t_contact_response& force_response) {
			add_lin_force(force_response.get_lin_force());
			add_ang_force(force_response.get_ang_force());
		}


		void t_rigid_body::print_state(unsigned int frame_num, const char* hdr_str) const {
			#if (ENABLE_RIGIDBODY_DEBUG_PRINTS == 1)
			const char* fmt_str = "%s[f=%u][id=%u] pos=<%.1f,%.6f,%.1f> dif=<%.1f,%.1f,%.1f>  lv=<%f,%f,%f> av=<%f,%f,%f>  lf=<%f,%f,%f> af=<%f,%f,%f>\n";

			const t_pos4f&  cp = get_position();
			const t_vec4f& dp = get_dif_position();
			const t_vec4f& lv = get_lin_velocity();
			const t_vec4f& av = get_ang_velocity();
			const t_vec4f& lf = get_lin_force();
			const t_vec4f& af = get_ang_force();

			printf(
				fmt_str, hdr_str,
				frame_num, get_id(),
				cp.x(),cp.y(),cp.z(),
				dp.x(),dp.y(),dp.z(),
				lv.x(),lv.y(),lv.z(),
				av.x(),av.y(),av.z(),
				lf.x(),lf.y(),lf.z(),
				af.x(),af.y(),af.z()
			);
			#else
			(void) frame_num;
			#endif
		}
	};
};

