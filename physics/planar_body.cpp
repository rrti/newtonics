#include "../global_defs.hpp"
#include "planar_body.hpp"
#include "rigid_body.hpp"

namespace newtonics {
	namespace lib_physics {
		bool t_planar_body::check_collision(const t_planar_body& collider, const t_rigid_body& collidee) {
			bool ret = false;

			const t_collision_mesh& collidee_mesh = collidee.get_collision_mesh();
			const t_mat44f& collidee_matrix = collidee.get_transform_mat();

			// note: this assumes an infinite plane, actual contacts can be empty
			// collidee mesh surfaces count as inclusive for collision against us
			for (unsigned int n = 0; (n < collidee_mesh.get_num_verts() && !ret); n++) {
				ret |= (collider.get_point_distance(collidee_matrix * collidee_mesh.get_vertex(n)) <= 0.0f);
			}

			return ret;
		}


		bool t_planar_body::handle_collision(
			const t_planar_body& collider,
			      t_rigid_body& collidee,
			std::vector<t_contact_vertex>& collider_contacts,
			std::vector<t_contact_vertex>& collidee_contacts
		) {
			const t_pos4f p = collidee.get_position();

			if (!collider.handle_clipping_object(collider, collidee, collider_contacts, collidee_contacts))
				return false;

			if (calc_contacts_wspace(collider, collidee, collider_contacts, collidee_contacts) != 0) {
				// do a number of passes so the reaction force is propagated
				for (unsigned int n = 0; handle_colliding_object(collider, collidee, collider_contacts, collidee_contacts, n); n++);
			} else {
				// if no contacts after clip-prevention, object is not on
				// the ground proper (i.e. outside the plane's boundaries)
				collidee.set_position(p);
			}

			return true;
		}

		bool t_planar_body::handle_clipping_object(
			const t_planar_body& collider,
			      t_rigid_body& collidee,
			const std::vector<t_contact_vertex>&,
			const std::vector<t_contact_vertex>&
		) {
			const t_collision_mesh& collidee_mesh = collidee.get_collision_mesh();
			const t_mat44f& collider_matrix = collider.get_transform_mat();
			const t_mat44f& collidee_matrix = collidee.get_transform_mat();

			// find the maximally-penetrating distance;
			// might also use the previous g_*_contacts
			float max_clip_dist = t_tuple4f::eps_scalar();

			for (unsigned int n = 0; n < collidee_mesh.get_num_verts(); n++) {
				max_clip_dist = std::min(max_clip_dist, collider.get_point_distance(collidee_matrix * collidee_mesh.get_vertex(n)));
			}

			if (max_clip_dist > 0.0f)
				return false;

			// NOTE:
			//   due to clipping an object can build up more speed than it would
			//   if the collision detection was non-discrete, and re-positioning
			//   lets the object incorrectly keep this momentum
			//   if the extra momentum also cancels the restitution coefficient
			//   reduction, the result is oscillation --> binary-search for time
			//   and do vel -= acc*t and pos -= vel?
			// push the body out by this distance
			//
			collidee.set_position(collidee.get_position() - (collider_matrix.get_y_vector() * max_clip_dist));
			return true;
		}

		bool t_planar_body::handle_colliding_object(
			const t_planar_body& collider,
			      t_rigid_body& collidee,
			const std::vector<t_contact_vertex>& collider_contacts,
			const std::vector<t_contact_vertex>& collidee_contacts,
			unsigned int n
		) {
			bool ret = false;


			const t_vec4f lf = collidee.get_lin_force();
			const t_vec4f cf = collider.calc_collision_force(collidee);
			const t_vec4f rf = collider.calc_reaction_force(collidee);
			const t_vec4f ff = collider.calc_friction_force(collidee);
			const t_vec4f ev = t_vec4f::eps_vector() * 10.0f;

			assert(!collider_contacts.empty());
			assert(!collidee_contacts.empty());

			#if (ENABLE_PLANEBODY_DEBUG_PRINTS == 1)
			collidee.print_state(n, "\t\t[pb::handle_col_ground::iter][pre]");
			#endif

			const std::array<t_vec4f, NUM_FORCE_TYPES> collidee_force_vectors = {cf, rf, ff};

			// only the reaction force needs to be added every iteration
			const unsigned int min_force_t = FORCE_TYPE_REAC * (n > 0);
			const unsigned int max_force_t = FORCE_TYPE_REAC * (n > 0) + (NUM_FORCE_TYPES - 1) * (n == 0);

			for (unsigned int t = min_force_t; t <= max_force_t; t++) {
				collidee.calc_contact_force_distrib(collidee_force_vectors[t], collidee_contacts, t);
			}


			#if (ENABLE_PLANEBODY_DEBUG_PRINTS == 1)
			collidee.print_state(n, "\t\t[pb::handle_col_ground::iter][pst]");
			#endif

			// this is an epsilon-tolerance comparison: true does not imply zero
			if ((collidee.get_lin_force()).equals(t_vec4f::zero_vector(), ev))
				collidee.set_lin_force(t_vec4f::zero_vector());

			ret |= (!lf.equals(collidee.get_lin_force(), ev));
			ret &= (n < MAX_FORCE_EXCHANGE_ITERS);

			return ret;
		}


		unsigned int t_planar_body::calc_contacts_wspace(
			const t_planar_body& collider,
			const t_rigid_body& collidee,
			std::vector<t_contact_vertex>& collider_contacts,
			std::vector<t_contact_vertex>& collidee_contacts
		) {
			const t_collision_mesh& p_mesh = collider.get_collision_mesh();
			const t_collision_mesh& q_mesh = collidee.get_collision_mesh();
			const t_mat44f& p_matrix = collider.get_transform_mat();
			const t_mat44f& q_matrix = collidee.get_transform_mat();

			// ground uses a plane mesh (via set_mesh_scales)
			t_collision_mesh::calc_contacts_wspace(p_mesh, q_mesh, p_matrix, q_matrix, collider_contacts, collidee_contacts);

			return (collider_contacts.size() + collidee_contacts.size());
		}


		t_vec4f t_planar_body::calc_reaction_force(const t_rigid_body& collidee) const {
			#if (ENABLE_REACTION_FORCES == 1)
			const t_vec4f& n = m_wspace_transform.get_y_vector();
			const t_vec4f& f = collidee.get_lin_force();

			return (n * -std::min(0.0f, f.inner_product(n)));
			#else
			return (t_vec4f::zero_vector());
			#endif
		}

		t_vec4f t_planar_body::calc_collision_force(const t_rigid_body& collidee) const {
			#if (ENABLE_COLLISION_FORCES == 1)
			const t_vec4f& n = m_wspace_transform.get_y_vector();
			const t_vec4f& p = collidee.calc_lin_momentum();

			// project pre-momentum onto normal, then reflect it
			// note: also acts as kinetic friction when OGRC < 1
			const t_vec4f& p0 = p.normalize() * -std::min(0.0f, n.inner_product(p));
			const t_vec4f  p1 = p0.reflect(n) * (collidee.get_force_coeffs()).y();
			// const t_vec4f  p0 = n *  std::min(0.0f, n.inner_product(p));
			// const t_vec4f  p1 = n * -std::min(0.0f, n.inner_product(p)) * (collidee.get_force_coeffs()).y();
			const t_vec4f  cf = p1 - p0;

			// first two can fail numerically when OGRC = 1
			assert((p0.sq_len() * 0.99f) <= p.sq_len());
			assert((p1.sq_len() * 0.99f) <= p0.sq_len());
			assert(cf.inner_product(n) >= 0.0f);

			return cf;
			#else
			return (t_vec4f::zero_vector());
			#endif
		}

		t_vec4f t_planar_body::calc_friction_force(const t_rigid_body& collidee) const {
			#if (ENABLE_FRICTION_FORCES == 1)
			const t_vec4f& n  = m_wspace_transform.get_y_vector();
			const t_vec4f& lf = collidee.get_lin_force();
			const t_vec4f& rf = calc_reaction_force(collidee);

			// calculate how much of <lf> is parallel to surface
			// even a force pushing partially "up" should still be
			// opposed!
			const t_vec4f lfp = lf - (n * n.inner_product(lf));

			// maximum static friction (scale) is rf.len() * OBJECT_GROUND_STATIC_FRICTION_COEFF
			// static friction depends on (and increases with) the force being exerted on the object
			const float lfp_scale = lfp.len();
			const float ffp_scale = rf.len() * (collidee.get_force_coeffs()).w();

			if (lfp_scale > 0.0f && lfp_scale <= ffp_scale)
				return -lfp;
			#endif
			return (collidee.get_lin_force() * 0.0f);
		}
	};
};

