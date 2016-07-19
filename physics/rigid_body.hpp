#ifndef PHYSICS_RIGID_BODY_HDR
#define PHYSICS_RIGID_BODY_HDR

#include <vector>
#include "../math/point.hpp"
#include "../math/vector.hpp"
#include "../math/matrix.hpp"
#include "collision_mesh.hpp"
#include "phys_state.hpp"

namespace newtonics {
	namespace lib_physics {
		using namespace lib_math;

		struct t_contact_vertex;
		struct t_contact_response;
		struct t_planar_body;

		enum {
			FORCE_TYPE_COLL = 0, // collision
			FORCE_TYPE_REAC = 1, // reaction
			FORCE_TYPE_FRIC = 2, // friction
			NUM_FORCE_TYPES = 3,
		};

		enum {
			PREV_ATTRIB_IDX = 0,
			CURR_ATTRIB_IDX = 1,
		};


		struct t_rigid_body {
		public:
			t_rigid_body(unsigned int id): m_id(id) {
				m_mass = 0.0f;
				m_rmoi = 0.0f;

				m_ignore_lin_forces = false;
				m_ignore_ang_forces = false;
				m_ignore_lin_motion = false;
				m_ignore_ang_motion = false;
			}

			unsigned int get_id() const { return m_id; }

			// execute a frame/tick/timestep; quantities are in units per frame
			// (which makes the simulation more stable) so the <dt> is redundant
			// note that <dt> is always a fixed value
			void update(float dt) {
				//! m_physical_state[PREV_ATTRIB_IDX] = m_physical_state[CURR_ATTRIB_IDX];
				//! m_physical_state[CURR_ATTRIB_IDX] = m_physical_state[CURR_ATTRIB_IDX].integrate(dt);

				store_transform();
				store_forces();


				// calculate rotation-axis and MOI
				// split torque into rotation-axis (xyz) and magnitude (w) parts
				// calculate and set moment of inertia with respect to this axis
				set_rmoi(m_collision_mesh.calc_moi_scalar(m_center_of_mass, calc_rot_axis()));

				// integrate velocity; without (faked) friction objects would keep moving forever
				set_lin_velocity(calc_lin_velocity(dt));
				set_ang_velocity(calc_ang_velocity(dt));

				add_position(get_lin_velocity()); // integrate position (pos += lin_vel)
				add_rotation(get_ang_velocity()); // integrate direction (dir += ang_vel)

				// clear force-vectors
				clear_forces();
			}


			const t_collision_mesh& get_collision_mesh() const { return m_collision_mesh; }
			      t_collision_mesh& get_collision_mesh()       { return m_collision_mesh; }
			const t_mat44f& get_transform_mat(unsigned int idx = CURR_ATTRIB_IDX) const { return m_transform_mat[idx]; }
			      t_mat44f& get_transform_mat(unsigned int idx = CURR_ATTRIB_IDX)       { return m_transform_mat[idx]; }

			t_tup4f get_pos_and_radius() const {
				const t_bounding_box& bbox = m_collision_mesh.get_bounding_box();

				const t_pos4f& p = get_position();
				const t_tup4f  t = t_tup4f(p.x(), p.y(), p.z(), bbox.get_radius());
				return t;
			}

			t_pos4f get_position(unsigned int idx = CURR_ATTRIB_IDX) const {
				const t_mat44f& mat = m_transform_mat[idx];
				const t_vec4f& vec = mat.get_t_vector();
				return (t_pos4f(vec.x(), vec.y(), vec.z()));
			}
			t_mat44f get_rotation(unsigned int idx = CURR_ATTRIB_IDX) const {
				const t_vec4f& x = m_transform_mat[idx].get_x_vector();
				const t_vec4f& y = m_transform_mat[idx].get_y_vector();
				const t_vec4f& z = m_transform_mat[idx].get_z_vector();
				return (t_mat44f(x, y, z));
			}
			t_vec4f get_dif_position() const {
				const t_pos4f cur_pos = get_position(CURR_ATTRIB_IDX);
				const t_pos4f prv_pos = get_position(PREV_ATTRIB_IDX);

				return (cur_pos - prv_pos);
			}

			const t_vec4f& get_lin_velocity(unsigned int idx = CURR_ATTRIB_IDX) const { return m_lin_velocity[idx]; }
			const t_vec4f& get_ang_velocity(unsigned int idx = CURR_ATTRIB_IDX) const { return m_ang_velocity[idx]; }

			const t_pos4f& get_center_of_mass() const { return m_center_of_mass; }

			const t_vec4f& get_lin_force(unsigned int idx = CURR_ATTRIB_IDX) const { return m_lin_force[idx]; }
			const t_vec4f& get_ang_force(unsigned int idx = CURR_ATTRIB_IDX) const { return m_ang_force[idx]; }

			const t_tup4f& get_force_coeffs() const { return m_force_coeffs; }
			const t_tup4f& get_misc_coeffs() const { return m_misc_coeffs; }

			// derived quantities
			//   angmom = (rotation_axis * rmoi) * angular_vel
			//   torque = (rotation_axis * rmoi) * angular_acc
			//
			t_vec4f calc_lin_momentum(unsigned int idx = CURR_ATTRIB_IDX) const { return (get_lin_velocity(idx) * m_mass); }
			t_vec4f calc_ang_momentum(unsigned int idx = CURR_ATTRIB_IDX) const { return (get_ang_velocity(idx) * m_rmoi); }
			t_vec4f calc_lin_acceleration(unsigned int idx = CURR_ATTRIB_IDX) const { return (get_lin_force(idx) / m_mass); }
			t_vec4f calc_ang_acceleration(unsigned int idx = CURR_ATTRIB_IDX) const { return (get_ang_force(idx) / m_rmoi); }

			// returns a vector parallel to the axis of rotation, the
			// length of which describes how difficult it is to rotate
			// along <axis> (the MOI; length should be equal to return
			// value of mesh.calc_moi_scalar(com, axis) as well as to
			// this.calc_moi_scalar(axis))
			t_vec4f calc_moi_vector(const t_vec4f& axis) const { return (m_inertia_tensor * axis); }
			t_vec4f calc_rot_axis(unsigned int idx = CURR_ATTRIB_IDX) const { return (m_ang_force[idx].normalize()); }

			t_vec4f calc_lin_velocity(float dt) const { return ((get_lin_velocity() * m_misc_coeffs.x() + calc_lin_acceleration()) * dt); }
			t_vec4f calc_ang_velocity(float dt) const { return ((get_ang_velocity() * m_misc_coeffs.y() + calc_ang_acceleration()) * dt); }


			// for a given point in local-space, calculate its linear
			// velocity (in world-space) from the rotation velocities
			// around our local axes
			t_vec4f calc_point_lin_velocity(const t_pos4f& point_ls) const {
				const t_vec4f x_lin_vel = calc_point_lin_velocity(point_ls, t_vec4f::x_axis_vector()) * m_ang_velocity[CURR_ATTRIB_IDX].x();
				const t_vec4f y_lin_vel = calc_point_lin_velocity(point_ls, t_vec4f::y_axis_vector()) * m_ang_velocity[CURR_ATTRIB_IDX].y();
				const t_vec4f z_lin_vel = calc_point_lin_velocity(point_ls, t_vec4f::z_axis_vector()) * m_ang_velocity[CURR_ATTRIB_IDX].z();
				return (x_lin_vel + y_lin_vel + z_lin_vel);
			}

			t_vec4f calc_point_lin_velocity(const t_pos4f& point_ls, const t_vec4f& axis_ls) const {
				const t_vec4f pvr = point_ls - m_center_of_mass; // COM-vector

				const t_vec4f pvp = pvr - (pvr * axis_ls); // nullify a-component of pvr
				// const t_vec4f pvp = pvr - axis_ls * (pvr.inner_product(axis_ls)); // nullify a-component of pvr

				const t_vec4f lta = pvp.outer_product(axis_ls); // calculate local-tangent axis
				const t_vec4f wta = m_transform_mat[CURR_ATTRIB_IDX] * lta; // calculate world-tangent axis
				return wta;
			}


			t_vec4f calc_air_drag_force(float drag_coeff = 1.0f, float air_dens = 1.225f) const {
				t_vec4f drag_force;
				t_vec4f lin_vel = get_lin_velocity();

				const t_bounding_box& bbox = m_collision_mesh.get_bounding_box();

				const float surf_area = 2.0f * M_PI * lib_math::square(bbox.get_radius());
				const float dyn_press = 0.5f * air_dens * get_sq_lin_speed();

				drag_force = -lin_vel.signs() * drag_coeff * dyn_press * surf_area;
				// clamp a=F/m
				// drag_force = drag_force.clamp(-lin_vel.abs(), lin_vel.abs());

				return (drag_force * m_misc_coeffs.z());
			}


			bool ignore_lin_forces() const { return m_ignore_lin_forces; }
			bool ignore_ang_forces() const { return m_ignore_ang_forces; }
			bool ignore_lin_motion() const { return m_ignore_lin_motion; }
			bool ignore_ang_motion() const { return m_ignore_ang_motion; }


			static bool check_collision(const t_rigid_body&, const t_planar_body&);
			static bool check_collision(const t_rigid_body& p, const t_rigid_body& q, t_coltest_params* tp = nullptr);

			static t_vec4f get_separating_axis(const t_rigid_body& p, const t_rigid_body& q, const t_coltest_params& tp);

			// calculate the time (between frames) at which objects come into
			// contact; assumes they are currently colliding but were not the
			// previous frame
			// returns PAST time (w.r.t present), not FUTURE time (w.r.t past)
			static std::pair<float, float> calc_collision_time_rev(t_rigid_body& p, t_rigid_body& q);

			// push objects apart along their world-space contact normals if
			// they are available (e.g. data from a previous time-step); when
			// called, objects must already be colliding
			static bool handle_clipping_objects(
				t_rigid_body& p,
				t_rigid_body& q,
				t_vec4f& p_sep_axis,
				t_vec4f& q_sep_axis,
				std::vector<t_contact_vertex>& pq_contacts_ws,
				std::vector<t_contact_vertex>& qp_contacts_ws
			);



			// for two contacting bodies, generate the forces and moments
			// both exert on one another to prevent them from intersecting
			static bool apply_reaction_forces(
				t_rigid_body& p,
				t_rigid_body& q,
				const std::vector<t_contact_vertex>& pq_contacts_ws,
				const std::vector<t_contact_vertex>& qp_contacts_ws
			);

			// this does nothing unless at least one object is moving
			//
			// momentum is ALWAYS conserved during all types of collision
			// kinetic energy is ONLY conserved during elastic collisions
			// spectrum: 100% elastic : X% {in}elastic : 100% inelastic
			//
			static bool apply_collision_forces(
				t_rigid_body& p,
				t_rigid_body& q,
				const std::vector<t_contact_vertex>& pq_contacts_ws,
				const std::vector<t_contact_vertex>& qp_contacts_ws
			);

			static void apply_friction_forces(
				t_rigid_body& p,
				t_rigid_body& q,
				const std::vector<t_contact_vertex>& pq_contacts_ws,
				const std::vector<t_contact_vertex>& qp_contacts_ws
			);


			static std::pair<t_vec4f, t_vec4f> calc_collision_forces(
				t_rigid_body& p,
				t_rigid_body& q,
				const t_mat44f& p_inv_mat,
				const t_mat44f& q_inv_mat,
				const t_contact_vertex& p_contact_ws,
				const t_contact_vertex& q_contact_ws
			);

			static std::pair<t_vec4f, t_vec4f> calc_friction_forces(
				t_rigid_body& p,
				t_rigid_body& q,
				const t_contact_vertex& p_contact_ws,
				const t_contact_vertex& q_contact_ws
			);


			// these calculate the linear and angular force (moment) resulting from a force acting through a contact-point
			t_contact_response calc_contact_response_raw(
				const t_vec4f& contact_force_ws,
				const t_pos4f& contact_point_ws,
				unsigned int force_type
			) const;
			t_contact_response calc_contact_response(
				const t_vec4f& contact_force_ws,
				const t_contact_vertex& contact_vertex_ws,
				unsigned int force_type
			) const;


			std::vector<t_vec4f> calc_contact_force_weights(
				const std::vector<t_contact_vertex>& contacts_ws
			) const;

			// find the distribution of a COM-force acting over multiple contacts
			void calc_contact_force_distrib(
				const t_vec4f& contact_force_ws,
				const std::vector<t_contact_vertex>& contacts_ws,
				unsigned int force_type
			);


			// contact forces (can) cause both linear and angular acceleration
			//   angular: component of force orthogonal to contact_vector
			//   linear: component of force parallel to contact_vector
			//
			// any non-parallel contact force can be split into components (parallel
			// and orthogonal to contact normal) which when combined yield identical
			// linear and angular accelerations
			//
			void add_clamped_contact_force(
				const t_vec4f& contact_force_ws,
				const t_contact_vertex& contact_vertex_ws,
				const unsigned int force_type
			);
			void add_contact_response(const t_contact_response& force_response);


			float get_kinetic_energy() const { return (0.5f * m_mass * get_sq_lin_speed()); }
			float get_sq_lin_speed() const { return ((get_lin_velocity()).sq_len()); }
			float get_sq_ang_speed() const { return ((get_ang_velocity()).sq_len()); }
			float get_lin_speed() const { return (std::sqrt(get_sq_lin_speed())); }
			float get_ang_speed() const { return (std::sqrt(get_sq_ang_speed())); }
			float get_mass() const { return m_mass; }
			float get_rmoi() const { return m_rmoi; }

			// same as length(calc_moi_vector(axis)) but without
			// incurring a sqrt call, though not necessarily faster
			//
			// because norm(a) is 1, so is that of transpose(a) * a
			// the 1x4 * 4x1 vector product yields a scalar equal to
			// dot(a, a) calculated in a different order, and the DP
			// of any unit-length vector with itself is 1)
			//
			// if <a> is a row-vector, a * a^T is a 1x1 scalar and a^T * a 4x4 matrix
			// if <a> is a col-vector, a * a^T is a 4x4 matrix and a^T * a 1x1 scalar
			//
			// (a^T=1x4 row-vector) * (M=4x4 matrix) * (a=4x1 col-vector)
			// should be a scalar but the compiler thinks a vec4 is left
			// because of the overloaded "vect::operator * (const vect&)"
			// if we write "return (a.transpose() * M * a);"
			float calc_moi_scalar(const t_vec4f& axis) const { return ((axis * m_inertia_tensor).inner_product(axis)); }


			void set_collision_mesh(const t_collision_mesh& m) { m_collision_mesh = m; }
			void set_inertia_tensor(const t_mat44f& m) { m_inertia_tensor = m; }
			void set_transform_mat(const t_mat44f& m, unsigned int idx = CURR_ATTRIB_IDX) { m_transform_mat[idx] = m; }

			void set_position(const t_pos4f& p, unsigned int idx = CURR_ATTRIB_IDX) { if (!m_ignore_lin_motion) { m_transform_mat[idx].set_t_vector((                    p).to_vector()); } }
			void add_position(const t_vec4f& v, unsigned int idx = CURR_ATTRIB_IDX) { if (!m_ignore_lin_motion) { m_transform_mat[idx].set_t_vector((get_position(idx) + v).to_vector()); } }
			void add_rotation(const t_vec4f& r, unsigned int idx = CURR_ATTRIB_IDX) { if (!m_ignore_ang_motion) { m_transform_mat[idx].rotate_xyz_int_ref(r); } }

			void set_lin_velocity(const t_vec4f& v, unsigned int idx = CURR_ATTRIB_IDX) { v.sanity_assert(); m_lin_velocity[idx] = v; }
			void set_ang_velocity(const t_vec4f& v, unsigned int idx = CURR_ATTRIB_IDX) { v.sanity_assert(); m_ang_velocity[idx] = v; }


			void store_transform() {
				// make a copy of our current transform
				m_transform_mat[PREV_ATTRIB_IDX] = m_transform_mat[CURR_ATTRIB_IDX];
			}
			void store_forces() {
				m_lin_force[PREV_ATTRIB_IDX] = m_lin_force[CURR_ATTRIB_IDX];
				m_ang_force[PREV_ATTRIB_IDX] = m_ang_force[CURR_ATTRIB_IDX];
			}
			void clear_forces(unsigned int idx = CURR_ATTRIB_IDX) {
				set_lin_force(t_vec4f::zero_vector(), idx);
				set_ang_force(t_vec4f::zero_vector(), idx);
			}


			void set_lin_force(const t_vec4f& f, unsigned int idx = CURR_ATTRIB_IDX) { assert(f.sq_len() < 1000000.0f); f.sanity_assert(); m_lin_force[idx]  = (f * (1 - m_ignore_lin_forces)); }
			void set_ang_force(const t_vec4f& f, unsigned int idx = CURR_ATTRIB_IDX) { assert(f.sq_len() < 1000000.0f); f.sanity_assert(); m_ang_force[idx]  = (f * (1 - m_ignore_ang_forces)); }
			void add_lin_force(const t_vec4f& f, unsigned int idx = CURR_ATTRIB_IDX) { assert(f.sq_len() < 1000000.0f); f.sanity_assert(); m_lin_force[idx] += (f * (1 - m_ignore_lin_forces)); }
			void add_ang_force(const t_vec4f& f, unsigned int idx = CURR_ATTRIB_IDX) { assert(f.sq_len() < 1000000.0f); f.sanity_assert(); m_ang_force[idx] += (f * (1 - m_ignore_ang_forces)); }

			void set_center_of_mass(const t_pos4f& c) { m_center_of_mass = c; }

			void set_force_coeffs(const t_tup4f& c) { m_force_coeffs = c; }
			void set_misc_coeffs(const t_tup4f& c) { m_misc_coeffs = c; }

			void set_mass(float m) { m_mass = m; }
			void set_rmoi(float m) { m_rmoi = m; }
			void set_ignore_lin_forces(bool b) { m_ignore_lin_forces = b; }
			void set_ignore_ang_forces(bool b) { m_ignore_ang_forces = b; }
			void set_ignore_lin_motion(bool b) { m_ignore_lin_motion = b; }
			void set_ignore_ang_motion(bool b) { m_ignore_ang_motion = b; }

			void print_state(unsigned int frame_num, const char* hdr_str = "") const;

		protected:
			unsigned int m_id;

			t_collision_mesh m_collision_mesh;
			//! t_phys_state m_physical_state[2];

			t_mat44f m_transform_mat[2]; // wspace
			t_mat44f m_inertia_tensor; // constant

			#if 1
			// in units per frame
			t_vec4f m_lin_velocity[2];
			// in radians per frame; contains the rotation-axis (its
			// magnitude is the rotational *speed* around this axis)
			t_vec4f m_ang_velocity[2];

			// net linear force vector acting on us
			t_vec4f m_lin_force[2];
			// net angular force vector (moment) acting on us
			t_vec4f m_ang_force[2];
			#endif

			// point through which NON-contact forces (e.g. gravity) act
			// contact forces act exclusively through contact points on
			// surface of object
			// stored in OBJECT-space, so must be transformed
			t_pos4f m_center_of_mass;

			t_tup4f m_force_coeffs;
			t_tup4f m_misc_coeffs;


			// in kilograms
			float m_mass;
			// moment of inertia; this is always relative to *current* rot-axis
			float m_rmoi;

			bool m_ignore_lin_forces;
			bool m_ignore_ang_forces;
			bool m_ignore_lin_motion; // ignore {set,add}_position, but not set_transform_mat
			bool m_ignore_ang_motion; // ignore {set,add}_rotation
		};
	};
};

#endif

