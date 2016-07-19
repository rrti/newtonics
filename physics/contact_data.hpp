#ifndef PHYSICS_CONTACT_DATA_HDR
#define PHYSICS_CONTACT_DATA_HDR

#include "../math/point.hpp"
#include "../math/vector.hpp"

namespace newtonics {
	namespace lib_physics {
		using namespace lib_math;

		struct t_contact_vertex {
		public:
			t_contact_vertex(
				const t_pos4f& point = t_pos4f::zero_point(),
				const t_vec4f& normal = t_vec4f::zero_vector()
			) {
				set_point(point);
				set_normal(normal);
			}

			bool operator == (const t_contact_vertex& contact) const {
				return (m_point.equals(contact.get_point()) && m_normal.equals(contact.get_normal()));
			}

			#if 0
			bool operator < (const t_contact_vertex& contact) const {
				const t_pos4f& p = contact.get_point();
				const t_vec4f& n = contact.get_normal();

				#if 0
				// this does not satisfy weak partial ordering constraints
				// (so can not be used as criterion for e.g. set insertion)
				if ((m_point.x() - p.x()) <= -t_tup4f::eps_scalar())
					return true;
				if ((m_point.y() - p.y()) <= -t_tup4f::eps_scalar())
					return true;
				if ((m_point.z() - p.z()) <= -t_tup4f::eps_scalar())
					return true;

				if ((m_normal.x() - n.x()) <= -t_tup4f::eps_scalar())
					return true;
				if ((m_normal.y() - n.y()) <= -t_tup4f::eps_scalar())
					return true;
				if ((m_normal.z() - n.z()) <= -t_tup4f::eps_scalar())
					return true;

				#else

				if (m_point.x() < p.x())
					return true;
				if (m_point.y() < p.y())
					return true;
				if (m_point.z() < p.z())
					return true;

				if (m_normal.x() < n.x())
					return true;
				if (m_normal.y() < n.y())
					return true;
				if (m_normal.z() < n.z())
					return true;
				#endif

				return false;
			}
			#endif

			void set_point(const t_pos4f& point) { m_point = point; }
			void set_normal(const t_vec4f& normal) { m_normal = normal; }

			const t_pos4f& get_point() const { return m_point; }
			const t_vec4f& get_normal() const { return m_normal; }

		private:
			// both in world-space
			t_pos4f m_point;
			t_vec4f m_normal;
		};



		struct t_shared_contact_vertex {
		public:
			t_shared_contact_vertex(
				const t_pos4f& contact_point = t_pos4f::zero_point(),
				const t_vec4f& collider_normal = t_vec4f::zero_vector(),
				const t_vec4f& collidee_normal = t_vec4f::zero_vector()
			) {
				set_contact_point(contact_point);
				set_collider_normal(collider_normal);
				set_collidee_normal(collidee_normal);
			}

			bool operator == (const t_shared_contact_vertex& contact) const {
				if (m_contact_point != contact.get_contact_point())
					return false;
				if (m_collider_normal != contact.get_collider_normal())
					return false;
				if (m_collidee_normal != contact.get_collidee_normal())
					return false;

				return true;
			}

			void set_contact_point(const t_pos4f& point) { m_contact_point = point; }
			void set_collider_normal(const t_vec4f& normal) { m_collider_normal = normal; }
			void set_collidee_normal(const t_vec4f& normal) { m_collidee_normal = normal; }

			const t_pos4f& get_contact_point() const { return m_contact_point; }
			const t_vec4f& get_collider_normal() const { return m_collider_normal; }
			const t_vec4f& get_collidee_normal() const { return m_collidee_normal; }

		private:
			// all in world-space
			t_pos4f m_contact_point;
			t_vec4f m_collider_normal;
			t_vec4f m_collidee_normal;
		};


		struct t_contact_response {
		public:
			t_contact_response(
				const t_vec4f& lin_force = t_vec4f::zero_vector(),
				const t_vec4f& ang_force = t_vec4f::zero_vector(),
				const t_vec4f& ang_lever = t_vec4f::zero_vector()
			) {
				set_lin_force(lin_force);
				set_ang_force(ang_force);
				set_ang_lever(ang_lever);
			}

			void set_lin_force(const t_vec4f& f) { m_lin_force = f; }
			void set_ang_force(const t_vec4f& f) { m_ang_force = f; }
			void set_ang_lever(const t_vec4f& l) { m_ang_lever = l; }

			const t_vec4f& get_lin_force() const { return m_lin_force; }
			const t_vec4f& get_ang_force() const { return m_ang_force; }
			const t_vec4f& get_ang_lever() const { return m_ang_lever; }

			// if lever is used, call this when applying the response
			t_vec4f calc_rot_vector() const { return (m_ang_force.outer_product(m_ang_lever)); }

		private:
			// part of contact-force that will cause linear acceleration
			t_vec4f m_lin_force;

			// part of contact-force that will cause angular acceleration
			//
			// note: if ang_lever is unused this is NOT the angular force
			// (component of F orthogonal to its moment-arm aka lever) but
			// rather the angular force crossed with the arm (M x F) which
			// is the axis around which object wants to rotate as a result
			// of ang_force (scaled by magnitude of angular force)
			t_vec4f m_ang_force;
			t_vec4f m_ang_lever;
		};
	};
};

#endif

