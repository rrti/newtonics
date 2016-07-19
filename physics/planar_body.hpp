#ifndef PHYSICS_PLANAR_BODY_HDR
#define PHYSICS_PLANAR_BODY_HDR

#include <vector>

#include "../math/plane.hpp"
#include "../math/vector.hpp"
#include "../math/matrix.hpp"
#include "contact_data.hpp"
#include "collision_mesh.hpp"

namespace newtonics {
	namespace lib_physics {
		using namespace lib_math;

		struct t_rigid_body;
		struct t_planar_body {
		public:
			void set_mesh_scales(const t_vec4f& s) { t_collision_mesh::create_plane_mesh(m_plane_mesh, s); }
			void set_transform_mat(const t_matrix44f& m) { m_wspace_transform = m; }

			const t_collision_mesh& get_collision_mesh() const { return m_plane_mesh; }
			const t_matrix44f& get_transform_mat() const { return m_wspace_transform; }

			float get_point_distance(const t_pos4f& p) const {
				const t_vec4f  v = p.to_vector();
				const t_vec4f& y = m_wspace_transform.get_y_vector();
				const t_vec4f& t = m_wspace_transform.get_t_vector();

				// t.y is the origin offset
				return (v.inner_product(y) - t.y());
			}

			static bool check_collision(const t_planar_body&, const t_rigid_body&);
			static bool handle_collision(
				const t_planar_body&,
				      t_rigid_body&,
				std::vector<t_contact_vertex>&,
				std::vector<t_contact_vertex>&
			);

			static bool handle_clipping_object(
				const t_planar_body&,
				      t_rigid_body&,
				const std::vector<t_contact_vertex>&,
				const std::vector<t_contact_vertex>&
			);
			static bool handle_colliding_object(
				const t_planar_body&,
				      t_rigid_body&,
				const std::vector<t_contact_vertex>&,
				const std::vector<t_contact_vertex>&,
				unsigned int
			);

			static unsigned int calc_contacts_wspace(
				const t_planar_body&,
				const t_rigid_body&,
				std::vector<t_contact_vertex>&,
				std::vector<t_contact_vertex>&
			);


			// these assume that check_collision() was called and returned true
			t_vec4f calc_reaction_force(const t_rigid_body&) const;
			t_vec4f calc_collision_force(const t_rigid_body&) const;
			t_vec4f calc_friction_force(const t_rigid_body&) const;

			t_mat44f calc_wspace_transform(const t_vec4f& n) const {
				t_mat44f m;

				// n.xyz is the normal, n.w is the offset
				const t_vec4f y = t_vec4f::xyz_axis_vector() * n;
				const t_vec4f t = t_vec4f::y_axis_vector() * n.w();

				m.set_t_vector(t);
				m.set_xz_vectors_from_y_vector(y);
				return m;
			}

		private:
			// stores normal in y-column, offset in t-column
			t_mat44f m_wspace_transform;
			// holds simple two-triangle geometry
			t_collision_mesh m_plane_mesh;
		};
	};
};

#endif

