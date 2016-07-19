#ifndef PHYSICS_COLLISION_MESH_HDR
#define PHYSICS_COLLISION_MESH_HDR

#include <cassert>
#include <cstddef>
#include <vector>
#include <limits>

#include "../math/point.hpp"
#include "../math/vector.hpp"
#include "../math/matrix.hpp"
#include "../math/idx_triangle.hpp"
#include "bounding_box.hpp"
#include "contact_data.hpp"

namespace newtonics {
	namespace lib_math {
		struct t_ray;
		struct t_plane;
	};

	namespace lib_physics {
		using namespace lib_math;


		enum {
			COL_TEST_TYPE_SAT_INNER = 0,
			COL_TEST_TYPE_SAT_OUTER = 1,
			COL_TEST_TYPE_BOUND_BOX = 2,
		};

		struct t_collision_mesh;
		struct t_coltest_params {
		public:
			t_coltest_params() {
				m_meshes[0] = nullptr;
				m_meshes[1] = nullptr;
				m_matrices[0] = nullptr;
				m_matrices[1] = nullptr;

				m_test_bits.x() = 1;
				m_test_bits.y() = 1;
				m_test_bits.z() = 1;
				m_test_bits.w() = 0;

				m_test_type.x() = -1u; // invalid
				m_test_type.y() = -1u; // invalid
			}

			const t_collision_mesh* m_meshes[2];
			const t_mat44f* m_matrices[2];

			// which collision-tests should be performed
			t_tup4ui m_test_bits;
			// .x := COL_TEST_TYPE_*, .y := obj_indx {0,1}
			t_tup4ui m_test_type;

			t_tup4ui m_sep_inds[2];
			t_vec4f m_sep_vecs[2];
			t_vec4f m_sep_dsts[2];
		};


		struct t_collision_mesh {
		public:
			t_collision_mesh(unsigned int num_verts = 8, unsigned int num_polys = 12) {
				m_vertices.reserve(num_verts);
				m_vnormals.reserve(num_verts);
				m_polygons.reserve(num_polys);
				m_vmasses.reserve(num_verts);
			}

			void clear() {
				m_vertices.clear();
				m_vnormals.clear();
				m_polygons.clear();
				m_vmasses.clear();
			}

			size_t get_num_verts() const { return (m_vertices.size()); }
			size_t get_num_polys() const { return (m_polygons.size()); }

			const std::vector<t_pos4f>& get_verts() const { return m_vertices; }
			const std::vector<t_idx_tri>& get_polys() const { return m_polygons; }
			const std::vector<float>& get_masses() const { return m_vmasses; }

			const t_pos4f& get_vertex(unsigned int idx) const { return m_vertices[idx]; }
			const t_vec4f& get_normal(unsigned int idx) const { return m_vnormals[idx]; }
			const t_idx_tri& get_tri(unsigned int idx) const { return m_polygons[idx]; }
			float get_mass(unsigned int idx) const { return m_vmasses[idx]; }

			const t_bounding_box& get_bounding_box() const { return m_bounding_box; }
			      t_bounding_box& get_bounding_box()       { return m_bounding_box; }


			void add_vertex(const t_pos4f& v) { m_vertices.push_back(v            ); }
			void add_normal(const t_vec4f& n) { m_vnormals.push_back(n.normalize()); }
			void add_poly(const t_idx_tri& t) { m_polygons.push_back(t            ); }
			void add_mass(float w) { m_vmasses.push_back(w); }

			// these should not be used before all geometry has been added
			void set_vertex(unsigned int idx, const t_pos4f& v) { m_vertices[idx] = v; }
			void set_normal(unsigned int idx, const t_vec4f& n) { m_vnormals[idx] = n; }
			void set_poly(unsigned int idx, const t_idx_tri& t) { m_polygons[idx] = t; }
			void set_mass(unsigned int idx, float w) { m_vmasses[idx] = w; }

			void set_vertex_normals();

			// <point> must be in mesh-space!
			bool point_inside_mesh(const t_pos4f& point_mspace) const;

			static bool check_collision(t_coltest_params* tp, unsigned int t_type, unsigned int p_indx, unsigned int q_indx);
			static bool check_collision_inner(
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
			);
			static bool check_collision_outer(
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
			);

			static unsigned int calc_contacts_wspace(
				const t_collision_mesh& collider_mesh,
				const t_collision_mesh& collidee_mesh,
				const t_mat44f& collider_mspace_to_wspace_mat,
				const t_mat44f& collidee_mspace_to_wspace_mat,
				std::vector<t_contact_vertex>& collider_contacts,
				std::vector<t_contact_vertex>& collidee_contacts
			);


			t_pos4f calc_ray_intersection(const t_ray& wspace_ray, const t_mat44f& wspace_to_mspace_mat) const;
			// note: COM can never be outside a convex mesh!
			t_pos4f calc_center_of_mass(float mass) const;

			t_mat44f calc_inertia_tensor() const;

			// this version assumes the COM is at the coordinate origin
			float calc_moi_scalar_com(const t_vec4f& axis) const;
			// calculates the MOI relative to a point (e.g. the COM)
			float calc_moi_scalar(const t_pos4f& point) const;
			// calculates the MOI relative to axis through a point
			float calc_moi_scalar(const t_pos4f& point, const t_vec4f& axis) const;

			float calc_mass() const;
			float calc_bounding_box();

			static void create_plane_mesh(t_collision_mesh& mesh, const t_vec4f& scales);
			static void create_cube_mesh(t_collision_mesh& mesh, const t_vec4f& scales);
			static void create_pyramid_mesh(t_collision_mesh& mesh, const t_vec4f& scales);

		private:
			std::vector<t_pos4f> m_vertices;
			std::vector<t_vec4f> m_vnormals;
			// indexed triangles referencing m_vertices
			std::vector<t_idx_tri> m_polygons;
			// per-vertex masses to determine COM
			std::vector<float> m_vmasses;

			t_bounding_box m_bounding_box;
		};
	};
};

#endif

