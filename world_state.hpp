#ifndef NEWTONICS_WORLD_STATE_HDR
#define NEWTONICS_WORLD_STATE_HDR

#include <vector>

#include "physics/planar_body.hpp"
#include "physics/rigid_body.hpp"
#include "misc/object_grid.hpp"

using namespace newtonics;
using namespace newtonics::lib_math;
using namespace newtonics::lib_physics;

typedef std::vector<t_contact_vertex> t_contacts_vector;
typedef std::pair<t_contacts_vector, t_contacts_vector> t_contacts_pair;

struct t_arg_parser;
struct t_physics_scene;

struct t_world_state {
public:
	void load_scene(unsigned int scene_index, unsigned int num_objects);
	void reload_scene(size_t scene_offset);

	void create_objects(unsigned int num_objects);
	void update_objects(unsigned int frame_num, float delta_time);

	void create_object(unsigned int object_num);
	void delete_object(unsigned int object_num);
	void create_ground(unsigned int object_num);

	void run_headless(const t_arg_parser* arg_parser);

	const t_rigid_body* trace_ray(const t_ray raw_ray, t_pos4f* hit_pos = nullptr) const;

	bool find_colliding_objects(unsigned int frame_num);
	bool find_colliding_objects_sap(unsigned int frame_num);

	void set_object_pair_mapping(unsigned int frame_num);
	void clear_collision_contacts(unsigned int frame_num);
	void calc_collision_contacts(unsigned int frame_num);
	void handle_clipping_objects(unsigned int frame_num);
	void handle_colliding_objects(unsigned int frame_num);

	void apply_common_forces(unsigned int frame_num);
	void apply_friction_forces(unsigned int frame_num);
	bool apply_collision_forces(unsigned int pass_num);
	bool apply_reaction_forces(unsigned int pass_num);

	static void init_global_scenes();
	static void exec_global_tests();

	const std::vector<t_rigid_body>& get_rigid_objects() const { return m_rigid_objects; }
	      std::vector<t_rigid_body>& get_rigid_objects()       { return m_rigid_objects; }

	const std::vector<t_contacts_pair>& get_object_contacts(unsigned int idx = CURR_ATTRIB_IDX) const { return m_object_contacts[idx]; }
	      std::vector<t_contacts_pair>& get_object_contacts(unsigned int idx = CURR_ATTRIB_IDX)       { return m_object_contacts[idx]; }

	const t_physics_scene* get_scene() const { return m_global_scene; }
	// const t_planar_body& get_ground() const { return m_ground_object; }
	const t_rigid_body& get_ground() const { assert(!m_rigid_objects.empty()); return (m_rigid_objects.back()); }

private:
	std::vector<t_rigid_body> m_rigid_objects;
	std::vector<t_grid_object> m_proxy_objects; // proxy-objects tracked by grid

	std::vector<t_contacts_pair> m_object_contacts[2]; // indexed by OBJ_PID; prev. and curr. ctacts for a pair of objs (p,q)
	std::vector<t_contacts_pair> m_ground_contacts[2]; // indexed by OBJ_ID; prev. and curr. ctacts for an obj p and ground

	std::vector< std::pair<t_vec4f, t_vec4f> > m_object_pair_sep_axes; // indexed by OBJ_PID


	// maps OBJ_ID's to vectors of OBJ_PID's (currently unused)
	std::vector< std::vector<unsigned int> > m_object_pair_mapping;

	// stores the indices and impact-times of all colliding object pairs; a
	// pair-index (OBJ_PID) is simply a_id * m_rigid_objects.size() + b_id
	// [if a=p and b=q]  [if a=q and b=p]
	//   i = q * N + p     j = p * N + q
	//   p = i % N         q = j % N
	//   q = i / N         p = j / N
	std::vector<unsigned int> m_object_pair_indices;

	// stores the frame-number at which a pair was last tested for collision
	// and the frame-number at which its contact vertices were last calculated
	std::vector<unsigned int> m_last_col_test_frames;
	std::vector<unsigned int> m_last_ctc_calc_frames;

	t_physics_scene* m_global_scene;
	t_object_grid m_object_grid;
	// t_planar_body m_ground_object;
};

#endif

