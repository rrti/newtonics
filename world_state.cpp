#include <cassert>
#include <cstdio> // printf
#include <cstdlib> // atoi

#include <algorithm> // random_shuffle
#include <vector>

#include "global_defs.hpp"
#include "world_state.hpp"
#include "physics_scene.hpp"
#include "math/lib.hpp"
#include "misc/arg_parser.hpp"
#include "misc/interval_tree.hpp"
#include "misc/scoped_timer.hpp"
#include "misc/uniform_rng.hpp"

static std::array<t_physics_scene, 8> global_scenes;



void t_world_state::exec_global_tests() {
	#if 0
	if (false) {
		t_mat44f m1;
		t_mat44f m2;
		t_vec4f alh;
		t_vec4f arh;

		// m1.rotate_z_ref(M_PI * 0.333f);
		// m1.print();

		m1.rotate_xyz_int_ref(t_vec4f(0.1f, 0.2f, 0.3f)); alh = m1.get_angles_lh(); // lh
		m1.transpose_ref(); arh = m1.get_angles_rh();

		printf("[%s] arh=<%f,%f,%f> alh=<%f,%f,%f>\n", __func__, arh.x(),arh.y(),arh.z(), alh.x(),alh.y(),alh.z());

		m1.set_unit_values();
		m1.rotate_z_ref(M_PI * 0.333f);
		m2.rotate_axis_ref(t_vec4f::z_axis_vector(), M_PI * 0.333f);

		assert(m1 == m2);
		assert(alh == arh);

		assert(signum( 1.0f) ==  1.0f);
		assert(signum( 0.1f) ==  1.0f);
		assert(signum( 0.0f) ==  0.0f);
		assert(signum(-0.1f) == -1.0f);
		assert(signum(-1.0f) == -1.0f);
	}
	#endif

	#if 0
	if (false) {
		unsigned int j = 0;
		unsigned int k = 0;

		{
			t_uniform_real_rng<float> rng;
			t_scoped_timer scoped_timer("signum_br");

			for (unsigned int n = 0; n < 1000u*1000u*1000u; n++) {
				j += signum_br(rng());
			}
		}

		{
			t_uniform_real_rng<float> rng;
			t_scoped_timer scoped_timer("signum");

			for (unsigned int n = 0; n < 1000u*1000u*1000u; n++) {
				k += signum(rng());
			}
		}

		printf("[%s] j=%u k=%u\n", __func__, j, k);
	}
	#endif



	if (true) {
		// contact-finding; collision-force calculation
		printf("[%s]\n", __func__);

		std::vector<t_rigid_body> objects;
		std::vector<t_contact_vertex> pq_contacts;
		std::vector<t_contact_vertex> qp_contacts;

		const t_physics_scene* scene = &global_scenes[0];

		const t_vec4f obj_spacing = t_vec4f(2.0f, 0.0f, 0.0f);
		// const t_vec4f obj_spacing = scene->object_spacing_vector;
		const t_tup4f obj_coeffs = scene->object_force_coeffs;

		for (float restit_coeff: {1.0f, 0.5f, 0.0f}) {
			objects.clear();
			pq_contacts.clear();
			qp_contacts.clear();

			for (unsigned int n = 0; n < 2; n++) {
				objects.push_back(t_rigid_body(n));

				t_rigid_body& obj = objects.back();
				t_collision_mesh& mesh = obj.get_collision_mesh();
				t_mat44f mat;

				t_collision_mesh::create_cube_mesh(mesh, t_vec4f(1.0f, 1.0f, 1.0f, 100.0f));

				mat.set_t_vector(obj_spacing * n);
				obj.set_transform_mat(mat);
				obj.set_mass(mesh.calc_mass());
				obj.set_center_of_mass(mesh.calc_center_of_mass(obj.get_mass()));
				obj.set_inertia_tensor(mesh.calc_inertia_tensor());

				obj.set_force_coeffs(t_tup4f(restit_coeff, obj_coeffs.y(), obj_coeffs.z(), obj_coeffs.w()));
				obj.set_misc_coeffs(scene->object_misc_coeffs);

				obj.set_lin_velocity(t_vec4f(0.1f * (n == 0), 0.0f, 0.0f, 0.0f));
			}

			t_rigid_body& p = objects[0];
			t_rigid_body& q = objects[1];

			// only true with the default spacing for scenes[0]
			assert(t_rigid_body::check_collision(p, q));

			{
				const t_collision_mesh& p_mesh = p.get_collision_mesh();
				const t_collision_mesh& q_mesh = q.get_collision_mesh();
				const t_mat44f& p_matrix = p.get_transform_mat();
				const t_mat44f& q_matrix = q.get_transform_mat();

				t_collision_mesh::calc_contacts_wspace(p_mesh, q_mesh, p_matrix, q_matrix, pq_contacts, qp_contacts);
			}

			assert(!pq_contacts.empty());
			assert(!qp_contacts.empty());
			printf("\trestit_coeff=%f  #pq_contacts=%lu #qp_contacts=%lu\n", restit_coeff, pq_contacts.size(), qp_contacts.size());

			for (size_t n = 0; n < pq_contacts.size(); n++) {
				const t_pos4f  pt = pq_contacts[n].get_point();
				const t_vec4f nvp = pq_contacts[n].get_normal();
				const t_vec4f nvq = qp_contacts[n].get_normal();

				printf("\t\tpt=<%f,%f,%f> nvp=<%f,%f,%f> nvq=<%f,%f,%f>\n", pt.x(),pt.y(),pt.z(), nvp.x(),nvp.y(),nvp.z(), nvq.x(),nvq.y(),nvq.z());
			}

			for (unsigned int n = 0; n < 2; n++) {
				t_rigid_body& obj = objects[n];

				obj.set_lin_velocity(obj.get_lin_velocity(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
				obj.set_ang_velocity(obj.get_ang_velocity(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);

				obj.store_forces();
				obj.clear_forces(CURR_ATTRIB_IDX);

				obj.print_state(0, "\t\t[ws::rgs][pre_apply_cf]");
			}

			// all correct to within epsilon-tolerance, but not when "garbage" contacts
			// are included (in that case *both* tests behave like restit_coeff=0, even
			// though momentum exchange should be independent of them)
			//
			// at high object velocities objects can become interpenetrating within one
			// frame to such a degree that only parallel faces "survive" --> if we embed
			// each object in a bounding-box we can fix this but only if the actual geometry
			// is also a box (!), otherwise binary-search *or* rely on a previous separating
			// axis if available
			for (unsigned int n = 0; n < MAX_FORCE_EXCHANGE_ITERS && t_rigid_body::apply_collision_forces(p, q, pq_contacts, qp_contacts); n++) {
				p.print_state(0, "\t\t\t[ws::rgs][apply_cf_iter]");
				q.print_state(0, "\t\t\t[ws::rgs][apply_cf_iter]");
			}

			for (unsigned int n = 0; n < 2; n++) {
				t_rigid_body& obj = objects[n];

				obj.set_lin_velocity(obj.get_lin_velocity(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
				obj.set_ang_velocity(obj.get_ang_velocity(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);

				obj.set_lin_force(obj.get_lin_force(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
				obj.set_ang_force(obj.get_ang_force(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
			}

			p.print_state(0, "\t\t[ws::rgs][pst_apply_cf]");
			q.print_state(0, "\t\t[ws::rgs][pst_apply_cf]");
		}


		if (true) {
			// contact-force weight derivation
			printf("\n");
			printf("[%s][pre_weights]\n", __func__);

			// note: normals are not needed
			pq_contacts.clear();
			pq_contacts.push_back(t_contact_vertex(t_pos4f(-1.0f, -1.0f, 0.0f))); // bottom-left
			pq_contacts.push_back(t_contact_vertex(t_pos4f(-0.5f, -1.0f, 0.0f)));
			pq_contacts.push_back(t_contact_vertex(t_pos4f( 0.0f, -1.0f, 0.0f))); // bottom-middle
			pq_contacts.push_back(t_contact_vertex(t_pos4f( 0.5f, -1.0f, 0.0f)));
			pq_contacts.push_back(t_contact_vertex(t_pos4f(+1.0f, -1.0f, 0.0f))); // bottom-right

			const std::vector<t_vec4f>& weights = objects[0].calc_contact_force_weights(pq_contacts);
			// objects[0].calc_contact_force_distrib(contact_force, pq_contacts, FORCE_TYPE_COLL);

			// weights.back() is the abs-sum (normalizer)
			assert(weights.size() == (pq_contacts.size() + 1));

			for (size_t n = 0; n < weights.size(); n++) {
				printf("\tidx=%lu wgt=<%f,%f,%f>\n", n, weights[n].x(), weights[n].y(), weights[n].z());
			}

			printf("[%s][pst_weights]\n", __func__);
		}


		if (true) {
			// angular to linear velocity conversion
			objects[0].set_ang_velocity(t_vec4f(0.0f, 0.0f, 0.5f));

			printf("\n");
			printf("[%s][pre_ang_to_lin]\n", __func__);

			const t_vec4f plv_rgt_mid = objects[0].calc_point_lin_velocity(t_pos4f(1.0f, 0.0f, 0.0f));
			const t_vec4f plv_top_mid = objects[0].calc_point_lin_velocity(t_pos4f(1.0f, 1.0f, 0.0f));
			const t_vec4f plv_top_rgt = objects[0].calc_point_lin_velocity(t_pos4f(1.0f, 1.0f, 1.0f));

			printf("[%s][pst_ang_to_lin]\n", __func__);
			printf("\tplv_rgt_mid=<%f,%f,%f>\n", plv_rgt_mid.x(), plv_rgt_mid.y(), plv_rgt_mid.z());
			printf("\tplv_top_mid=<%f,%f,%f>\n", plv_top_mid.x(), plv_top_mid.y(), plv_top_mid.z());
			printf("\tplv_top_rgt=<%f,%f,%f>\n", plv_top_rgt.x(), plv_top_rgt.y(), plv_top_rgt.z());

			objects[0].set_ang_velocity(t_vec4f(0.0f, 0.0f, 0.0f));
		}


		if (false) {
			const t_collision_mesh& mesh = objects[1].get_collision_mesh();
			const t_mat44f& mat = objects[1].get_transform_mat();
			const t_mat44f& inv = mat.invert_affine();

			const t_ray raw_ray(t_pos4f(1.606302f, 0.892725f, 8.722910f), t_vec4f(0.044854f, -0.101708f, -0.993803f));
			const t_ray inf_ray(raw_ray.pos(), raw_ray.dir());

			printf("\n");
			printf("[%s][pre_trace]\n", __func__);

			// trace_ray(t_ray(t_pos4f(3.1f, 0.0f, 2.1f), t_vec4f(0.8f, 0.0f, 0.2f)));
			const t_pos4f hit_pos = mesh.calc_ray_intersection(inf_ray, inv);

			printf("[%s][pst_trace] hit_pos=<%f,%f,%f,w=%f>\n", __func__, hit_pos.x(),hit_pos.y(),hit_pos.z(),hit_pos.w());
		}


		if (false) {
			// MOI calculation; angular responses
			t_uniform_real_rng<float> rng;

			{
				const t_rigid_body& object = objects[0];
				const t_collision_mesh& mesh = object.get_collision_mesh();

				// const t_vec4f raw_axis = t_vec4f::xyz_axis_vector();
				const t_vec4f raw_axis = t_vec4f((rng() * 2) - 1, (rng() * 2) - 1, (rng() * 2) - 1);
				const t_vec4f rot_axis = raw_axis.normalize();
				const t_vec4f moi_axis = object.calc_moi_vector(rot_axis);

				// FIXME? moi_scalar_m is always 1 less than moi_scalar_o
				const float moi_scalar_o = object.calc_moi_scalar(rot_axis);
				const float moi_scalar_m = mesh.calc_moi_scalar(object.get_center_of_mass(), rot_axis);

				printf("\n");
				printf("[%s][R1] moi_scalar{_o=%f,_m=%f} moi_axis=<%f,%f,%f> (length=%f)\n", __func__, moi_scalar_o, moi_scalar_m, moi_axis.x(),moi_axis.y(),moi_axis.z(), moi_axis.len());
			}

			{
				t_rigid_body& object = objects[0];
				object.clear_forces();
				object.set_lin_velocity(t_vec4f::zero_vector());
				object.set_ang_velocity(t_vec4f::zero_vector());
				object.print_state(0, "\t[ws::rgs][pre_add_ctc_force]");

				const t_pos4f contact_point_ws_a  = t_pos4f(-1.0f,  1.0f,   0.0f,  1.0f);
				const t_pos4f contact_point_ws_b  = t_pos4f( 0.0f,  1.0f,   1.0f,  1.0f);
				const t_vec4f contact_normal_ws_a = t_vec4f(-1.0f,  0.0f,   0.0f,  0.0f);
				const t_vec4f contact_normal_ws_b = t_vec4f( 0.0f,  0.0f,   1.0f,  0.0f);
				const t_vec4f contact_force_ws_a  = t_vec4f(10.0f,  0.0f,   0.0f,  0.0f);
				const t_vec4f contact_force_ws_b  = t_vec4f( 0.0f,  0.0f, -10.0f,  0.0f);

				const t_contact_vertex contact_vertex_ws_a = t_contact_vertex(contact_point_ws_a, contact_normal_ws_a);
				const t_contact_vertex contact_vertex_ws_b = t_contact_vertex(contact_point_ws_b, contact_normal_ws_b);
				const t_contact_response contact_response_a = object.calc_contact_response(contact_force_ws_a, contact_vertex_ws_a, FORCE_TYPE_COLL);
				const t_contact_response contact_response_b = object.calc_contact_response(contact_force_ws_b, contact_vertex_ws_b, FORCE_TYPE_COLL);

				const t_vec4f lin_res_a = contact_response_a.get_lin_force();
				const t_vec4f lin_res_b = contact_response_b.get_lin_force();
				const t_vec4f ang_res_a = contact_response_a.get_ang_force();
				const t_vec4f ang_res_b = contact_response_b.get_ang_force();

				// needs "#if 1" in calc_contact_response_raw
				object.add_contact_response(contact_response_a);
				object.add_contact_response(contact_response_b);
				object.print_state(0, "\t[ws::rgs][pst_add_ctc_f::pre_up]");
				object.update(1.0f);
				object.print_state(0, "\t[ws::rgs][pst_add_ctc_f::pst_up]");

				const t_vec4f rot_axis = object.calc_rot_axis(PREV_ATTRIB_IDX); // only before update()

				const float rot_moi_int = object.get_rmoi();
				const float rot_moi_ext = object.calc_moi_scalar(rot_axis);

				printf("[%s][R2] lin_res_{a,b}={<%f,%f,%f>, <%f,%f,%f>} ang_res_{a,b}={<%f,%f,%f>, <%f,%f,%f>}  rot_axis=<%f,%f,%f,w=%f> rot_moi=%f (%f)\n", __func__,
					lin_res_a.x(), lin_res_a.y(), lin_res_a.z(),
					lin_res_b.x(), lin_res_b.y(), lin_res_b.z(),
					ang_res_a.x(), ang_res_a.y(), ang_res_a.z(),
					ang_res_b.x(), ang_res_b.y(), ang_res_b.z(),
					rot_axis.x(), rot_axis.y(), rot_axis.z(), rot_axis.w(),
					rot_moi_int, rot_moi_ext);
				printf("\n");
			}
		}

		if (false) {
			// vector interpolation
			printf("\n");
			printf("[%s]\n", __func__);

			const t_vec4f vz = t_vec4f::x_axis_vector(); // source
			const t_vec4f vw = t_vec4f::z_axis_vector(); // target

			constexpr unsigned int N = 45; // steps
			constexpr unsigned int K =  3; // exponent

			for (unsigned int n = 0; n <= N; n++) {
				const float raw_alpha = n / (N * 1.0f);
				const float pow_alpha = std::pow(raw_alpha, K * 1.0f);

				// inv-transform source and target to axis-aligned vectors
				// interpolate in local (aligned to canonical axes) space,
				// transform back to world
				// expensive; the matrix could be reused
				// const t_vec4f& vi = t_vec4f::slerp(vz, vw, pow_alpha);
				const t_vec4f& vi = lib_math::vec_slerp(vz, vw, pow_alpha);

				printf("\tn=%u alpha=%f vi=<%f,%f,%f>\n", n, pow_alpha, vi.x(), vi.y(), vi.z());
			}
		}
	}
}

void t_world_state::init_global_scenes() {
	global_scenes[0].object_spacing_vector = t_vec4f(5.0f, 0.0f, 0.0f); // queue (spaced)
	global_scenes[1].object_spacing_vector = t_vec4f(2.0f, 0.0f, 0.0f); // queue

	global_scenes[2].object_spacing_vector = t_vec4f(0.0f, 5.0f, 0.0f); // tower (spaced)
	global_scenes[3].object_spacing_vector = t_vec4f(0.0f, 2.0f, 0.0f); // tower

	global_scenes[4].object_spacing_vector = t_vec4f(0.5f, 5.0f, 0.0f); // staircase (spaced)
	global_scenes[5].object_spacing_vector = t_vec4f(0.5f, 2.0f, 0.0f); // staircase

	global_scenes[6].object_spacing_vector = t_vec4f(0.5f, 5.0f, 0.0f); // twisting staircase (spaced)
	global_scenes[7].object_spacing_vector = t_vec4f(0.5f, 2.0f, 0.0f); // twisting staircase

	// global_scenes[0].ground_normal_vector = t_vec4f(0.2f, 0.8f, 0.0f);
	// global_scenes[3].ground_normal_vector = t_vec4f(0.2f, 0.8f, 0.0f);

	for (size_t n = 0; n < global_scenes.size(); n++)
		global_scenes[n].ground_normal_vector.normalize_ref();

	// kill gravity and lower the ground s.t. objects do not come in contact with it
	global_scenes[0].ground_offset_vector = t_vec4f(0.0f, -20.0f, 0.0f);
	global_scenes[0].object_gravity_accell *= 0.0f;
	#if 0
	// in world-space +x is right, +y is up, +z is down
	// global_scenes[0].object_tst_lin_accell = t_vec4f::xyz_axis_vector() * TIME_STEP_SIZE_S * TIME_STEP_SIZE_S * 0.0f;
	global_scenes[0].object_spacing_vector = t_vec4f(6.0f, 0.0f, 2.0f) * 1.0f;
	global_scenes[0].object_mesh_scales = t_vec4f(4.0f, 1.0f, 1.0f, 100.0f);
	global_scenes[0].object_ignore_lin_forces =  true;
	global_scenes[0].object_ignore_ang_forces = false;
	#endif

	global_scenes[3].object_tst_lin_accell *= 0.0f;
	global_scenes[3].object_tst_ang_accell *= 0.0f;
	global_scenes[5].object_tst_lin_accell *= 0.0f;
	global_scenes[5].object_tst_ang_accell *= 0.0f;
	global_scenes[7].object_tst_lin_accell *= 0.0f;
	global_scenes[7].object_tst_ang_accell *= 0.0f;

	global_scenes[6].object_axis_angles.y() = 45.0f * (M_PI / 180.0f);
	global_scenes[7].object_axis_angles.y() = 45.0f * (M_PI / 180.0f);
}



void t_world_state::load_scene(unsigned int scene_index, unsigned int num_objects) {
	const unsigned int num_scenes = (global_scenes.size());
	const unsigned int def_argval = (num_objects == 0);

	scene_index = std::min(scene_index, num_scenes);
	num_objects = std::min(lerp(num_objects, global_scenes[scene_index].num_objects, def_argval), MAX_OBJECT_COUNT);

	m_global_scene = &global_scenes[scene_index];

	#if (ENABLE_PHYSSTATE_DEBUG_PRINTS == 1)
	printf("[%s] scene_gravity=<%f,%f,%f>\n", __func__, m_global_scene->object_gravity_accell.x(), m_global_scene->object_gravity_accell.y(), m_global_scene->object_gravity_accell.z());
	#endif

	// NOTE:
	//   if an object (e.g. in a vertical stack) is too far outside grid-bounds, it
	//   will not be inserted into any cells and skipped during collision-pair tests
	m_object_grid = t_object_grid(m_global_scene->object_grid_divs, m_global_scene->object_grid_mins, m_global_scene->object_grid_maxs);

	create_objects(num_objects);
	create_ground(num_objects);
}

void t_world_state::reload_scene(size_t scene_offset) {
	// reconstruct grid from scratch
	m_object_grid = t_object_grid(m_global_scene->object_grid_divs, m_global_scene->object_grid_mins, m_global_scene->object_grid_maxs);

	const size_t cur_scene_idx = m_global_scene - &global_scenes[0];
	const size_t nxt_scene_idx = (cur_scene_idx + scene_offset) % (global_scenes.size());

	m_global_scene = &global_scenes[nxt_scene_idx];

	// use previous number of created objects to repopulate scene
	create_objects(m_rigid_objects.size() - 1);
	create_ground(m_rigid_objects.size());
}

void t_world_state::create_objects(unsigned int num_objects) {
	assert(m_global_scene != nullptr);

	m_rigid_objects.clear();
	m_proxy_objects.clear();
	m_rigid_objects.reserve(MAX_OBJECT_COUNT);
	m_proxy_objects.reserve(MAX_OBJECT_COUNT);

	// NOTE: if num_objects gets larger, these should probably be hashmaps
	m_object_contacts[CURR_ATTRIB_IDX].clear();
	m_object_contacts[PREV_ATTRIB_IDX].clear();
	m_object_contacts[CURR_ATTRIB_IDX].resize(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT);
	m_object_contacts[PREV_ATTRIB_IDX].resize(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT);

	m_object_pair_sep_axes.clear();
	m_object_pair_sep_axes.resize(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT);

	m_ground_contacts[CURR_ATTRIB_IDX].clear();
	m_ground_contacts[PREV_ATTRIB_IDX].clear();
	m_ground_contacts[CURR_ATTRIB_IDX].resize(MAX_OBJECT_COUNT);
	m_ground_contacts[PREV_ATTRIB_IDX].resize(MAX_OBJECT_COUNT);

	m_object_pair_mapping.clear();
	m_object_pair_indices.clear();
	m_object_pair_mapping.resize(MAX_OBJECT_COUNT);
	m_object_pair_indices.reserve(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT);

	for (unsigned int n = 0; n < MAX_OBJECT_COUNT; n++) {
		m_object_pair_mapping[n].clear();
		m_object_pair_mapping[n].reserve(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT);
	}

	m_last_col_test_frames.clear();
	m_last_ctc_calc_frames.clear();
	m_last_col_test_frames.resize(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT, 0);
	m_last_ctc_calc_frames.resize(MAX_OBJECT_COUNT * MAX_OBJECT_COUNT, 0);

	for (unsigned int n = 0; n < num_objects; n++) {
		create_object(n);
	}

	#if (RANDOMIZE_OBJECT_ORDER == 1)
	// randomize objects to test independence on ordering, excluding ground
	std::random_shuffle(m_rigid_objects.begin(), m_rigid_objects.end());
	#endif

	for (size_t n = 0; n < m_proxy_objects.size(); n++) {
		t_rigid_body& rb = m_rigid_objects[n];
		t_grid_object& go = m_proxy_objects[n];

		go.set_unique_id(rb.get_id());

		if (m_object_grid.add_object(&go, rb.get_pos_and_radius()) == 0) {
			assert(false);
		}
	}
}

void t_world_state::create_object(unsigned int object_num) {
	m_rigid_objects.push_back(t_rigid_body(object_num));
	m_proxy_objects.push_back(t_grid_object(object_num));

	t_rigid_body& object = m_rigid_objects.back();
	t_collision_mesh& mesh = object.get_collision_mesh();
	t_mat44f matrix;

	t_collision_mesh::create_cube_mesh(mesh, m_global_scene->object_mesh_scales);

	matrix.rotate_xyz_int_ref(m_global_scene->object_axis_angles * object_num);
	matrix.set_t_vector(m_global_scene->object_spacing_vector * object_num + m_global_scene->object_offset_vector);

	object.set_transform_mat(matrix);
	object.store_transform();
	object.set_mass(mesh.calc_mass());
	object.set_center_of_mass(mesh.calc_center_of_mass(object.get_mass()));
	object.set_inertia_tensor(mesh.calc_inertia_tensor());

	object.set_force_coeffs(m_global_scene->object_force_coeffs);
	object.set_misc_coeffs(m_global_scene->object_misc_coeffs);

	object.set_ignore_lin_forces(m_global_scene->object_ignore_lin_forces);
	object.set_ignore_ang_forces(m_global_scene->object_ignore_ang_forces);
}

void t_world_state::create_ground(unsigned int object_num) {
	#if 0
	// use a plane to represent the ground: interpenetration-resolution
	// between meshes interferes with regular object (position) updates
	m_ground_object.set_mesh_scales(m_global_scene->ground_mesh_scales);
	m_ground_object.set_transform_mat(m_ground_object.calc_wspace_transform(m_global_scene->ground_normal_vector + m_global_scene->ground_offset_vector));
	#endif

	m_rigid_objects.push_back(t_rigid_body(object_num));
	m_proxy_objects.push_back(t_grid_object(object_num));

	t_rigid_body& object = m_rigid_objects.back();
	t_collision_mesh& mesh = object.get_collision_mesh();
	t_mat44f matrix;

	t_collision_mesh::create_cube_mesh(mesh, m_global_scene->ground_mesh_scales);

	matrix.set_t_vector(m_global_scene->ground_offset_vector);
	matrix.rotate_xyz_int_ref(m_global_scene->ground_axis_angles);

	object.set_transform_mat(matrix);
	// do not set infinite mass; ground is not excluded from force calculations
	object.set_mass(m_global_scene->ground_mesh_scales.w());

	// ground ignores all forces and motion
	object.set_ignore_lin_forces(true);
	object.set_ignore_ang_forces(true);
	object.set_ignore_lin_motion(true);
	object.set_ignore_ang_motion(true);

	// NOTE: will insert into far too many cells, needs a (mins,maxs) variant
	m_object_grid.add_object(&m_proxy_objects.back(), object.get_pos_and_radius());
}

void t_world_state::update_objects(unsigned int frame_num, float delta_time) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	// forces that are not the result of contacts
	apply_common_forces(frame_num);
	// clear any expired contact data (not just for current pairs)
	clear_collision_contacts(frame_num);

	if (find_colliding_objects(frame_num)) {
		set_object_pair_mapping(frame_num);
		calc_collision_contacts(frame_num);
		handle_clipping_objects(frame_num);
		// all object-to-object contact forces
		handle_colliding_objects(frame_num);
	}

	for (size_t k = 0; k < (m_rigid_objects.size() - 1); k++) {
		#if 0
		// needs to be done after propagating reaction-forces among
		// colliding pairs, since the ground is not a "real" object
		// and will not be part of such pairs
		if (t_planar_body::check_collision(m_ground_object, m_rigid_objects[k])) {
			#if (ENABLE_PLANEBODY_DEBUG_PRINTS == 1)
			m_rigid_objects[k].print_state(frame_num, "\t[pb::handle_col_ground][pre]");
			#endif

			std::vector<t_contact_vertex>& collider_contacts = m_ground_contacts[CURR_ATTRIB_IDX][k].first;
			std::vector<t_contact_vertex>& collidee_contacts = m_ground_contacts[CURR_ATTRIB_IDX][k].second;

			t_planar_body::handle_collision(m_ground_object, m_rigid_objects[k], collider_contacts, collidee_contacts);

			// must hold until objects fall off the ground-plane
			// assert((m_rigid_objects[k].get_position()).y() >= -t_tup4f::eps_scalar());

			#if (ENABLE_PLANEBODY_DEBUG_PRINTS == 1)
			m_rigid_objects[k].print_state(frame_num, "\t[pb::handle_col_ground][pst]");
			#endif
		}
		#endif

		m_object_grid.del_object(&m_proxy_objects[k]);
		m_rigid_objects[k].print_state(frame_num, "\t[ws::update_objs][pre]");
		m_rigid_objects[k].update(delta_time);
		m_rigid_objects[k].print_state(frame_num, "\t[ws::update_objs][pst]");
		m_object_grid.add_object(&m_proxy_objects[k], m_rigid_objects[k].get_pos_and_radius());
	}
}


bool t_world_state::find_colliding_objects(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	// first find all pairs of colliding objects and their
	// (world-space) points of contact, apply forces later
	//
	// if N=3 then there are (N*(N-1))/2 = (3*2)/2 = 3 unique
	// pairs; index of a pair (p,q) can be calculated directly
	// by f(p,q) = q * num_objects + p
	//
	//  <0,0>*  <0,1>   <0,2>   <0,3>
	//  <1,0>*  <1,1>*  <1,2>   <1,3>
	//  <2,0>*  <2,1>*  <2,2>*  <2,3>
	//  <3,0>*  <3,1>*  <3,2>*  <3,3>*
	//
	// pairs marked* are skipped cq. not considered for testing
	// note that collision detection should be fully symmetric:
	// if A collides with B, then B collides with A
	//
	// forces need to be applied to BOTH objects but we do not
	// want to calculate them twice for each pair, so q starts
	// at p + 1 (not at 0) --> filters out the non-unique ones
	//
	// all counters start at 0, but so does frame
	frame_num += 1;

	for (size_t k = 0; k < m_rigid_objects.size(); k++) {
		const t_rigid_body& p_obj = m_rigid_objects[k];
		const unsigned int p_id = p_obj.get_id();

		const t_vec4f& p_vel = p_obj.get_lin_velocity();
		const t_tup4f&  p_pos = p_obj.get_pos_and_radius();


		// get all objects in the vicinity of p (taking velocity
		// into account; note that the grid touches every object
		// only once if it exists in multiple cells)
		t_object_grid::t_object_vector object_vector;
		m_object_grid.get_objects_in_cube_vel(t_pos4f(p_pos.x(), p_pos.y(), p_pos.z()), p_vel.abs(), object_vector);

		for (const t_grid_object* proxy_object: object_vector) {
			const unsigned int  q_id = proxy_object->get_unique_id();
			const unsigned int pq_id = p_id * m_rigid_objects.size() + q_id;
			const unsigned int qp_id = q_id * m_rigid_objects.size() + p_id;

			const t_rigid_body& q_obj = m_rigid_objects[q_id];

			if (p_id == q_id)
				continue;

			assert(q_id == q_obj.get_id());
			assert(pq_id != qp_id);

			// check if object-pair (p,q) was already tested as (q,p) this frame
			if (m_last_col_test_frames[qp_id] == frame_num)
				continue;
			if (m_last_col_test_frames[pq_id] == frame_num)
				continue;

			// store the frame-number for test(p,q) and its mirror test(q,p)
			m_last_col_test_frames[qp_id] = frame_num;
			m_last_col_test_frames[pq_id] = frame_num;

			// if we do not already have a sep-axis, calculate one ASAP
			// (turns off the early-out test until a valid one is found)
			t_coltest_params tp;
			tp.m_test_bits.x() = 1 - (m_object_pair_sep_axes[pq_id].first == t_vec4f::zero_vector() || m_object_pair_sep_axes[pq_id].second == t_vec4f::zero_vector());
			tp.m_test_bits.w() = 0;

			if (!t_rigid_body::check_collision(p_obj, q_obj, &tp)) {
				m_object_pair_sep_axes[pq_id].first  =  t_rigid_body::get_separating_axis(p_obj, q_obj, tp);
				m_object_pair_sep_axes[pq_id].second = -t_rigid_body::get_separating_axis(p_obj, q_obj, tp);

				// if the early-out is disabled (x()=0) and the objects are
				// found *not* to be colliding, we should *always* get back
				// a valid axis
				assert(tp.m_test_bits.x() == 1 || m_object_pair_sep_axes[pq_id].first  != t_vec4f::zero_vector());
				assert(tp.m_test_bits.x() == 1 || m_object_pair_sep_axes[pq_id].second != t_vec4f::zero_vector());
				continue;
			}

			// store the (p,q) object-pair for later handling (p*N+q)
			m_object_pair_indices.push_back(pq_id);
		}

		object_vector.clear();
	}

	#if (ENABLE_PHYSSTATE_DEBUG_PRINTS == 1)
	printf("[ws::%s][f=%u] #col_pairs=%lu\n", __func__, frame_num - 1, m_object_pair_indices.size());
	#endif

	return (!m_object_pair_indices.empty());
}

bool t_world_state::find_colliding_objects_sap(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	frame_num += 1;

	#if 0
	const t_vec4f axes[] = {x_axis_vector(), y_axis_vector(), z_axis_vector()};
	#else
	// defaults to identity (wx,wy,wz)
	const t_mat44f axes;
	#endif

	struct t_object_proj {
	public:
		t_object_proj(unsigned int id = -1u, float pos = 0.0f): m_obj_id(id), m_obj_pos(pos) {}
		bool operator < (const t_object_proj& p) const { return (m_obj_pos < p.pos()); }

		unsigned int id() const { return m_obj_id; }
		float pos() const { return m_obj_pos; }

	private:
		unsigned int m_obj_id;
		// projected position (min or max)
		float m_obj_pos;
	};

	struct t_object_pair {
		unsigned int p_id;
		unsigned int q_id;
	};


	std::vector<t_object_proj> obj_projections[3];
	std::vector<t_object_pair> col_candidates[3];

	obj_projections[0].reserve(m_rigid_objects.size() * 2);
	obj_projections[1].reserve(m_rigid_objects.size() * 2);
	obj_projections[2].reserve(m_rigid_objects.size() * 2);
	col_candidates[0].reserve(m_rigid_objects.size());
	col_candidates[1].reserve(m_rigid_objects.size());
	col_candidates[2].reserve(m_rigid_objects.size());

	// project objects onto coordinate axes
	// while sweeping, whenever we encounter a min-coordinate of p
	// *not* directly followed by the max-coordinate of p there is
	// a potential colliding pair --> still has O(N^2) pathological
	// performance
	for (size_t n = 0; n < m_rigid_objects.size(); n++) {
		const unsigned int p_id = m_rigid_objects[n].get_id();

		const t_rigid_body& p_obj = m_rigid_objects[p_id];
		const t_collision_mesh& p_mesh = p_obj.get_collision_mesh();
		const t_bounding_box& p_bbox = p_mesh.get_bounding_box();
		const t_pos4f& p_pos = p_obj.get_position();

		for (unsigned int d = 0; d < 3; d++) {
			// project position onto current axis-vector
			const t_vec4f pro_vec = axes.get_col_vector(d); // axes[d];
			const t_pos4f pro_pos = (p_pos.to_vector() * pro_vec).to_point();

			// find the (bounding-sphere) extrema along this axis
			const t_pos4f& max_pos = pro_pos + (pro_vec * p_bbox.get_radius());
			const t_pos4f& min_pos = pro_pos - (pro_vec * p_bbox.get_radius());

			obj_projections[d].push_back({p_id, min_pos[d]});
			obj_projections[d].push_back({p_id, max_pos[d]});
		}
	}

	// sort projections by projected position
	std::sort(obj_projections[0].begin(), obj_projections[0].end());
	std::sort(obj_projections[1].begin(), obj_projections[1].end());
	std::sort(obj_projections[2].begin(), obj_projections[2].end());

	// add sentinels
	obj_projections[0].push_back({-1u, M_FINF});
	obj_projections[1].push_back({-1u, M_FINF});
	obj_projections[2].push_back({-1u, M_FINF});


	for (unsigned int d = 0; d < 3; d++) {
		// last pair is (p, q=sentinel)
		const size_t num_projections_ex = obj_projections[d].size() - 1;

		for (size_t n = 0; n < num_projections_ex; n++) {
			const t_object_proj& p = obj_projections[d][n];

			// without the inner loop this would only consider *immediately*
			// adjacent objects and consequently there are edge-cases where
			// overlaps will be missed; if the RHS (q) of this pair entirely
			// spans the LHS (p), then it can also intersect p's predecessor(s)
			// even if (p,q) itself has no overlap
			//
			// using more sets of projection axes presents only half a solution,
			// while sweeping left/right *always* works and is unlikely to cause
			// O(N^2) performance
			//
			// alternative: add objects into bins based on projected position and
			// size, then for each bin for each object do overlap tests (not good,
			// this will degrade if all objects project into the same bin)
			//
			for (size_t k = n + 1; k < num_projections_ex; k++) {
				const t_object_proj& q = obj_projections[d][k];

				// if the min and max projection match, then by definition the
				// following pair ([k]=max, [k+1]=nxt) will not and we need to
				// skip it
				if (p.id() == q.id()) {
					n++; break;
				}

				col_candidates[d].push_back({p.id(), q.id()});
			}
		}

		// if no overlaps exist along this dimension, bail early
		if (col_candidates[d].empty())
			return false;
	}

	#if 0
	// TODO:
	//   inner loop must be over the dimensions, s.t. a pair can be pre-checked
	//   in all three and rejected *before* the actual collision test takes place
	//
	// each dimension contains an equal number of projections, just pick [0]
	const size_t num_projections_ex = obj_projections[0].size() - 1;

	for (size_t n = 0; n < num_projections_ex; n++) {
		for (size_t k = n + 1; k < num_projections_ex; k++) {
			unsigned int num_dim_overlaps = 0;

			t_object_proj p;
			t_object_proj q;

			for (unsigned int d = 0; d < 3; d++) {
				p = obj_projections[d][n];
				q = obj_projections[d][k];

				if (p.id() == q.id())
					break;

				num_dim_overlaps += 1;
			}

			if (num_dim_overlaps == 3) {
				// do not need to store the candidates per dimension
				col_candidates[0].push_back({p.id(), q.id()});
			} else {
				break;
			}
		}
	}
	#endif

	// for each unique pair of candidates, test for actual collision
	for (unsigned int d = 0; d < 3; d++) {
		for (size_t n = 0; n < col_candidates[d].size(); n++) {
			const t_object_pair& pair = col_candidates[d][n];

			const unsigned int p_id = pair.p_id;
			const unsigned int q_id = pair.q_id;

			const unsigned int pq_id = p_id * m_rigid_objects.size() + q_id;
			const unsigned int qp_id = q_id * m_rigid_objects.size() + p_id;

			assert(p_id != q_id);
			assert(pq_id != qp_id);

			const t_rigid_body& p_obj = m_rigid_objects[p_id];
			const t_rigid_body& q_obj = m_rigid_objects[q_id];

			if (m_last_col_test_frames[qp_id] == frame_num)
				continue;
			if (m_last_col_test_frames[pq_id] == frame_num)
				continue;

			m_last_col_test_frames[qp_id] = frame_num;
			m_last_col_test_frames[pq_id] = frame_num;

			t_coltest_params tp;
			tp.m_test_bits.x() = 1 - (m_object_pair_sep_axes[pq_id].first == t_vec4f::zero_vector() || m_object_pair_sep_axes[pq_id].second == t_vec4f::zero_vector());
			tp.m_test_bits.w() = 0;

			// note: here we could also simply store axes[d]
			if (!t_rigid_body::check_collision(p_obj, q_obj, &tp)) {
				m_object_pair_sep_axes[pq_id].first  =  t_rigid_body::get_separating_axis(p_obj, q_obj, tp);
				m_object_pair_sep_axes[pq_id].second = -t_rigid_body::get_separating_axis(p_obj, q_obj, tp);

				assert(tp.m_test_bits.x() == 1 || m_object_pair_sep_axes[pq_id].first  != t_vec4f::zero_vector());
				assert(tp.m_test_bits.x() == 1 || m_object_pair_sep_axes[pq_id].second != t_vec4f::zero_vector());
				continue;
			}

			m_object_pair_indices.push_back(pq_id);
		}
	}

	#if (ENABLE_PHYSSTATE_DEBUG_PRINTS == 1)
	printf("[ws::%s][f=%u] #col_pairs=%lu\n", __func__, frame_num - 1, m_object_pair_indices.size());
	#endif

	return (!m_object_pair_indices.empty());
}


void t_world_state::set_object_pair_mapping(unsigned int frame_num) {
	(void) frame_num;

	// construct a mapping from object ID's to collision pair ID's
	for (size_t n = 0; n < m_object_pair_indices.size(); n++) {
		const unsigned int pair = m_object_pair_indices[n];
		const unsigned int p_id = pair / m_rigid_objects.size();
		const unsigned int q_id = pair % m_rigid_objects.size();

		assert(p_id < MAX_OBJECT_COUNT);
		assert(q_id < MAX_OBJECT_COUNT);

		m_object_pair_mapping[p_id].push_back(pair);
		m_object_pair_mapping[q_id].push_back(pair);
	}
}


void t_world_state::handle_clipping_objects(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	assert(!m_object_pair_indices.empty());

	// use either contact normals or some previous separating axis
	// (or even binary-search rollback) to resolve clipping objects
	// by iteratively moving them apart
	(void) frame_num;

	unsigned int pair_cnt = 1;
	unsigned int pass_num = 0;

	while (pair_cnt != 0) {
		pair_cnt = 0;

		for (size_t n = 0; n < m_object_pair_indices.size(); n++) {
			const unsigned int pair = m_object_pair_indices[n];
			const unsigned int p_id = pair / m_rigid_objects.size();
			const unsigned int q_id = pair % m_rigid_objects.size();

			const unsigned int pq_id = p_id * m_rigid_objects.size() + q_id;
			// const unsigned int qp_id = q_id * m_rigid_objects.size() + p_id;

			assert(pq_id == pair);

			t_vec4f& p_axis = m_object_pair_sep_axes[pq_id].first;
			t_vec4f& q_axis = m_object_pair_sep_axes[pq_id].second;

			t_rigid_body& p = m_rigid_objects[p_id];
			t_rigid_body& q = m_rigid_objects[q_id];

			// any data from previous frame(s) will be cleared automatically
			// by reference; this allows using any contacts from the previous
			// time-step as seed for resolving interpenetrations
			//
			// contacts[pair].1st := ws-points of q_mesh in contact with p_mesh [1]
			// contacts[pair].2nd := ws-points of p_mesh in contact with q_mesh [2]
			//
			t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
			t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

			pair_cnt += t_rigid_body::handle_clipping_objects(p, q, p_axis, q_axis, pq_contacts, qp_contacts);
		}

		pass_num += 1;
		pair_cnt *= (pass_num < lib_math::square(m_rigid_objects.size()));
	}
}


void t_world_state::clear_collision_contacts(unsigned int frame_num) {
	#if 0
	for (size_t n = 0; n < m_object_pair_indices.size(); n++) {
		const unsigned int pair = m_object_pair_indices[n];
	#endif
	for (size_t pair = 0; pair < m_object_contacts[CURR_ATTRIB_IDX].size(); pair++) {
		const unsigned int p_id = pair / m_rigid_objects.size();
		const unsigned int q_id = pair % m_rigid_objects.size();

		// m_object_contacts has room for MAX_OBJECT_COUNT^2 total pairs
		// active pair-ID's are composed of object-ID's within the range
		// [0, m_objects.size() - 1], so abort when p_id becomes illegal
		if (p_id >= m_rigid_objects.size())
			break;
		if (q_id >= m_rigid_objects.size())
			continue;

		// get the age (in frames) of (p,q)'s last contact-data
		const unsigned int pq_contact_age = frame_num - m_last_ctc_calc_frames[q_id * m_rigid_objects.size() + p_id];
		const unsigned int qp_contact_age = frame_num - m_last_ctc_calc_frames[p_id * m_rigid_objects.size() + q_id];

		assert(pq_contact_age == qp_contact_age);

		t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
		t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

		if (MAX_OBJECT_CONTACT_AGE == 0 || pq_contact_age >= MAX_OBJECT_CONTACT_AGE) {
			pq_contacts.clear();
			qp_contacts.clear();
		}
	}
}

void t_world_state::calc_collision_contacts(unsigned int frame_num) {
	for (size_t n = 0; n < m_object_pair_indices.size(); n++) {
		const unsigned int pair = m_object_pair_indices[n];
		const unsigned int p_id = pair / m_rigid_objects.size();
		const unsigned int q_id = pair % m_rigid_objects.size();

		const unsigned int pq_id = p_id * m_rigid_objects.size() + q_id;
		const unsigned int qp_id = q_id * m_rigid_objects.size() + p_id;

		assert(pq_id == pair);

		// NOTE:
		//   handle_clipping_objects which runs after us will just clear
		//   the contacts if they are non-empty, i.e. when no separation
		//   axis currently exists --> no point calculating them if one
		//   already does
		if (m_object_pair_sep_axes[pq_id].first != t_vec4f::zero_vector())
			continue;

		t_rigid_body& p = m_rigid_objects[p_id];
		t_rigid_body& q = m_rigid_objects[q_id];

		t_contacts_vector& pq_prv_contacts = m_object_contacts[PREV_ATTRIB_IDX][pair].first;
		t_contacts_vector& qp_prv_contacts = m_object_contacts[PREV_ATTRIB_IDX][pair].second;
		t_contacts_vector& pq_cur_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
		t_contacts_vector& qp_cur_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

		// temporarily move current data out of the way
		pq_prv_contacts.swap(pq_cur_contacts);
		qp_prv_contacts.swap(qp_cur_contacts);

		// calculate new contact-data among colliding pairs
		{
			const t_collision_mesh& p_mesh = p.get_collision_mesh();
			const t_collision_mesh& q_mesh = q.get_collision_mesh();
			const t_mat44f& p_matrix = p.get_transform_mat();
			const t_mat44f& q_matrix = q.get_transform_mat();

			t_collision_mesh::calc_contacts_wspace(p_mesh, q_mesh, p_matrix, q_matrix, pq_cur_contacts, qp_cur_contacts);
		}

		if (pq_cur_contacts.empty()) {
			// no new contacts, switch back to current data (if any)
			pq_prv_contacts.swap(pq_cur_contacts);
			qp_prv_contacts.swap(qp_cur_contacts);
		} else {
			// reset axes if we have new contacts (which hopefully
			// will resolve the collision s.t. we can find new axes)
			m_object_pair_sep_axes[pq_id].first  = t_vec4f::zero_vector();
			m_object_pair_sep_axes[pq_id].second = t_vec4f::zero_vector();

			// new contact-data, update timestamps
			m_last_ctc_calc_frames[qp_id] = frame_num;
			m_last_ctc_calc_frames[pq_id] = frame_num;

			pq_prv_contacts.clear();
			qp_prv_contacts.clear();
		}

		assert(pq_prv_contacts.empty());
		assert(qp_prv_contacts.empty());
	}


	#if (ENABLE_PHYSSTATE_DEBUG_PRINTS == 1)
	printf("[ws::%s][f=%u] #col_pairs=%lu\n", __func__, frame_num, m_object_pair_indices.size());
	#endif

	#if (ENABLE_PHYSSTATE_DEBUG_PRINTS == 123)
	for (size_t n = 0; n < m_object_pair_indices.size(); n++) {
		const unsigned int pair = m_object_pair_indices[n];
		const unsigned int p_id = pair / m_rigid_objects.size();
		const unsigned int q_id = pair % m_rigid_objects.size();

		t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
		t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

		printf("\tpair=<p=%u,q=%u> #pq_contacts=%lu #qp_contacts=%lu\n", p_id, q_id, pq_contacts.size(), qp_contacts.size());

		for (size_t k = 0; k < pq_contacts.size(); k++) {
			const t_contact_vertex& cv = pq_contacts[k];
			const t_pos4f& pnt = cv.get_point();
			const t_vec4f& nrm = cv.get_normal();

			printf("\t\t[p2q] k=%lu cv=<P=<%f,%f,%f> N=<%f,%f,%f>>\n", k, pnt.x(),pnt.y(),pnt.z(), nrm.x(),nrm.y(),nrm.z());
		}
		for (size_t k = 0; k < qp_contacts.size(); k++) {
			const t_contact_vertex& cv = qp_contacts[k];
			const t_pos4f& pnt = cv.get_point();
			const t_vec4f& nrm = cv.get_normal();

			printf("\t\t[q2p] k=%lu cv=<P=<%f,%f,%f> N=<%f,%f,%f>>\n", k, pnt.x(),pnt.y(),pnt.z(), nrm.x(),nrm.y(),nrm.z());
		}
	}
	#endif
}

void t_world_state::handle_colliding_objects(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	#if (ENABLE_FRICTION_FORCES == 1)
	apply_friction_forces(frame_num);
	#else
	(void) frame_num;
	#endif

	#if (ENABLE_REACTION_FORCES == 1)
	// first propagate action-reaction chains
	for (unsigned int n = 0; apply_reaction_forces(n); n++);
	#endif

	#if (ENABLE_COLLISION_FORCES == 1)
	// derive collision forces from momentum changes
	for (unsigned int n = 0; apply_collision_forces(n); n++);
	#endif

	// clear the indices of colliding pairs
	m_object_pair_indices.clear();
}



void t_world_state::apply_common_forces(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	for (size_t k = 0; k < (m_rigid_objects.size() - 1); k++) {
		m_rigid_objects[k].add_lin_force(m_global_scene->object_gravity_accell * m_rigid_objects[k].get_mass());
		m_rigid_objects[k].add_lin_force(m_rigid_objects[k].calc_air_drag_force());
	}

	// continuous acceleration (TODO: physics_scene parameter)
	// m_rigid_objects[0].add_lin_force(m_global_scene->object_tst_lin_accell * m_rigid_objects[0].get_mass() * (frame_num > 0));

	// time-limited acceleration
	// m_rigid_objects[0].add_lin_force(m_global_scene->object_tst_lin_accell * m_rigid_objects[0].get_mass() * (frame_num < 1000));

	// instantaneous acceleration
	m_rigid_objects[0].add_lin_force(m_global_scene->object_tst_lin_accell * 100.0f * m_rigid_objects[0].get_mass() * (frame_num == 0));
	m_rigid_objects[0].add_ang_force(m_global_scene->object_tst_ang_accell * 100.0f * m_rigid_objects[0].get_mass() * (frame_num == 0));
}

void t_world_state::apply_friction_forces(unsigned int frame_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	(void) frame_num;

	for (size_t k = 0; k < m_object_pair_indices.size(); k++) {
		#if 0
		if (std::max(m_object_pair_timings[k].first, m_object_pair_timings[k].second) < max_rollback_time)
			continue;
		#endif

		const unsigned int pair = m_object_pair_indices[k];
		const unsigned int p_id = pair / m_rigid_objects.size();
		const unsigned int q_id = pair % m_rigid_objects.size();

		t_rigid_body& p = m_rigid_objects[p_id];
		t_rigid_body& q = m_rigid_objects[q_id];

		t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
		t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

		t_rigid_body::apply_friction_forces(p, q, pq_contacts, qp_contacts);
	}
}

bool t_world_state::apply_collision_forces(unsigned int pass_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	bool fwd_pass = false;
	bool bck_pass = false;

	if (pass_num == 0) {
		// first pass; copy current velocities and forces to temporaries
		for (size_t k = 0; k < (m_rigid_objects.size() - 1); k++) {
			t_rigid_body& obj = m_rigid_objects[k];

			obj.set_lin_velocity(obj.get_lin_velocity(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);
			obj.set_ang_velocity(obj.get_ang_velocity(CURR_ATTRIB_IDX), PREV_ATTRIB_IDX);

			obj.store_forces();
			obj.clear_forces(CURR_ATTRIB_IDX);
		}
	}

	if (true) {
		// forward pass
		for (size_t k = 0; k < m_object_pair_indices.size(); k++) {
			#if 0
			// ignore collisions that occurred closer to the present than max_rollback
			if (std::max(m_object_pair_timings[k].first, m_object_pair_timings[k].second) < max_rollback_time)
				continue;
			#endif

			const unsigned int pair = m_object_pair_indices[k];
			const unsigned int p_id = pair / m_rigid_objects.size();
			const unsigned int q_id = pair % m_rigid_objects.size();

			t_rigid_body& p = m_rigid_objects[p_id];
			t_rigid_body& q = m_rigid_objects[q_id];

			t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
			t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

			fwd_pass |= t_rigid_body::apply_collision_forces(p, q, pq_contacts, qp_contacts);
		}

		fwd_pass &= (pass_num < MAX_FORCE_EXCHANGE_ITERS);
	}

	if (fwd_pass) {
		// backward pass
		for (size_t k = 0; k < m_object_pair_indices.size(); k++) {
			#if 0
			// ignore collisions that occurred closer to the present than max_rollback
			if (std::max(m_object_pair_timings[k].first, m_object_pair_timings[k].second) < max_rollback_time)
				continue;
			#endif

			const unsigned int pair = m_object_pair_indices[(m_object_pair_indices.size() - 1) - k];
			const unsigned int p_id = pair / m_rigid_objects.size();
			const unsigned int q_id = pair % m_rigid_objects.size();

			t_rigid_body& p = m_rigid_objects[p_id];
			t_rigid_body& q = m_rigid_objects[q_id];

			t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
			t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

			bck_pass |= t_rigid_body::apply_collision_forces(p, q, pq_contacts, qp_contacts);
		}

		bck_pass &= (pass_num < MAX_FORCE_EXCHANGE_ITERS);
	}

	if (!fwd_pass && !bck_pass) {
		// last pass; copy back the temporaries
		for (size_t k = 0; k < (m_rigid_objects.size() - 1); k++) {
			t_rigid_body& obj = m_rigid_objects[k];

			obj.set_lin_velocity(obj.get_lin_velocity(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
			obj.set_ang_velocity(obj.get_ang_velocity(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);

			obj.set_lin_force(obj.get_lin_force(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
			obj.set_ang_force(obj.get_ang_force(PREV_ATTRIB_IDX), CURR_ATTRIB_IDX);
		}
	}

	return (fwd_pass || bck_pass);
}

bool t_world_state::apply_reaction_forces(unsigned int pass_num) {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	bool fwd_pass = false;
	bool bck_pass = false;

	if (true) {
		// forward pass
		for (size_t k = 0; k < m_object_pair_indices.size(); k++) {
			#if 0
			if (std::max(m_object_pair_timings[k].first, m_object_pair_timings[k].second) < max_rollback_time)
				continue;
			#endif

			const unsigned int pair = m_object_pair_indices[k];
			const unsigned int p_id = pair / m_rigid_objects.size();
			const unsigned int q_id = pair % m_rigid_objects.size();

			t_rigid_body& p = m_rigid_objects[p_id];
			t_rigid_body& q = m_rigid_objects[q_id];

			t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
			t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

			fwd_pass |= t_rigid_body::apply_reaction_forces(p, q, pq_contacts, qp_contacts);
		}

		fwd_pass &= (pass_num < MAX_FORCE_EXCHANGE_ITERS);
	}

	if (fwd_pass) {
		// backward pass
		// this is present in case the order in which we found colliding
		// pairs would lead to an inefficient force propagation pattern
		for (size_t k = 0; k < m_object_pair_indices.size(); k++) {
			#if 0
			if (std::max(m_object_pair_timings[k].first, m_object_pair_timings[k].second) < max_rollback_time)
				continue;
			#endif

			const unsigned int pair = m_object_pair_indices[(m_object_pair_indices.size() - 1) - k];
			const unsigned int p_id = pair / m_rigid_objects.size();
			const unsigned int q_id = pair % m_rigid_objects.size();

			t_rigid_body& p = m_rigid_objects[p_id];
			t_rigid_body& q = m_rigid_objects[q_id];

			t_contacts_vector& pq_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].first;
			t_contacts_vector& qp_contacts = m_object_contacts[CURR_ATTRIB_IDX][pair].second;

			bck_pass |= t_rigid_body::apply_reaction_forces(p, q, pq_contacts, qp_contacts);
		}

		// we should never need more passes than we have objects
		// (in practice a handful are sufficient even for stacks
		// but with USE_AVERAGE_CONTACT_VERTEX = 0 the number is
		// unpredictable)
		bck_pass &= (pass_num < MAX_FORCE_EXCHANGE_ITERS);
	}

	return (fwd_pass || bck_pass);
}



void t_world_state::run_headless(const t_arg_parser* arg_parser) {
	const unsigned int num_objects = arg_parser->get_int_arg("num_objects");
	const unsigned int num_sframes = arg_parser->get_int_arg("max_sframes");
	const unsigned int scene_index = arg_parser->get_int_arg("scene_index");

	load_scene(scene_index, num_objects);

	for (unsigned int n = 0; n < num_sframes; n++) {
		update_objects(n, 1.0f);
	}
}



const t_rigid_body* t_world_state::trace_ray(const t_ray raw_ray, t_pos4f* hit_pos) const {
	const t_rigid_body* hit_obj = nullptr;

	// if no intersection, .w remains 0
	t_pos4f tmp_int = t_pos4f::null_point();
	t_pos4f cur_int = t_pos4f::null_point();

	// this ensures the first actual intersection is always closer
	cur_int.x() = M_FINF;

	std::vector<unsigned int> tests;
	std::vector<t_object_grid::t_grid_cell> cells;

	if (m_object_grid.get_cells_on_ray(raw_ray, cells) == 0)
		return hit_obj;

	tests.resize(m_rigid_objects.size(), 0);

	// make the ray infinitely long
	const t_ray inf_ray = {raw_ray.pos(), raw_ray.dir()};

	for (const t_object_grid::t_grid_cell& cell: cells) {
		const t_object_grid::t_object_vector& cell_objects = cell.get_objects();

		for (const t_grid_object* object: cell_objects) {
			if (tests[object->get_unique_id()] != 0)
				continue;

			tests[object->get_unique_id()] += 1;

			// convert proxy back to actual object
			const t_rigid_body& body = m_rigid_objects[object->get_unique_id()];
			const t_collision_mesh& mesh = body.get_collision_mesh();
			const t_mat44f& mat = body.get_transform_mat();
			const t_mat44f& inv = mat.invert_affine();

			// .w can be 0 or 1; ignored for distance-calculation
			tmp_int = mesh.calc_ray_intersection(inf_ray, inv);

			if (tmp_int.w() == 0.0f)
				continue;
			if ((tmp_int - inf_ray.pos()).sq_len() > (cur_int - inf_ray.pos()).sq_len())
				continue;

			cur_int = tmp_int;
			hit_obj = &body;
		}
	}

	if (hit_pos != nullptr)
		*hit_pos = cur_int;

	return hit_obj;
}

