#include <cstdio> // printf

#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL/SDL.h>

#include "global_defs.hpp"
#include "program_clock.hpp"
#include "program_state.hpp"
#include "world_state.hpp"
#include "misc/arg_parser.hpp"
#include "misc/highres_clock.hpp"
#include "misc/opengl_wrapper.hpp"
#include "misc/scoped_timer.hpp"


enum {
	CAM_ACTION_MOVE_FWD = 0, // move forward
	CAM_ACTION_MOVE_BCK = 1, // move backward
	CAM_ACTION_STRA_LFT = 2, // strafe left
	CAM_ACTION_STRA_RGT = 3, // strafe right
	CAM_ACTION_STRA_UP  = 4, // strafe up
	CAM_ACTION_STRA_DWN = 5, // strafe down
	CAM_ACTION_TURN_LFT = 6, // yaw left
	CAM_ACTION_TURN_RGT = 7, // yaw right
	CAM_ACTION_TURN_DWN = 8, // pitch down
	CAM_ACTION_TURN_UP  = 9, // pitch up
};

enum {
	GUI_ACTION_QUIT_SIM  = 10,
	GUI_ACTION_FULL_SCR  = 11,
	GUI_ACTION_TRCK_NXT  = 12,
	GUI_ACTION_TRCK_PRV  = 13,
	GUI_ACTION_DRAW_GRND = 14,
	GUI_ACTION_DRAW_OBJS = 15,
	GUI_ACTION_DRAW_CTCS = 16,
	GUI_ACTION_ATTR_INDX = 17,
};

enum {
	SIM_ACTION_INCR_TICK_R = 18,
	SIM_ACTION_DECR_TICK_R = 19,
	SIM_ACTION_INVERT_TIME = 20,
	SIM_ACTION_EXEC_PAUSE  = 21,
	SIM_ACTION_EXEC_FRAME  = 22,
	SIM_ACTION_RESET_SCENE = 23,
	SIM_ACTION_CYCLE_SCENE = 24,
};

static const t_key_state key_defs[] = {
	{SDLK_w    , CAM_ACTION_MOVE_FWD, 0,  1.0f},
	{SDLK_s    , CAM_ACTION_MOVE_BCK, 0, -1.0f},
	{SDLK_a    , CAM_ACTION_STRA_LFT, 0, -1.0f},
	{SDLK_d    , CAM_ACTION_STRA_RGT, 0,  1.0f},
	{SDLK_e    , CAM_ACTION_STRA_UP , 0,  1.0f},
	{SDLK_q    , CAM_ACTION_STRA_DWN, 0, -1.0f},
	{SDLK_LEFT , CAM_ACTION_TURN_LFT, 0, -1.0f},
	{SDLK_RIGHT, CAM_ACTION_TURN_RGT, 0,  1.0f},
	{SDLK_UP   , CAM_ACTION_TURN_DWN, 0, -1.0f},
	{SDLK_DOWN , CAM_ACTION_TURN_UP , 0,  1.0f},

	{SDLK_ESCAPE      , GUI_ACTION_QUIT_SIM , 0, 0.0f},
	{SDLK_F1          , GUI_ACTION_FULL_SCR , 0, 0.0f},
	{SDLK_t           , GUI_ACTION_TRCK_NXT , 0, 0.0f},
	{SDLK_y           , GUI_ACTION_TRCK_PRV , 0, 0.0f},
	{SDLK_BACKSLASH   , GUI_ACTION_DRAW_GRND, 0, 0.0f},
	{SDLK_LEFTBRACKET , GUI_ACTION_DRAW_OBJS, 0, 0.0f},
	{SDLK_RIGHTBRACKET, GUI_ACTION_DRAW_CTCS, 0, 0.0f},
	{SDLK_RETURN      , GUI_ACTION_ATTR_INDX, 0, 0.0f},

	{SDLK_EQUALS   , SIM_ACTION_INCR_TICK_R, 0, 0.0f},
	{SDLK_MINUS    , SIM_ACTION_DECR_TICK_R, 0, 0.0f},
	{SDLK_BACKSPACE, SIM_ACTION_INVERT_TIME, 0, 0.0f},
	{SDLK_p        , SIM_ACTION_EXEC_PAUSE , 0, 0.0f},
	{SDLK_o        , SIM_ACTION_EXEC_FRAME , 0, 0.0f},
	{SDLK_l        , SIM_ACTION_RESET_SCENE, 0, 0.0f},
	{SDLK_i        , SIM_ACTION_CYCLE_SCENE, 0, 0.0f},
};



struct t_camera {
public:
	void move(const t_vec4f& vec, float tgt_move = 1.0f) {
		pos += (vec           );
		tgt += (vec * tgt_move);
	}

	void strafe(const t_vec4f& vec, float tgt_move, float tgt_dist) {
		move(vec, tgt_move);
		set(pos, tgt, tgt_dist);
	}

	// yaw, pitch
	void turn(const t_vec4f& vec, float tgt_move, float tgt_dist) {
		set(tgt += (vec * tgt_move), pos, tgt_dist);
	}

	void set(t_pos4f& dst, const t_pos4f& src, float tgt_dist) {
		dst = src + (dst - src).normalize() * tgt_dist;
	}

	t_vec4f get_tgt_vec() const { return  (tgt - pos);              }
	t_vec4f get_tgt_dir() const { return ((tgt - pos).normalize()); }

public:
	t_pos4f  pos;
	t_pos4f  tgt;
	t_vec4f upv;
};






struct t_opengl_light {
public:
	t_pos4f pos_world;
	t_tup4f amb_color;
	t_tup4f dif_color;
};

struct t_opengl_material {
	t_tup4f amb_refl;
	t_tup4f dif_refl;
};

struct t_render_mesh {
	std::vector<t_gl_wrapper_state::t_vertex> m_verts;
	std::vector<t_gl_wrapper_state::t_index> m_inds;
};


static const t_opengl_light default_light = {
	{2.0f,  3.0f, 4.0f, 0.0f}, // w=0 (directional)
	{0.3f,  0.3f, 0.3f, 1.0f},
	{1.0f,  1.0f, 1.0f, 1.0f},
};

static const t_opengl_material default_material = {
	{1.0f, 1.0f, 1.0f, 1.0f},
	{1.0f, 1.0f, 1.0f, 1.0f},
};
static const t_opengl_material ground_material = {
	{0.5f, 0.5f, 0.75f, 1.0f},
	{0.5f, 0.5f, 0.75f, 1.0f},
};
static const t_opengl_material mesh_materials[] = {
	{{1.0f, 0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f, 1.0f}}, // red
	{{0.0f, 1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f, 1.0f}}, // green
	{{0.0f, 0.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 1.0f, 1.0f}}, // blue
	{{1.0f, 1.0f, 0.0f, 1.0f}, {1.0f, 1.0f, 0.0f, 1.0f}}, // yellow
	{{1.0f, 0.0f, 1.0f, 1.0f}, {1.0f, 0.0f, 1.0f, 1.0f}}, // purple
	{{0.0f, 1.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 1.0f, 1.0f}}, // cyan
	{{1.0f, 1.0f, 1.0f, 1.0f}, {1.0f, 1.0f, 1.0f, 1.0f}}, // white
	{{0.0f, 0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 0.0f, 1.0f}}, // black
	{{0.5f, 0.5f, 0.5f, 1.0f}, {0.5f, 0.5f, 0.5f, 1.0f}}, // grey
};

static t_camera default_camera = {
	{10.0f + 5.0f, 5.0f, 10.0f + 5.0f, 1.0f},
	{ 0.0f + 0.0f, 0.0f,  0.0f + 0.0f, 1.0f},
	{ 0.0f,        1.0f,  0.0f,        0.0f},
};

// [0] := ground, [1] := cube
static t_render_mesh render_meshes[2];




// already using namespace newtonics via world_state.hpp, etc
static const t_rigid_body* tracked_object = nullptr;
static const t_rigid_body* rtraced_object = nullptr;

static bool draw_ground   = true;
static bool draw_objects  = true;
static bool draw_contacts = true;
static bool prev_attribs  = true; // when paused



static void create_render_mesh(const t_collision_mesh& col_mesh, t_render_mesh& rdr_mesh) {
	const std::vector<t_idx_tri>& idx_tris = col_mesh.get_polys();

	rdr_mesh.m_verts.resize(col_mesh.get_num_verts());
	rdr_mesh.m_inds.resize(col_mesh.get_num_polys() * 3);

	t_gl_wrapper_state::t_vertex* verts = &rdr_mesh.m_verts[0];
	t_gl_wrapper_state::t_index* inds = &rdr_mesh.m_inds[0];

	for (size_t i = 0; i < col_mesh.get_num_verts(); i++) {
		const t_pos4f& p = col_mesh.get_vertex(i);

		verts[i].pxyzw[0] = p.x();
		verts[i].pxyzw[1] = p.y();
		verts[i].pxyzw[2] = p.z();
		verts[i].pxyzw[3] = p.w();
	}

	for (size_t i = 0; i < idx_tris.size(); i++) {
		const t_tup4ui& idx = idx_tris[i].get_indices();
		const t_vec4f& n = idx_tris[i].get_normal();

		for (unsigned int k = 0; k < 3; k++) {
			t_gl_wrapper_state::t_vertex& v = verts[ idx[k] ];

			v.nxyz[0] = n.x();
			v.nxyz[1] = n.y();
			v.nxyz[2] = n.z();
		}

		// account for windings
		if ((i & 1) != 0) {
			inds[i * 3 + 0] = idx[2];
			inds[i * 3 + 1] = idx[1];
			inds[i * 3 + 2] = idx[0];
		} else {
			inds[i * 3 + 0] = idx[0];
			inds[i * 3 + 1] = idx[1];
			inds[i * 3 + 2] = idx[2];
		}
	}
}

#if 0
static void render_mesh_imm_raw(const t_collision_mesh& mesh) {
	const std::vector<t_idx_tri>& idx_tris = mesh.get_polys();

	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i < idx_tris.size(); i++) {
		const t_tup4ui& tri_indices = idx_tris[i].get_indices();
		const t_vec4f& tri_normal = idx_tris[i].get_normal();

		const t_pos4f& v0 = mesh.get_vertex(tri_indices[0]);
		const t_pos4f& v1 = mesh.get_vertex(tri_indices[1]);
		const t_pos4f& v2 = mesh.get_vertex(tri_indices[2]);

		glNormal3f(tri_normal.x(), tri_normal.y(), tri_normal.z());
		glVertex3f(v0.x(), v0.y(), v0.z());
		glVertex3f(v1.x(), v1.y(), v1.z());
		glVertex3f(v2.x(), v2.y(), v2.z());
	}

	glEnd();
}
#endif

#if 0
static void render_mesh_ind_raw(const t_collision_mesh& mesh) {
	const std::vector<t_idx_tri>& idx_tris = mesh.get_polys();

	t_gl_wrapper_state::t_vertex* verts = gl.GetArrayPointer();
	t_gl_wrapper_state::t_index* inds = gl.GetIndexPointer();

	// not the right memory layout
	// const unsigned int* inds = (idx_tris[0].get_indices()).xyzw();

	// construct verts and inds on-the-fly, using the GL state buffers
	// can not loop over the verts directly; we want triangle normals
	for (size_t i = 0; i < idx_tris.size(); i++) {
		const t_idx_tri& tri = idx_tris[i];
		const t_tup4ui& idx = tri.get_indices();
		const t_vec4f& n = tri.get_normal();

		for (unsigned int k = 0; k < 3; k++) {
			const t_pos4f& v = mesh.get_vertex(idx[k]);

			verts[ idx[k] ].pxyzw[0] = v.x();
			verts[ idx[k] ].pxyzw[1] = v.y();
			verts[ idx[k] ].pxyzw[2] = v.z();
			verts[ idx[k] ].pxyzw[3] = v.w();
			verts[ idx[k] ].nxyz [0] = n.x();
			verts[ idx[k] ].nxyz [1] = n.y();
			verts[ idx[k] ].nxyz [2] = n.z();
		}

		// account for windings
		if ((i & 1) != 0) {
			inds[i * 3 + 0] = idx[2];
			inds[i * 3 + 1] = idx[1];
			inds[i * 3 + 2] = idx[0];
		} else {
			inds[i * 3 + 0] = idx[0];
			inds[i * 3 + 1] = idx[1];
			inds[i * 3 + 2] = idx[2];
		}
	}

	gl.DrawElements(GL_TRIANGLES, mesh.get_num_polys() * 3, verts, inds);
}
#endif




t_program_state::t_program_state(const t_arg_parser* arg_parser):
	m_prog_clock(nullptr),
	m_world_state(nullptr),

	m_window_size_x(0),
	m_window_size_y(0),

	m_num_sim_frames(0),
	m_max_sim_frames(0),
	m_num_draw_calls(0),

	m_is_running(true),
	m_is_paused(true),
	m_is_visible(true),
	m_is_inited(false),

	m_has_mouse_focus(true),
	m_has_keyboard_focus(true)
{
	m_window_size_x = arg_parser->get_int_arg("win_size_x");
	m_window_size_y = arg_parser->get_int_arg("win_size_y");

	m_num_sim_frames = 0;
	m_max_sim_frames = arg_parser->get_int_arg("max_sframes");

	printf("[pr::%s]\n", __func__);
	printf("\tENABLE_FULLSCREEN_MODE    = %u\n", ENABLE_FULLSCREEN_MODE);
	printf("\tENABLE_HIGHRES_CLOCK      = %u\n", ENABLE_HIGHRES_CLOCK);
	printf("\tENABLE_SCOPED_TIMERS      = %u\n", ENABLE_SCOPED_TIMERS);
	printf("\tENABLE_TIMER_THREAD       = %u\n", ENABLE_TIMER_THREAD);
	printf("\n");
	printf("\tTIME_STEP_SIZE_MS         = %u\n", TIME_STEP_SIZE_MS);
	printf("\tTIME_STEP_SIZE_S          = %f\n", TIME_STEP_SIZE_S);
	printf("\n");
	printf("\tENABLE_REACTION_FORCES    = %u\n", ENABLE_REACTION_FORCES);
	printf("\tENABLE_COLLISION_FORCES   = %u\n", ENABLE_COLLISION_FORCES);
	printf("\tENABLE_FRICTION_FORCES    = %u\n", ENABLE_FRICTION_FORCES);
	printf("\n");
	printf("\tCLAMP_FORCE_RESPONSES     = %u\n", CLAMP_FORCE_RESPONSES);
	printf("\tRANDOMIZE_OBJECT_ORDER    = %u\n", RANDOMIZE_OBJECT_ORDER);
	printf("\n");
	printf("\tMAX_OBJECT_COUNT          = %u\n", MAX_OBJECT_COUNT);
	printf("\tMAX_OBJECT_CONTACT_AGE    = %u\n", MAX_OBJECT_CONTACT_AGE);
	printf("\tMAX_FORCE_EXCHANGE_ITERS  = %u\n", MAX_FORCE_EXCHANGE_ITERS);
	printf("\n");
}

void t_program_state::init(t_world_state* wstate, t_program_clock* pclock) {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("[pr::%s] failed to initialize SDL: %s\n", __func__, SDL_GetError());
		exit(1);
	}

	// make sure SDL_Quit always gets called
	atexit(SDL_Quit);

	#if (ENABLE_FULLSCREEN_MODE == 1)
	const uint32_t vid_flags = SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_OPENGLBLIT | SDL_FULLSCREEN;
	#else
	const uint32_t vid_flags = SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_OPENGLBLIT;
	#endif

	if (!create_window(16, vid_flags)) {
		printf("[pr::%s] failed to create screen surface: %s\n", __func__, SDL_GetError());
		exit(1);
	}

	SDL_WM_SetCaption("newtonics", nullptr);
	// setting this prevents having to poll the key-state array
	SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY, SDL_DEFAULT_REPEAT_INTERVAL);

	{
		glClearColor(0.2f, 0.2f, 0.2f, 0.5f);
		glClearDepth(1.0f);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glLightfv(GL_LIGHT0, GL_AMBIENT, default_light.amb_color.xyzw());
		glLightfv(GL_LIGHT0, GL_DIFFUSE, default_light.dif_color.xyzw());
		glShadeModel(GL_FLAT);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	}

	m_world_state = wstate;
	m_prog_clock = pclock;
	m_prog_clock->reset(GET_CURRENT_TICK());

	{
		m_key_states.fill({0, 0, 0, 0.0f});
		m_key_rbinds.fill(nullptr);

		assert(SDLK_LAST < m_key_states.size());

		// initialize key bindings
		for (size_t n = 0; n < (sizeof(key_defs) / sizeof(t_key_state)); n++) {
			const t_key_state& def_key = key_defs[n];

			const uint32_t sdl_symb = def_key.symb;
			const uint32_t mod_symb = def_key.flag;

			t_key_state*  sdl_key = &m_key_states[sdl_symb];
			t_key_state** mod_key = &m_key_rbinds[mod_symb];

			sdl_key->symb = mod_symb;
			sdl_key->sign = def_key.sign;

			*mod_key = sdl_key;
		}
	}

	{
		// create render-meshes; all rigid objects share the same
		// const t_planar_body& ground_obj = m_world_state->get_ground();
		const t_rigid_body& ground_obj = m_world_state->get_ground();
		const t_rigid_body& rigid_obj = m_world_state->get_rigid_objects()[0];

		create_render_mesh(ground_obj.get_collision_mesh(), render_meshes[0]);
		create_render_mesh(rigid_obj.get_collision_mesh(), render_meshes[1]);
	}

	m_is_inited = true;
}

void t_program_state::kill() {
	SDL_Event event;
	event.type = SDL_QUIT;

	if (SDL_PushEvent(&event) == -1) {
		printf("[pr::%s] failed to push SDL_QUIT event: %s\n", __func__, SDL_GetError());
		exit(1);
	}
}


void t_program_state::update_keys(const uint8_t* sdl_key_state) {
	// always points to the same memory
	assert(sdl_key_state != nullptr);

	// copy states of all SDL keys we are interested in
	for (size_t n = 0; n < (sizeof(key_defs) / sizeof(t_key_state)); n++) {
		m_key_states[ key_defs[n].symb ].flag = sdl_key_state[ key_defs[n].symb ];
	}
}

void t_program_state::handle_key_press(const SDL_keysym& ks) {
	assert(ks.sym < m_key_states.size());

	// these keys do not need to be reacted to every frame
	switch (m_key_states[ks.sym].symb) {
		case GUI_ACTION_QUIT_SIM: {
			kill();
		} break;
		case GUI_ACTION_FULL_SCR: {
			toggle_fscr();
		} break;

		case GUI_ACTION_TRCK_NXT: {
			update_tracked_object(m_world_state->get_rigid_objects(), 1);
		} break;
		case GUI_ACTION_TRCK_PRV: {
			update_tracked_object(m_world_state->get_rigid_objects(), -1);
		} break;

		case GUI_ACTION_DRAW_GRND: {
			draw_ground = !draw_ground;
		} break;
		case GUI_ACTION_DRAW_OBJS: {
			draw_objects = !draw_objects;
		} break;
		case GUI_ACTION_DRAW_CTCS: {
			draw_contacts = !draw_contacts;
		} break;
		case GUI_ACTION_ATTR_INDX: {
			prev_attribs = !prev_attribs;
		} break;

		case SIM_ACTION_EXEC_PAUSE: {
			// toggle pause
			m_prog_clock->scale_acc_ticks(m_is_paused = !m_is_paused);
		} break;
		case SIM_ACTION_EXEC_FRAME: {
			// single-step a frame
			m_prog_clock->add_acc_ticks(m_prog_clock->get_time_step_size() * m_is_paused);
		} break;

		case SIM_ACTION_RESET_SCENE:
		case SIM_ACTION_CYCLE_SCENE: {
			update_scene(m_key_states[ks.sym].symb == SIM_ACTION_CYCLE_SCENE);
		} break;

		case SIM_ACTION_INCR_TICK_R:
		case SIM_ACTION_DECR_TICK_R: {
			m_prog_clock->scale_tick_rate(lerp(2.0f, 0.5f, float(m_key_states[ks.sym].symb == SIM_ACTION_DECR_TICK_R)));

			#if (ENABLE_PROGSTATE_DEBUG_PRINTS == 1)
			printf("[pr::%s][+] tick_rate=%f\n", __func__, m_prog_clock->get_tick_rate_mult());
			#endif
		} break;
		case SIM_ACTION_INVERT_TIME: {
			// does not really work, forces need to inverted instead
			m_prog_clock->scale_time_step(-1.0f);
		} break;

		default: {
		} break;
	}
}

void t_program_state::handle_key_release(const SDL_keysym& ks) {
	switch (ks.sym) {
		default: {
		} break;
	}
}


void t_program_state::update_camera(t_camera& cam, float dif_time_ssecs) {
	float tgt_move = 1.0f;
	float tgt_dist = 0.0f;

	if (tracked_object != nullptr) {
		const t_rigid_body& object = *tracked_object;
		const t_collision_mesh& mesh = object.get_collision_mesh();
		const t_bounding_box& bbox = mesh.get_bounding_box();

		// focus on tracked object
		cam.tgt  =  object.get_position();
		cam.tgt += (object.get_lin_velocity() * m_prog_clock->get_tick_extr_fact());

		tgt_move = 0.0f;
		tgt_dist = bbox.get_diameter();
	}

	rtraced_object = m_world_state->trace_ray(t_ray(cam.pos, cam.tgt));
	#if (0 == 1 && ENABLE_PROGSTATE_DEBUG_PRINTS == 1)
	printf("[pr::%s] rtraced_object_id=%u\n", __func__, (rtraced_object != nullptr)? rtraced_object->get_id(): -1u);
	#endif

	// TODO
	// cam.update(dif_time_ssecs, tgt_move, tgt_dist);

	const t_vec4f tgt_vec   = cam.get_tgt_vec();
	const t_vec4f tgt_dir   = cam.get_tgt_dir();
	// const t_vec4f pan_vec_h = t_vec4f::x_axis_vector() * dif_time_ssecs;
	const t_vec4f pan_vec_v = t_vec4f::y_axis_vector() * dif_time_ssecs;
	const t_vec4f pan_vec_h = (tgt_vec.outer_product(t_vec4f::y_axis_vector())).normalize() * dif_time_ssecs;
	// const t_vec4f pan_vec_v = (tgt_vec.outer_product(t_vec4f::x_axis_vector())).normalize() * dif_time_ssecs;

	const float cur_dist = tgt_vec.inner_product(tgt_dir);

	// move forward
	if (m_key_rbinds[CAM_ACTION_MOVE_FWD]->flag != 0)
		cam.move(tgt_dir * m_key_rbinds[CAM_ACTION_MOVE_FWD]->sign * dif_time_ssecs * (cur_dist > tgt_dist));
	// move backward
	if (m_key_rbinds[CAM_ACTION_MOVE_BCK]->flag != 0)
		cam.move(tgt_dir * m_key_rbinds[CAM_ACTION_MOVE_BCK]->sign * dif_time_ssecs);

	// strafe left; orbits at fixed distance when tracking
	if (m_key_rbinds[CAM_ACTION_STRA_LFT]->flag != 0)
		cam.strafe(pan_vec_h * m_key_rbinds[CAM_ACTION_STRA_LFT]->sign, tgt_move, cur_dist);
	// strafe right; orbits at fixed distance when tracking
	if (m_key_rbinds[CAM_ACTION_STRA_RGT]->flag != 0)
		cam.strafe(pan_vec_h * m_key_rbinds[CAM_ACTION_STRA_RGT]->sign, tgt_move, cur_dist);

	// move down (vertical strafe; semi-orbits when tracking)
	if (m_key_rbinds[CAM_ACTION_STRA_DWN]->flag != 0)
		cam.move(t_vec4f::y_axis_vector() * m_key_rbinds[CAM_ACTION_STRA_DWN]->sign * dif_time_ssecs, tgt_move);
	// move up (vertical strafe; semi-orbits when tracking)
	if (m_key_rbinds[CAM_ACTION_STRA_UP]->flag != 0)
		cam.move(t_vec4f::y_axis_vector() * m_key_rbinds[CAM_ACTION_STRA_UP]->sign * dif_time_ssecs, tgt_move);

	// own-axis yaw (disabled when tracking)
	if (m_key_rbinds[CAM_ACTION_TURN_LFT]->flag != 0)
		cam.turn(pan_vec_h * m_key_rbinds[CAM_ACTION_TURN_LFT]->sign, tgt_move, cur_dist);
	if (m_key_rbinds[CAM_ACTION_TURN_RGT]->flag != 0)
		cam.turn(pan_vec_h * m_key_rbinds[CAM_ACTION_TURN_RGT]->sign, tgt_move, cur_dist);

	// own-axis pitch (disabled when tracking)
	if (m_key_rbinds[CAM_ACTION_TURN_DWN]->flag != 0)
		cam.turn(pan_vec_v * m_key_rbinds[CAM_ACTION_TURN_DWN]->sign, tgt_move, cur_dist);
	if (m_key_rbinds[CAM_ACTION_TURN_UP]->flag != 0)
		cam.turn(pan_vec_v * m_key_rbinds[CAM_ACTION_TURN_UP]->sign, tgt_move, cur_dist);
}

void t_program_state::update_scene(uint32_t scene_offset) {
	// reset or cycle scene
	m_world_state->reload_scene(scene_offset);
	m_prog_clock->reset(GET_CURRENT_TICK());

	m_num_sim_frames = 0;
	m_num_draw_calls = 0;

	// leave to user
	// m_is_paused = true;

	tracked_object = nullptr;
}

void t_program_state::update_tracked_object(const std::vector<newtonics::lib_physics::t_rigid_body>& objects, int32_t dir) {
	const size_t init_id = (objects.size() - 1) * (dir == -1);

	// enable tracking
	if (tracked_object == nullptr) {
		tracked_object = &objects[init_id];
		return;
	}

	const size_t next_id = tracked_object->get_id() + dir;

	// cycle to next object or disable
	if (next_id >= objects.size()) {
		tracked_object = nullptr;
		return;
	}

	tracked_object = &objects[next_id];
}

bool t_program_state::tick_clock(t_program_clock* clock) {
	return ((m_num_draw_calls *= (1 - clock->tick(m_is_paused))) == 0);
}

void t_program_state::update() {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	#if (ENABLE_TIMER_THREAD == 0)
	if (tick_clock(m_prog_clock))
		print_fps();
	#endif

	update_keys(SDL_GetKeyState(nullptr));
	update_camera(default_camera, m_prog_clock->get_dif_time_ssecs() * 0.33f);

	const uint64_t time_step_size = m_prog_clock->get_time_step_size();
	const uint64_t num_time_steps = m_prog_clock->get_num_time_steps();

	// fixed-timestep physics; ticks are not accumulated while paused
	for (uint64_t n = 0; (n < num_time_steps  &&  m_num_sim_frames < m_max_sim_frames); n++) {
		#if (ENABLE_PROGSTATE_DEBUG_PRINTS == 1)
		printf("[pr::%s][f=%u]\n", __func__, m_num_sim_frames);
		#endif

		m_world_state->update_objects(m_num_sim_frames, 1.0f);
		m_prog_clock->sub_acc_ticks(time_step_size);

		#if (ENABLE_PROGSTATE_DEBUG_PRINTS == 1)
		printf("\n");
		#endif

		m_num_sim_frames += 1;
	}

	m_num_draw_calls += 1;
}

void t_program_state::render() {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	// when paused, grab the previous position and contact-data
	// this is more useful when single-stepping since after each
	// simulated frame the objects are at pos(t)+vel(t) while the
	// contacts were calculated at pos(t), ie. one object::update
	// step behind
	const unsigned int attrib_index = lerp(int(CURR_ATTRIB_IDX), int(PREV_ATTRIB_IDX), int(prev_attribs));

	// const t_planar_body& ground_body = m_world_state->get_ground();
	const t_rigid_body& ground_body = m_world_state->get_ground();
	const std::vector<t_rigid_body>& objects = m_world_state->get_rigid_objects();
	const std::vector<t_contacts_pair>& contacts = m_world_state->get_object_contacts(attrib_index);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// set up-vector along y-dir s.t. world +x and +y correspond to window +x and +y
	gluLookAt(
		default_camera.pos.x(), default_camera.pos.y(), default_camera.pos.z(),
		default_camera.tgt.x(), default_camera.tgt.y(), default_camera.tgt.z(),
		default_camera.upv.x(), default_camera.upv.y(), default_camera.upv.z());
	glLightfv(GL_LIGHT0, GL_POSITION, default_light.pos_world.xyzw());

	gl.EnableClientArrays();

	if (draw_ground) {
		glPushMatrix();
		glMultMatrixf((ground_body.get_transform_mat()).get_values());

			// draw the ground mesh, filled
			glMaterialfv(GL_FRONT, GL_AMBIENT, ground_material.amb_refl.xyzw());
			glMaterialfv(GL_FRONT, GL_DIFFUSE, ground_material.dif_refl.xyzw());

			gl.DrawElements(GL_TRIANGLES, render_meshes[0].m_inds.size(), &render_meshes[0].m_verts[0], &render_meshes[0].m_inds[0]);

			// draw the ground mesh again, black-outlined
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(2.0f);

			glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[7].amb_refl.xyzw());
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[7].dif_refl.xyzw());

			glScalef(1.0f + (1.0f / 1024.0f), 1.0f + (1.0f / 1024.0f), 1.0f + (1.0f / 1024.0f));
			gl.DrawElements(GL_TRIANGLES, render_meshes[0].m_inds.size(), &render_meshes[0].m_verts[0], &render_meshes[0].m_inds[0]);

			glLineWidth(1.0f);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glPopMatrix();
	}

	if (draw_objects) {
		// draw object meshes excluding ground (objects.back()), filled
		for (uint32_t n = 0; n < (objects.size() - 1); n++) {
			glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[n & 7].amb_refl.xyzw());
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[n & 7].dif_refl.xyzw());

			const t_mat44f& mat = objects[n].get_transform_mat(attrib_index);
			const t_vec4f& vec = objects[n].get_lin_velocity(attrib_index) * m_prog_clock->get_tick_extr_fact() * (1 - m_is_paused);

			glPushMatrix();
			glMultMatrixf(mat.get_values());
			glTranslatef(vec.x(), vec.y(), vec.z());

			gl.DrawElements(GL_TRIANGLES, render_meshes[1].m_inds.size(), &render_meshes[1].m_verts[0], &render_meshes[1].m_inds[0]);

			glPopMatrix();
		}

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2.0f);
		glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[7].amb_refl.xyzw());
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[7].dif_refl.xyzw());

		// draw object meshes again, black-outlined
		for (uint32_t n = 0; n < (objects.size() - 1); n++) {
			const t_mat44f& mat = objects[n].get_transform_mat(attrib_index);
			const t_vec4f& vec = objects[n].get_lin_velocity(attrib_index) * m_prog_clock->get_tick_extr_fact() * (1 - m_is_paused);

			glPushMatrix();
			glMultMatrixf(mat.get_values());
			glTranslatef(vec.x(), vec.y(), vec.z());
			glScalef(1.0f + (1.0f / 1024.0f), 1.0f + (1.0f / 1024.0f), 1.0f + (1.0f / 1024.0f));

			gl.DrawElements(GL_TRIANGLES, render_meshes[1].m_inds.size(), &render_meshes[1].m_verts[0], &render_meshes[1].m_inds[0]);

			glPopMatrix();
		}

		glLineWidth(1.0f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	if (draw_contacts) {
		glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[8].amb_refl.xyzw());
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[8].dif_refl.xyzw());

		for (uint32_t pair = 0; pair < contacts.size(); pair++) {
			// typeof(contacts         ) := std::vector<t_contacts_pair>
			// typeof(contacts[i]      ) := std::pair<t_contacts_vector, t_contacts_vector> (t_contacts_pair)
			// typeof(contacts[i].first) := std::vector<t_contact_vertex> (t_contacts_vector)
			const t_contacts_pair& c = contacts[pair];
			const t_contacts_vector& vf = c.first;
			const t_contacts_vector& vs = c.second;

			if (!vf.empty()) {
				// glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[7].amb_refl.xyzw());
				// glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[7].dif_refl.xyzw());

				// collider contacts (black/grey)
				for (uint32_t n = 0; n < vf.size(); n++) {
					const t_contact_vertex& cv = vf[n];
					const t_pos4f& p0 = cv.get_point();

					glPushMatrix();
					glTranslatef(p0.x(), p0.y(), p0.z());
					glScalef(0.04f, 0.04f, 0.04f);
					gl.DrawElements(GL_TRIANGLES, render_meshes[1].m_inds.size(), &render_meshes[1].m_verts[0], &render_meshes[1].m_inds[0]);
					glPopMatrix();
				}
			}
			if (!vs.empty()) {
				// glMaterialfv(GL_FRONT, GL_AMBIENT, mesh_materials[6].amb_refl.xyzw());
				// glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_materials[6].dif_refl.xyzw());

				// collidee contacts (white/grey, same positions as collider)
				for (uint32_t n = 0; n < vs.size(); n++) {
					const t_contact_vertex& cv = vs[n];
					const t_pos4f& p0 = cv.get_point();

					glPushMatrix();
					glTranslatef(p0.x(), p0.y(), p0.z());
					glScalef(0.04f, 0.04f, 0.04f);
					gl.DrawElements(GL_TRIANGLES, render_meshes[1].m_inds.size(), &render_meshes[1].m_verts[0], &render_meshes[1].m_inds[0]);
					glPopMatrix();
				}
			}
		}
	}

	gl.DisableClientArrays();

	// flush the rendering pipeline
	// glFlush();
	SDL_GL_SwapBuffers();
}


void t_program_state::toggle_fscr() {
	SDL_Surface* draw_surf = nullptr;

	if ((draw_surf = SDL_GetVideoSurface()) == nullptr || (SDL_WM_ToggleFullScreen(draw_surf) != 1)) {
		printf("[pr::%s] unable to toggle fullscreen: %s\n", __func__, SDL_GetError());
	}
}

bool t_program_state::create_window(uint32_t bpp, uint32_t flags) {
	SDL_Surface* draw_surf = nullptr;

	if ((draw_surf = SDL_SetVideoMode(m_window_size_x, m_window_size_y, bpp, flags)) == nullptr)
		return false;

	reshape_viewport(m_window_size_x, m_window_size_y);
	return true;
}

void t_program_state::reshape_viewport(uint32_t width, uint32_t height) {
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0f, (width * 1.0f) / height, 0.25f, 1024.0f);
	// glFrustum(-1.0f, 1.0f,  -1.0f, 1.0f,  znear, zfar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


void t_program_state::loop() {
	t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);

	while (m_is_running) {
		SDL_Event event;

		// grab first event from queue and process it
		if (SDL_PollEvent(&event)) {
			switch (event.type) {
				case SDL_QUIT: {
					m_is_running = false;
				} break;

				case SDL_VIDEORESIZE: {
					reshape_viewport(event.resize.w, event.resize.h);
				} break;

				case SDL_ACTIVEEVENT: {
					// check if active-state changed (i.e. if app was iconified)
					if (event.active.state & SDL_APPACTIVE) {
						if (event.active.gain) {
							m_is_visible = true;
						} else {
							m_is_visible = false;
						}
					}

					// check if mouse-cursor left/entered window
					if (event.active.state & SDL_APPMOUSEFOCUS) {
						if (event.active.gain) {
							m_has_mouse_focus = true;
						} else {
							m_has_mouse_focus = false;
						}
					}

					// check if window gained/lost input focus
					if (event.active.state & SDL_APPINPUTFOCUS) {
						if (event.active.gain) {
							m_has_keyboard_focus = true;
						} else {
							m_has_keyboard_focus = false;
						}
					}
				} break;

				case SDL_KEYDOWN: {
					// struct SDL_KeyboardEvent { ... struct SDL_keysym { } ... }
					assert(event.key.state == SDL_PRESSED);
					handle_key_press(event.key.keysym);
				} break;
				case SDL_KEYUP: {
					assert(event.key.state == SDL_RELEASED);
					handle_key_release(event.key.keysym);
				} break;
			}
		} else {
			// no events to poll
			if (!m_is_visible) {
				SDL_WaitEvent(nullptr);
			} else {
				update();
				render();
			}
		}
	}
}

