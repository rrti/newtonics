#ifndef NEWTONICS_PROGRAM_STATE_HDR
#define NEWTONICS_PROGRAM_STATE_HDR

#include <cstdint>
#include <array>
#include <vector>

#include "global_defs.hpp"

struct SDL_Surface;
struct SDL_keysym;

struct t_arg_parser;
struct t_world_state;
struct t_program_clock;
struct t_camera;

struct t_key_state {
	uint32_t symb; // e.g. CAM_ACTION_*
	uint32_t flag; // state (down, up, ...)
	uint32_t xxxx;
	float sign;
};

namespace newtonics {
	namespace lib_physics {
		struct t_rigid_body;
	};
};


struct t_program_state {
public:
	t_program_state(const t_arg_parser* arg_parser);

	void init(t_world_state* wstate, t_program_clock* pclock);
	void loop();
	void kill();

	bool tick_clock(t_program_clock* clock);
	void print_fps() {
		#if (ENABLE_PROGSTATE_DEBUG_PRINTS == 123)
		// achtung, possibly not the main thread
		printf("[pr::%s][fps=%u]\n", __func__, m_num_draw_calls);
		#endif
	}

	bool is_running() const { return m_is_running; }
	bool is_paused() const { return m_is_paused; }
	bool is_inited() const { return m_is_inited; }

private:
	void update();
	void update_tracked_object(const std::vector<newtonics::lib_physics::t_rigid_body>& objects, int32_t dir);
	void update_scene(uint32_t scene_offset);
	void update_camera(t_camera& cam, float dif_time_ssecs);
	void render();

	void update_keys(const uint8_t* sdl_key_state);
	void handle_key_press(const SDL_keysym& ks);
	void handle_key_release(const SDL_keysym& ks);

	bool create_window(uint32_t bpp, uint32_t flags);
	void reshape_viewport(uint32_t width, uint32_t height);

	void toggle_fscr();

private:
	t_program_clock* m_prog_clock;
	t_world_state* m_world_state;

	// indexed by SDLK_* code
	std::array<t_key_state, 512> m_key_states;
	// indexed by *_ACTION_* code
	std::array<t_key_state*, 512> m_key_rbinds;

	// current window dimensions
	uint32_t m_window_size_x;
	uint32_t m_window_size_y;

	uint32_t m_num_sim_frames; // number of physics simulation steps so far
	uint32_t m_max_sim_frames; // number of physics simulation steps to run
	uint32_t m_num_draw_calls; // number of update() (i.e. draw()) calls made

	#if (ENABLE_TIMER_THREAD == 1)
	typedef std::atomic<bool> t_prog_state_bool;
	#else
	typedef bool t_prog_state_bool;
	#endif

	t_prog_state_bool m_is_running;
	t_prog_state_bool m_is_paused;
	t_prog_state_bool m_is_visible;
	t_prog_state_bool m_is_inited;

	bool m_has_mouse_focus;
	bool m_has_keyboard_focus;
};

#endif

