#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "global_defs.hpp"
#include "program_clock.hpp"
#include "program_state.hpp"
#include "world_state.hpp"
#include "misc/arg_parser.hpp"
#include "misc/highres_clock.hpp"
#include "misc/scoped_timer.hpp"

#if (ENABLE_TIMER_THREAD == 1)
// responsible for tick()'ing the clock
//
// note: this runs at full tilt, which means a float can not
// accurately represent the deltas when converted to seconds
static void clock_thread_func(t_program_state* p, t_program_clock* c) {
	while (!p->is_inited()) {
	}

	while (p->is_running()) {
		if (p->tick_clock(c))
			print_fps();

		#if 0
		boost::this_thread::sleep(boost::posix_time::microseconds(100));
		#endif
	}
}
#endif

int main(int argc, const char** argv) {
	// need to call ptr when ENABLE_HIGHRES_CLOCK is 1
	t_hres_clock::push_tick_rate();
	t_scoped_timer::init_timing_data();

	printf("[%s][clock=%s][ticks=%ld]\n", __func__, t_hres_clock::get_name(), t_hres_clock::get_ticks());
	srandom(t_hres_clock::get_ticks());

	{
		t_world_state::init_global_scenes();
		t_world_state::exec_global_tests();
		t_world_state world_state;

		t_scoped_timer scoped_timer(__PRETTY_FUNCTION__);
		t_arg_parser arg_parser;

		arg_parser.add_arg_obj(t_arg_obj("num_objects",    "0"));
		arg_parser.add_arg_obj(t_arg_obj("max_sframes", "1000"));
		arg_parser.add_arg_obj(t_arg_obj("scene_index",    "0"));
		arg_parser.add_arg_obj(t_arg_obj("win_size_x",   "800"));
		arg_parser.add_arg_obj(t_arg_obj("win_size_y",   "600"));
		arg_parser.add_arg_alias("num_objects", "--no");
		arg_parser.add_arg_alias("max_sframes", "--mf");
		arg_parser.add_arg_alias("scene_index", "--si");
		arg_parser.add_arg_alias("win_size_x", "--wsx");
		arg_parser.add_arg_alias("win_size_y", "--wsy");

		arg_parser.parse_args(argc, argv);
		world_state.load_scene(arg_parser.get_int_arg("scene_index"), arg_parser.get_int_arg("num_objects"));

		#ifndef RUN_HEADLESS
		t_program_state program_state(&arg_parser);
		t_program_clock program_clock;

		#if (ENABLE_TIMER_THREAD == 1)
		boost::thread clock_thread(boost::bind(&clock_thread_func, &program_state, &program_clock));
		#endif

		program_state.init(&world_state, &program_clock);
		program_state.loop();

		#if (ENABLE_TIMER_THREAD == 1)
		clock_thread.join();
		#endif

		#else

		world_state.run_headless(&arg_parser);
		#endif
	}

	t_scoped_timer::print_timing_data();
	t_hres_clock::pop_tick_rate();

	return 0;
}

