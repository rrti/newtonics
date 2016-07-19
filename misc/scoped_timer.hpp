#ifndef NEWTONICS_SCOPED_TIMER_HDR
#define NEWTONICS_SCOPED_TIMER_HDR

#include <cstdint>

struct t_scoped_timer {
public:
	struct t_timing_data {
	public:
		t_timing_data(): stime(0), ctime(0), calls(0), depth(0) {}
		t_timing_data(const t_timing_data& t): stime(t.stime), ctime(t.ctime), calls(t.calls), depth(t.depth) {}

		int64_t stime;
		int64_t ctime;
		int64_t calls;
		int64_t depth;
	};

	t_scoped_timer(const char* scope_name);
	~t_scoped_timer();

	const t_scoped_timer* get_parent_timer() const { return m_parent_timer; }
	const char* get_scope_name() const { return m_scope_name; }

	static t_timing_data get_timing_data(const char* scope_name);
	static void  init_timing_data();
	static void print_timing_data();

private:
	// allows walking the timer-stack more easily
	const t_scoped_timer* m_parent_timer;

	// name of the function-scope this timer was created in
	const char* m_scope_name;

	int64_t m_ctor_tick;
	int64_t m_dtor_tick;
};

#endif

