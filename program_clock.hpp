#ifndef NEWTONICS_PROGRAM_CLOCK_HDR
#define NEWTONICS_PROGRAM_CLOCK_HDR

#include "global_defs.hpp"

#include <cstdint>

#if (ENABLE_TIMER_THREAD == 1)
#include <atomic>
#endif

struct t_program_clock {
public:
	// initialized later (depends on SDL)
	t_program_clock() { reset(0); }

	void reset(uint64_t cur_tick) {
		m_cur_tick_count = cur_tick;
		m_prv_tick_count = cur_tick;
		m_dif_tick_count = 0;
		m_acc_tick_count = 0;
		m_tot_tick_count = 0;

		m_tick_rate_mult = 1.0f;
		m_tick_extr_fact = 0.0f;
		m_time_step_sign = 1.0f;
	}

	bool tick(bool is_paused);

	// "a = a * s" because atomic<float> has no operator*=
	void scale_acc_ticks(uint64_t s) { m_acc_tick_count = m_acc_tick_count * s; }
	void scale_tick_rate(float s) { m_tick_rate_mult = m_tick_rate_mult * s; }
	void scale_time_step(float s) { m_time_step_sign = m_time_step_sign * s; }

	void add_acc_ticks(uint64_t dt) { m_acc_tick_count += dt; }
	void sub_acc_ticks(uint64_t dt) { m_acc_tick_count -= dt; }

	uint64_t get_time_step_size() const { return (TIME_STEP_SIZE_NS / m_tick_rate_mult); } // in ns
	uint64_t get_num_time_steps() const { return (m_acc_tick_count / get_time_step_size()); }

	float get_tick_rate_mult() const { return m_tick_rate_mult; }
	float get_tick_extr_fact() const { return m_tick_extr_fact; }
	float get_time_step_sign() const { return m_time_step_sign; }

	float get_dif_time_msecs() const { return (  m_dif_tick_count   * 1e-6f); } // ns to ms
	float get_dif_time_ssecs() const { return (get_dif_time_msecs() * 1e-1f); } // ms to s/100

private:
	#if (ENABLE_TIMER_THREAD == 1)
	typedef std::atomic<uint64_t> t_tick_count;
	typedef std::atomic<   float> t_tick_delta;
	#else
	typedef uint64_t t_tick_count;
	typedef    float t_tick_delta;
	#endif

	t_tick_count m_cur_tick_count; // tick at start of current update_state()
	t_tick_count m_prv_tick_count; // tick at end of previous update_state()
	t_tick_count m_dif_tick_count; // delta between cur_* and prv_* ticks
	t_tick_count m_acc_tick_count; // running sum of dif_* tick values (constant when paused)
	t_tick_count m_tot_tick_count; // running sum of dif_* tick values (cleared every second)

	t_tick_delta m_tick_rate_mult; // multiplier for the number of sim-steps per second
	t_tick_delta m_tick_extr_fact; // extrapolation factor between last and next sim-step
	t_tick_delta m_time_step_sign; // +1 or -1
};

#endif

