#include <algorithm>

#include "program_clock.hpp"
#include "misc/highres_clock.hpp"

bool t_program_clock::tick(bool is_paused) {
	const uint64_t cur_tick = GET_CURRENT_TICK();

	m_cur_tick_count  =   cur_tick;
	m_dif_tick_count  = m_cur_tick_count - m_prv_tick_count; // in ns
	m_prv_tick_count  =   cur_tick;
	m_acc_tick_count += (m_dif_tick_count * (1 - is_paused));
	m_tot_tick_count += (m_dif_tick_count                  );

	// calculate the time extrapolation factor based on
	// how much delta-time we have accumulated since the
	// last simulation step
	m_tick_extr_fact = std::min(1.0f, (m_acc_tick_count * 1.0f) / get_time_step_size());
	m_tot_tick_count = m_tot_tick_count * (m_tot_tick_count < (1000lu * 1000lu * 1000lu));

	// when true, one wall-clock second has elapsed
	return (m_tot_tick_count == 0);
}

