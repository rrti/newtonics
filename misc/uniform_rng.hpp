#ifndef NEWTONICS_UNIFORM_RNG_HDR
#define NEWTONICS_UNIFORM_RNG_HDR

#include <chrono>
#include <random>

template<typename t_real = double> struct t_uniform_real_rng {
public:
	typedef  std::mt19937  m_sampler_type;
	// typedef  std::default_random_engine  m_sampler_type;

	typedef  typename std::uniform_real_distribution<t_real>  m_distrib_type;
	typedef  typename m_distrib_type::param_type  m_distrib_param_type;

	t_uniform_real_rng(uint64_t seed = 0) { set_seed(seed); }

	void set_seed(uint64_t seed) {
		if (seed == 0) {
			const std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
			const std::chrono::nanoseconds run_time = std::chrono::duration_cast<std::chrono::nanoseconds>(cur_time.time_since_epoch());

			m_sampler.seed(run_time.count());
		} else {
			m_sampler.seed(seed);
		}
	}

	void set_range(t_real min, t_real max) {
		m_distrib_param_type params = m_distrib.param();

		params.min(min);
		params.max(max);

		m_distrib.reset();
		m_distrib.param(params);
	}

	// same as std::generate, but with implicit generator
	template<typename t_iter> void gen_samples(t_iter first, t_iter last) {
		while (first != last) {
			(*first)++ = (*this)();
		}
	}

	t_real operator ()() { return (m_distrib(m_sampler)); }

private:
	m_sampler_type m_sampler;
	m_distrib_type m_distrib;
};



template<typename t_int = uint64_t> struct t_uniform_int_rng {
public:
	typedef  std::mt19937  m_sampler_type;
	// typedef  std::default_random_engine  m_sampler_type;

	typedef  typename std::uniform_int_distribution<t_int>  m_distrib_type;
	typedef  typename m_distrib_type::param_type  m_distrib_param_type;

	t_uniform_int_rng(uint64_t seed = 0) { set_seed(seed); }

	void set_seed(uint64_t seed) {
		if (seed == 0) {
			const std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
			const std::chrono::nanoseconds run_time = std::chrono::duration_cast<std::chrono::nanoseconds>(cur_time.time_since_epoch());

			m_sampler.seed(run_time.count());
		} else {
			m_sampler.seed(seed);
		}
	}

	void set_range(t_int min, t_int max) {
		m_distrib_param_type params = m_distrib.param();

		params.min(min);
		params.max(max);

		m_distrib.reset();
		m_distrib.param(params);
	}

	template<typename t_iter> void gen_samples(t_iter first, t_iter last) {
		while (first != last) {
			(*first)++ = (*this)();
		}
	}

	t_int operator ()() { return (m_distrib(m_sampler)); }

private:
	m_sampler_type m_sampler;
	m_distrib_type m_distrib;
};

#endif

