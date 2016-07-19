#ifndef NEWTONICS_HIGHRES_CLOCK_HDR
#define NEWTONICS_HIGHRES_CLOCK_HDR

// whether to use {boost,std}::chrono or native timers
// #define FORCE_CHRONO_TIMERS
// whether to use timers from boost::chrono or std::chrono
#define FORCE_BOOST_CHRONO

// needed for QPC which wants qword-aligned LARGE_INTEGER's
#define __FORCE_ALIGN_STACK__ __attribute__ ((force_align_arg_pointer))

// determines whether to use TGT or QPC on win32
#define USE_NATIVE_HIGHRES_MODE


#if (defined(FORCE_CHRONO_TIMERS))
	// only works if __cplusplus is defined properly by compiler
	// #if ((__cplusplus > 199711L) && !defined(FORCE_BOOST_CHRONO))
	#if (!defined(FORCE_BOOST_CHRONO))
		#include <chrono>
		namespace chrono { using namespace std::chrono; };
	#else
		#include <boost/chrono/include.hpp>
		namespace chrono { using namespace boost::chrono; };
	#endif
#else
	#if (defined(WIN32))
		#include <windows.h>
	#endif
#endif

#include <cstdint>
#include <ctime>



#if (defined(FORCE_CHRONO_TIMERS))
	struct t_chrono_clock {
	public:
		static void push_tick_rate() {}
		static void pop_tick_rate() {}

		static int64_t get_ticks() {
			// high_res_clock may just be an alias of system_clock
			// (so a HIGHRES_MODE conditional is not of value here)
			//
			// const chrono::system_clock::time_point cur_time = chrono::system_clock::now();
			const chrono::high_resolution_clock::time_point cur_time = chrono::high_resolution_clock::now();
			const chrono::nanoseconds run_time = chrono::duration_cast<chrono::nanoseconds>(cur_time.time_since_epoch());

			// number of ticks since chrono's epoch (note that there exist
			// very strong differences in behavior between the boost:: and
			// std:: versions with gcc 4.6.1)
			return (run_time.count());
		}

		static const char* get_name() {
			#if (defined(FORCE_BOOST_CHRONO))
			return "boost::chrono::high_resolution_clock";
			#else
			return   "std::chrono::high_resolution_clock";
			#endif
		}
	};

	typedef t_chrono_clock t_hres_clock;


#else


	#if (defined(WIN32))
		#if (!defined(FORCE_CHRONO_TIMERS))
		struct t_native_windows_clock {
		public:
			static void push_tick_rate() {
				#if (!defined(USE_NATIVE_HIGHRES_MODE))
				// set the number of milliseconds between interrupts
				// NOTE: this is a *GLOBAL* setting, not per-process
				timeBeginPeriod(1);
				#endif
			}
			static void pop_tick_rate() {
				#if (!defined(USE_NATIVE_HIGHRES_MODE))
				timeEndPeriod(1);
				#endif
			}

			__FORCE_ALIGN_STACK__
			static int64_t get_ticks() {
				#if (!defined(USE_NATIVE_HIGHRES_MODE))
				// timeGetTime is affected by time{Begin,End}Period whereas
				// GetTickCount is not ---> resolution of the former can be
				// configured but not for a specific process (they both read
				// from a shared counter that is updated by the system timer
				// interrupt)
				// it returns "the time elapsed since Windows was started"
				// (usually not a very large value so there is little risk
				// of overflowing)
				//
				// note: there is a GetTickCount64 but no timeGetTime64
				//
				return (nsecs_from_msecs<boost::uint32_t>(timeGetTime()));
				#else
				// NOTE:
				//   SDL 1.2 by default does not use QueryPerformanceCounter
				//   SDL 2.0 does, but code does not seem aware of its issues
				//
				//   QPC is an interrupt-independent (unlike timeGetTime & co)
				//   *virtual* timer that runs at a "fixed" frequency which is
				//   derived from hardware, but can be *severely* affected by
				//   thermal drift (heavy CPU load will change the precision)
				//
				//   more accurately QPC is an *interface* to either the TSC
				//   or the HPET or the ACPI timer, MS claims "it should not
				//   matter which processor is called" and setting the thread
				//   affinity is only necessary in case QPC picks TSC (which
				//   can happen if ACPI BIOS code is broken)
				//
				//      const DWORD_PTR oldMask = SetThreadAffinityMask(::GetCurrentThread(), 0);
				//      QueryPerformanceCounter(...);
				//      SetThreadAffinityMask(::GetCurrentThread(), oldMask);
				//
				//   TSC is not invariant and completely unreliable on multi-proc
				//   systems, but there exists an enhanced TSC on modern hardware
				//   which IS invariant (check CPUID 80000007H:EDX[8]) --> useful
				//   because reading TSC is much faster than an API call like QPC
				//   (note: it might work on single-proc multi-core systems since
				//   these are not SMP, cf. the POSIX clock_gettime docs)
				//
				//   the range of possible frequencies is extreme (KHz - GHz) and
				//   the hardware counter might only have a 32-bit register while
				//   QuadPart is a 64-bit integer --> no monotonicity guarantees
				//   (especially in combination with TSC if thread switches procs)
				//
				LARGE_INTEGER tick_freq;
				LARGE_INTEGER curr_tick;

				if (!QueryPerformanceFrequency(&tick_freq))
					return (nsecs_from_msecs<int64_t>(0));

				QueryPerformanceCounter(&curr_tick);

				// we want the raw tick (uncorrected for frequency)
				//
				// if clock ticks <freq> times per second, then the
				// total number of {milli,micro,nano}seconds elapsed
				// for any given tick is <tick> / <freq / resolution>
				// e.g. if freq = 15000Hz and tick = 5000, then
				//
				//        secs = 5000 / (15000 / 1e0) =                    0.3333333
				//   millisecs = 5000 / (15000 / 1e3) = 5000 / 15.000000 =       333
				//   microsecs = 5000 / (15000 / 1e6) = 5000 /  0.015000 =    333333
				//    nanosecs = 5000 / (15000 / 1e9) = 5000 /  0.000015 = 333333333
				//
				if (tick_freq.QuadPart >= int64_t(1e9)) return (nsecs_from_nsecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-9))));
				if (tick_freq.QuadPart >= int64_t(1e6)) return (nsecs_from_usecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-6))));
				if (tick_freq.QuadPart >= int64_t(1e3)) return (nsecs_from_msecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-3))));

				return (from_ssecs<int64_t>(std::max(0LL, curr_tick.QuadPart)));
				#endif
			}

			static const char* get_name() {
				#if (!defined(USE_NATIVE_HIGHRES_MODE))
				return "native::win32::TimeGetTime";
				#else
				return "native::win32::QueryPerformanceCounter";
				#endif
			}

		private:
			template<typename T>  static  int64_t nsecs_from_ssecs(const T  s) { return ( s * int64_t(1e9)); }
			template<typename T>  static  int64_t nsecs_from_msecs(const T ms) { return (ms * int64_t(1e6)); }
			template<typename T>  static  int64_t nsecs_from_usecs(const T us) { return (us * int64_t(1e3)); }
			template<typename T>  static  int64_t nsecs_from_nsecs(const T ns) { return (ns               ); }
		};

		typedef t_native_windows_clock t_native_clock;

		#endif

	#else

		struct t_native_posix_clock {
		public:
			static void push_tick_rate() {}
			static void pop_tick_rate() {}

			// CLOCK_{REALTIME,MONOTONIC,PROCESS_CPUTIME_ID,THREAD_CPUTIME_ID}
			static int64_t get_ticks(clockid_t clock_id = CLOCK_REALTIME) {
				return (get_time(clock_id));
			}

			static const char* get_name() {
				return "native::posix::gettime";
			}

		private:
			static int64_t get_time(clockid_t clock_id) {
				timespec ts;
				// number of nanoseconds since (unknown) epoch
				if (clock_gettime(clock_id, &ts) == 0)
					return (ts.tv_sec * 1000000000L + ts.tv_nsec);
				return -1;
			}

			static int64_t get_resolution(clockid_t clock_id) {
				timespec ts;
				// 1 := ns, 1000 := us, etc
				if (clock_getres(clock_id, &ts) == 0)
					return ts.tv_nsec;
				return -1;
			}
		};

		typedef t_native_posix_clock t_native_clock;

	#endif

	typedef t_native_clock t_hres_clock;

#endif






/*
struct t_hres_clock {
public:
	// NOTE:
	//   1e-x are double-precision literals but T can be float
	//   floats only provide ~6 decimal digits of precision so
	//   to_ssecs is inaccurate in that case
	template<typename T>  static  T to_ssecs(const int64_t ns) { return (ns * 1e-9); }
	template<typename T>  static  T to_msecs(const int64_t ns) { return (ns * 1e-6); }
	template<typename T>  static  T to_usecs(const int64_t ns) { return (ns * 1e-3); }
	template<typename T>  static  T to_nsecs(const int64_t ns) { return (ns       ); }

	// specializations
	// template<>  static  int64_t to_ssecs<int64_t>(const int64_t ns) { return (ns / int64_t(1e9)); }
	// ...
	static int64_t to_ssecs(const int64_t ns) { return (ns / int64_t(1e9)); }
	static int64_t to_msecs(const int64_t ns) { return (ns / int64_t(1e6)); }
	static int64_t to_usecs(const int64_t ns) { return (ns / int64_t(1e3)); }

	// these convert inputs to nanoseconds
	template<typename T>  static  int64_t nsecs_from_ssecs(const T  s) { return ( s * int64_t(1e9)); }
	template<typename T>  static  int64_t nsecs_from_msecs(const T ms) { return (ms * int64_t(1e6)); }
	template<typename T>  static  int64_t nsecs_from_usecs(const T us) { return (us * int64_t(1e3)); }
	template<typename T>  static  int64_t nsecs_from_nsecs(const T ns) { return (ns               ); }


	static void push_tick_rate() {
		#if USE_NATIVE_WINDOWS_CLOCK
			#if (USE_NATIVE_HIGHRES_MODE == 0)
			timeBeginPeriod(1);
			#endif
		#endif
	}
	static void pop_tick_rate() {
		#if USE_NATIVE_WINDOWS_CLOCK
			#if (USE_NATIVE_HIGHRES_MODE == 0)
			timeEndPeriod(1);
			#endif
		#endif
	}



	#if USE_NATIVE_WINDOWS_CLOCK
	__FORCE_ALIGN_STACK__
	static int64_t get_ticks_windows() {
		#if (defined(USE_NATIVE_HIGHRES_MODE))
			LARGE_INTEGER tick_freq;
			LARGE_INTEGER curr_tick;

			if (!QueryPerformanceFrequency(&tick_freq))
				return (nsecs_from_msecs<int64_t>(0));

			QueryPerformanceCounter(&curr_tick);

			if (tick_freq.QuadPart >= int64_t(1e9)) return (nsecs_from_nsecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-9))));
			if (tick_freq.QuadPart >= int64_t(1e6)) return (nsecs_from_usecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-6))));
			if (tick_freq.QuadPart >= int64_t(1e3)) return (nsecs_from_msecs<boost::uint64_t>(std::max(0.0, curr_tick.QuadPart / (tick_freq.QuadPart * 1e-3))));

			return (from_ssecs<int64_t>(std::max(0LL, curr_tick.QuadPart)));
		#else
			return (nsecs_from_msecs<boost::uint32_t>(timeGetTime()));
		#endif
	}
	#endif

	static int64_t get_ticks() {
		#if USE_NATIVE_WINDOWS_CLOCK
		return (get_ticks_windows());
		#else
		const chrono::high_resolution_clock::time_point cur_time = chrono::high_resolution_clock::now();
		const chrono::nanoseconds run_time = chrono::duration_cast<chrono::nanoseconds>(cur_time.time_since_epoch());

		return (run_time.count());
		#endif
	}

	static const char* get_name() {
		#if USE_NATIVE_WINDOWS_CLOCK
			#if (defined(USE_NATIVE_HIGHRES_MODE))
				return "win32::QueryPerformanceCounter";
			#else
				return "win32::TimeGetTime";
			#endif
		#else
			#ifdef TIME_USING_LIBCHRONO
			return "boost::chrono::high_resolution_clock";
			#endif
			#ifdef TIME_USING_STDCHRONO
			return "std::chrono::high_resolution_clock";
			#endif
		#endif
	}
};
*/


/*
struct t_clock_tick {
public:
	t_clock_tick(): m_tick(0) {}

	// common-case constructor
	template<typename T> explicit t_clock_tick(const T msecs) {
		set_tick(t_hres_clock::nsecs_from_msecs(msecs));
	}

	void set_tick(int64_t t) { m_tick = t; }
	int64_t get_tick() const { return m_tick; }

	t_clock_tick& operator += (const t_clock_tick ct)       { m_tick += ct.get_tick(); return *this; }
	t_clock_tick& operator -= (const t_clock_tick ct)       { m_tick -= ct.get_tick(); return *this; }
	t_clock_tick& operator %= (const t_clock_tick ct)       { m_tick %= ct.get_tick(); return *this;    }
	t_clock_tick  operator -  (const t_clock_tick ct) const { return (get_clock_tick(m_tick - ct.get_tick())); }
	t_clock_tick  operator +  (const t_clock_tick ct) const { return (get_clock_tick(m_tick + ct.get_tick())); }
	t_clock_tick  operator %  (const t_clock_tick ct) const { return (get_clock_tick(m_tick % ct.get_tick())); }

	bool operator <  (const t_clock_tick ct) const { return (m_tick <  ct.get_tick()); }
	bool operator >  (const t_clock_tick ct) const { return (m_tick >  ct.get_tick()); }
	bool operator <= (const t_clock_tick ct) const { return (m_tick <= ct.get_tick()); }
	bool operator >= (const t_clock_tick ct) const { return (m_tick >= ct.get_tick()); }


	// short-hands for to_*secs_t
	int64_t to_ssecs_i() const { return (to_ssecs_t<int64_t>()); }
	int64_t to_msecs_i() const { return (to_msecs_t<int64_t>()); }
	int64_t to_usecs_i() const { return (to_usecs_t<int64_t>()); }
	int64_t to_nsecs_i() const { return (to_nsecs_t<int64_t>()); }

	float to_ssecs_f() const { return (to_ssecs_t<float>()); }
	float to_msecs_f() const { return (to_msecs_t<float>()); }
	float to_usecs_f() const { return (to_usecs_t<float>()); }
	float to_nsecs_f() const { return (to_nsecs_t<float>()); }

	template<typename T> T to_ssecs_t() const { return (t_hres_clock::to_ssecs<T>(m_tick)); }
	template<typename T> T to_msecs_t() const { return (t_hres_clock::to_msecs<T>(m_tick)); }
	template<typename T> T to_usecs_t() const { return (t_hres_clock::to_usecs<T>(m_tick)); }
	template<typename T> T to_nsecs_t() const { return (t_hres_clock::to_nsecs<T>(m_tick)); }


	static t_clock_tick get_curr_time(bool init_call = false) {
		assert(get_epoch_time() != 0 || init_call);
		return (get_clock_tick(t_hres_clock::get_ticks()));
	}
	static t_clock_tick get_init_time() {
		assert(get_epoch_time() != 0);
		return (get_clock_tick(get_epoch_time()));
	}
	static t_clock_tick get_diff_time() {
		return (get_curr_time() - get_init_time());
	}

	static void set_epoch_time(const t_clock_tick ct) {
		assert(get_epoch_time() == 0);
		get_epoch_time() = ct.get_tick();
		assert(get_epoch_time() != 0);
	}

	// initial time (arbitrary epoch, e.g. program start)
	// all other time-points will be larger than this if
	// the clock is monotonically increasing
	static int64_t& get_epoch_time() {
		static int64_t epoch_time = 0;
		return epoch_time;
	}

	static t_clock_tick time_from_nsecs(const int64_t ns) { return (get_clock_tick(t_hres_clock::nsecs_from_nsecs(ns))); }
	static t_clock_tick time_from_usecs(const int64_t us) { return (get_clock_tick(t_hres_clock::nsecs_from_usecs(us))); }
	static t_clock_tick time_from_msecs(const int64_t ms) { return (get_clock_tick(t_hres_clock::nsecs_from_msecs(ms))); }
	static t_clock_tick time_from_ssecs(const int64_t  s) { return (get_clock_tick(t_hres_clock::nsecs_from_ssecs( s))); }

private:
	// convert integer to t_clock_tick object (n is interpreted as number of nanoseconds)
	static t_clock_tick get_clock_tick(const int64_t n) { t_clock_tick ct; ct.set_tick(n); return ct; }

private:
	int64_t m_tick;
};
*/

#endif

