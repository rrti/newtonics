#ifndef NEWTONICS_DEFINES_HDR
#define NEWTONICS_DEFINES_HDR

#if 1
#define printf(...) printf(__VA_ARGS__)
#else
#define printf(...)
#endif

#define MAKE_STRING_CONST(s) #s

#define ENABLE_FULLSCREEN_MODE 0
#define ENABLE_HIGHRES_CLOCK 1
#define ENABLE_SCOPED_TIMERS 1
#define ENABLE_TIMER_THREAD 0

#define ENABLE_PROGSTATE_DEBUG_PRINTS 1
#define ENABLE_PHYSSTATE_DEBUG_PRINTS 1
#define ENABLE_RIGIDBODY_DEBUG_PRINTS 1
#define ENABLE_PLANEBODY_DEBUG_PRINTS 0
#define ENABLE_COLLIMESH_DEBUG_PRINTS 0

#if (ENABLE_HIGHRES_CLOCK == 1)
#define GET_CURRENT_TICK() t_hres_clock::get_ticks()
#else
// SDL's clock only has millisecond-resolution
#define GET_CURRENT_TICK() (SDL_GetTicks() * 1e6f)
#endif

// N millisecs per step, so <1000/N> steps per second
// each step advances clock-time by <N * 1e-3> seconds
#define TIME_STEP_SIZE_MS 10
#define TIME_STEP_SIZE_NS (TIME_STEP_SIZE_MS *  1e6f)
#define TIME_STEP_SIZE_S  (TIME_STEP_SIZE_MS * 1e-3f)

// if 0, integrator reduces to plain Euler
#define ENABLE_RK4_INTEGRATION 1

#define ENABLE_REACTION_FORCES 1
#define ENABLE_COLLISION_FORCES 1
#define ENABLE_FRICTION_FORCES 0
#define ENABLE_ANGULAR_FORCES 1

#define CLAMP_FORCE_RESPONSES 1
#define RANDOMIZE_OBJECT_ORDER 0


#define MAX_OBJECT_COUNT 20u
#define MAX_OBJECT_CONTACT_AGE (0 * ((1000u / TIME_STEP_SIZE_MS) / 5u))
// keep old contacts around until new ones are calculated
// #define MAX_OBJECT_CONTACT_AGE -1u

#define MAX_FORCE_EXCHANGE_ITERS 40u

#endif

