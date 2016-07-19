#ifndef NEWTONICS_PHYSICS_SCENE_HDR
#define NEWTONICS_PHYSICS_SCENE_HDR

#include "global_defs.hpp"
#include "math/tuple.hpp"
#include "math/vector.hpp"

struct t_physics_scene {
public:
	t_physics_scene() {
		num_objects = 5;

		object_default_mass = 100.0f;

		object_ignore_lin_forces = false;
		object_ignore_ang_forces = false;

		// object_grid_divs = t_tup4i(1000, 10, 1000;
		object_grid_divs = t_tup4i(  50,  5,   50);
		object_grid_mins = t_vec4f(-500.0f, -500.0f, -500.0f);
		object_grid_maxs = t_vec4f( 500.0f,  500.0f,  500.0f);

		object_force_coeffs = t_tup4f(0.0f, 1.0f, 0.0f, 0.0f);
		object_misc_coeffs = t_tup4f(1.0f, 1.0f, 0.0f, 0.0f);

		// convert from m/s^2 to m/dt^2; quantities in units per frame
		//   10m / (  1s *   1s) = 10m / (   1000ms  *    1000ms  ) =
		//   10m / (100f * 100f) = 10m * ((1 / 100f) * (1 / 100f)))
		//
		// note that with gravity disabled objects can gain
		// random vertical velocity due to numerical drift
		//
		object_gravity_accell = ((t_vec4f(0.0f, 1.0f, 0.0f) * -10.0f) * (TIME_STEP_SIZE_S * TIME_STEP_SIZE_S)); // y_axis_vec
		object_tst_lin_accell = ((t_vec4f(1.0f, 0.0f, 0.0f) *  10.0f) * (TIME_STEP_SIZE_S * TIME_STEP_SIZE_S)); // x_axis_vec
		object_tst_ang_accell = ((t_vec4f(0.0f, 1.0f, 0.0f) *  10.0f) * (TIME_STEP_SIZE_S * TIME_STEP_SIZE_S));

		object_mesh_scales = t_vec4f(1.0f, 1.0f, 1.0f, object_default_mass);
		ground_mesh_scales = t_vec4f(50.0f, 1.0f, 50.0f, 1.0f);

		ground_offset_vector = t_vec4f(0.0f, -2.0f, 0.0f, 1.0f);
	}

public:
	unsigned int num_objects;

	float object_default_mass;

	bool object_ignore_lin_forces;
	bool object_ignore_ang_forces;

	t_vec4f object_axis_angles;
	t_vec4f ground_axis_angles;

	t_tup4i object_grid_divs;
	t_vec4f object_grid_mins;
	t_vec4f object_grid_maxs;

	// restitution_coeff={0: fully inelastic, 1: fully elastic}
	//   with c_r=0.0, every object in a connected chain will end up moving from force-propagation
	//   with c_r=1.0, only last object in a connected chain will end up moving from force-propagation
	//
	// .x := object_collis_restitution_coeff (0.0f)
	// .y := ground_collis_restitution_coeff (1.0f)
	// .z := object_static_friction_coeff (4.0f)
	// .w := ground_static_friction_coeff (5.0f)
	t_tup4f object_force_coeffs;
	// .x := lin_velocity_damping_coeff (1.0f)
	// .y := ang_velocity_damping_coeff (1.0f)
	// .z := air_drag_mult
	// .w := unused
	t_tup4f object_misc_coeffs;

	t_vec4f object_gravity_accell;
	t_vec4f object_tst_lin_accell;
	t_vec4f object_tst_ang_accell;

	t_vec4f object_mesh_scales;
	t_vec4f ground_mesh_scales;
	t_vec4f ground_normal_vector;

	t_vec4f object_spacing_vector;
	t_vec4f object_offset_vector;
	t_vec4f ground_offset_vector;
};

#endif

