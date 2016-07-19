#ifndef NEWTONICS_PHYS_STATE_HDR
#define NEWTONICS_PHYS_STATE_HDR

#include "../global_defs.hpp"
#include "../math/matrix.hpp"
#include "../math/point.hpp"
#include "../math/vector.hpp"

namespace newtonics {
	namespace lib_physics {
		using namespace lib_math;

		struct t_diff_state {
		public:
			t_diff_state(
				const t_vec4f& lin_vel = t_vec4f::zero_vector(),
				const t_vec4f& lin_acc = t_vec4f::zero_vector(),
				const t_vec4f& ang_vel = t_vec4f::zero_vector(),
				const t_vec4f& ang_acc = t_vec4f::zero_vector()
			) {
				set_lin_vel(lin_vel);
				set_lin_acc(lin_acc);
				set_ang_vel(ang_vel);
				set_ang_acc(ang_acc);
			}

			t_diff_state operator + (const t_diff_state& d) const {
				t_diff_state r;
				r.set_lin_vel(m_lin_vel + d.get_lin_vel());
				r.set_lin_acc(m_lin_acc + d.get_lin_acc());
				r.set_ang_vel(m_ang_vel + d.get_ang_vel());
				r.set_ang_acc(m_ang_acc + d.get_ang_acc());
				return r;
			}
			t_diff_state operator * (float dt) const {
				t_diff_state r;
				r.set_lin_vel(m_lin_vel * dt);
				r.set_lin_acc(m_lin_acc * dt);
				r.set_ang_vel(m_ang_vel * dt);
				r.set_ang_acc(m_ang_acc * dt);
				return r;
			}

			const t_vec4f& get_lin_vel() const { return m_lin_vel; }
			const t_vec4f& get_lin_acc() const { return m_lin_acc; }
			const t_vec4f& get_ang_vel() const { return m_ang_vel; }
			const t_vec4f& get_ang_acc() const { return m_ang_acc; }
			void set_lin_vel(const t_vec4f& v) { m_lin_vel = v; }
			void set_lin_acc(const t_vec4f& v) { m_lin_acc = v; }
			void set_ang_vel(const t_vec4f& v) { m_ang_vel = v; }
			void set_ang_acc(const t_vec4f& v) { m_ang_acc = v; }

		private:
			t_vec4f m_lin_vel; // linear velocity (v=dx/dt)
			t_vec4f m_lin_acc; // linear acceleration (a=dv/dt)
			t_vec4f m_ang_vel; // angular velocity
			t_vec4f m_ang_acc; // angular acceleration
		};


		struct t_phys_state {
		public:
			t_phys_state(
				const t_pos4f& cur_pos = t_pos4f::zero_point(),
				const t_vec4f& lin_vel = t_vec4f::zero_vector(),
				const t_vec4f& ang_vel = t_vec4f::zero_vector()
			) {
				set_cur_pos(cur_pos);
				set_lin_vel(lin_vel);
				set_ang_vel(ang_vel);
			}

			t_phys_state operator + (const t_diff_state& d) const {
				t_phys_state r;
				r.set_cur_pos(m_cur_pos + d.get_lin_vel());
				r.set_lin_vel(m_lin_vel + d.get_lin_acc());
				r.set_ang_vel(m_ang_vel + d.get_ang_acc());
				return r;
			}

			// for lerp()
			t_phys_state operator + (const t_phys_state& s) const {
				t_phys_state r;
				r.set_cur_pos(s.get_cur_pos() + m_cur_pos.to_vector());
				r.set_lin_vel(s.get_lin_vel() + m_lin_vel            );
				r.set_ang_vel(s.get_ang_vel() + m_ang_vel            );
				return r;
			}
			// for lerp()
			t_phys_state operator * (float s) const {
				t_phys_state r;
				r.set_cur_pos((m_cur_pos.to_vector() * s).to_point());
				r.set_lin_vel(m_lin_vel * s);
				r.set_ang_vel(m_ang_vel * s);
				return r;
			}


			// calculate the time-derivative of state <s>
			t_diff_state differentiate(const t_phys_state& s) const {
				t_diff_state d;
				d.set_lin_vel(s.get_lin_vel());
				d.set_lin_acc(calc_lin_acc(s));
				d.set_ang_acc(calc_ang_acc(s));
				return d;
			}

			// calculate next state at time <t+dt> from
			// current (using four sampled derivatives)
			t_phys_state integrate(float dt) const {
				const t_phys_state s = *this;

				#if (ENABLE_RK4_INTEGRATION == 1)
				// const t_diff_state a;
				const t_diff_state b = differentiate(s /*+ (a * dt * 0.0f)*/);
				const t_diff_state c = differentiate(s   + (b * dt * 0.5f)  );
				const t_diff_state d = differentiate(s   + (c * dt * 0.5f)  );
				const t_diff_state e = differentiate(s   + (d * dt * 1.0f)  );

				constexpr float sum_wgt = 1.0f / 6.0f;
				constexpr float mid_wgt = 2.0f;

				const t_diff_state f = (b + (c + d) * mid_wgt + e) * sum_wgt;
				#else
				const t_diff_state f = differentiate(s);
				#endif

				return (s + f * dt);
			}


			// TODO: actual direction vectors, rotation axis, mass, moi
			t_mat44f calc_transform() const {
				const t_vec4f x = t_vec4f::x_axis_vector();
				const t_vec4f y = t_vec4f::y_axis_vector();
				const t_vec4f z = t_vec4f::z_axis_vector();
				const t_vec4f t = m_cur_pos.to_vector(); // mat::set_t_vec enforces .w=1
				return (t_mat44f(x, y, z, t));
			}

			t_vec4f calc_lin_acc(const t_phys_state& s) const { (void) s; return (calc_spring_force(s)); }
			t_vec4f calc_ang_acc(const t_phys_state& s) const { (void) s; return (calc_spring_force(s)); }

			t_vec4f calc_spring_force(const t_phys_state& s) const {
				constexpr float k = 20.0f; // stiffness
				constexpr float b =  1.0f; // mass
				return (((s.get_cur_pos()).to_vector() * -k) - (s.get_lin_vel() * b));
			}

			const t_pos4f& get_cur_pos() const { return m_cur_pos; }
			const t_vec4f& get_lin_vel() const { return m_lin_vel; }
			const t_vec4f& get_ang_vel() const { return m_ang_vel; }
			void set_cur_pos(const t_pos4f& v) { m_cur_pos = v; }
			void set_lin_vel(const t_vec4f& v) { m_lin_vel = v; }
			void set_ang_vel(const t_vec4f& v) { m_ang_vel = v; }

			// linearly interpolate between states
			static t_phys_state lerp(const t_phys_state& prev_state, const t_phys_state& curr_state, float alpha) {
				const t_phys_state prev = prev_state * (1.0f - alpha);
				const t_phys_state curr = curr_state * (       alpha);
				return (prev + curr);
			}

		private:
			t_pos4f m_cur_pos; // position (x)
			t_vec4f m_lin_vel; // linear velocity (v=dx/dt)
			t_vec4f m_ang_vel; // angular velocity
		};
	};
};

#endif

