#ifndef NEWTONICS_BOUNDING_BOX_HDR
#define NEWTONICS_BOUNDING_BOX_HDR

#include "../math/point.hpp"
#include "../math/vector.hpp"

namespace newtonics {
	namespace lib_physics {
		using lib_math::t_pos4f;
		using lib_math::t_vec4f;

		struct t_bounding_box {
		public:
			t_bounding_box(): m_radius(0.0f) { clear(); }

			void clear() {
				m_mins = t_vec4f( 1e6f,  1e6f,  1e6f);
				m_maxs = t_vec4f(-1e6f, -1e6f, -1e6f);
			}

			void add_point(const t_pos4f& p) { add_point(p.to_vector()); }
			void add_point(const t_vec4f& p) {
				m_mins = p.min(m_mins);
				m_maxs = p.max(m_maxs);
			}

			float finalize() { return (m_radius = calc_radius()); }

			t_pos4f get_mpos() const { return ((m_mins + m_maxs) * 0.5f).to_point(); }
			t_vec4f get_size() const { return ((m_maxs - m_mins)                  ); }

			// size returns the full dimensions; we want its half-dims
			float calc_radius() const { return ((get_size()).len() * 0.5f); }
			float calc_diameter() const { return (calc_radius() * 2.0f); }

			// finalize before calling these
			float get_radius() const { return m_radius; }
			float get_diameter() const { return (get_radius() * 2.0f); }

		private:
			// bounding-box extrema
			t_vec4f m_mins;
			t_vec4f m_maxs;

			// radius of bounding-box
			float m_radius;
		};

		typedef t_bounding_box t_bbox;
	};
};

#endif

