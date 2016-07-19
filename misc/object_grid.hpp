#ifndef NEWTONICS_OBJECT_GRID_HDR
#define NEWTONICS_OBJECT_GRID_HDR

#include <cassert>
#include <vector>

#include "../math/ray.hpp"
#include "../math/tuple.hpp"
#include "../math/point.hpp"
#include "../math/vector.hpp"

using namespace newtonics;
using namespace newtonics::lib_math;



// light-weight proxy for objects inserted into grid
struct t_grid_object {
public:
	typedef typename std::pair<size_t, size_t> t_cell_data_pair;
	typedef typename std::vector<t_cell_data_pair> t_cell_data_vector;

	t_grid_object(size_t id = 0, size_t mn = 0) {
		set_unique_id(id);
		set_magic_num(mn);
	}

	bool operator == (const t_grid_object& o) const { return (m_unique_id == o.get_unique_id()); }
	bool operator != (const t_grid_object& o) const { return (m_unique_id != o.get_unique_id()); }
	bool operator  < (const t_grid_object& o) const { return (m_unique_id  < o.get_unique_id()); }

	void set_unique_id(size_t id) { m_unique_id = id; }
	void set_magic_num(size_t mn) { m_magic_num = mn; }

	void update_cell_data(const t_cell_data_pair& p) {
		for (auto it = m_cells.begin(); it != m_cells.end(); ++it) {
			if (it->first == p.first) {
				it->second = p.second;
				break;
			}
		}
	}

	size_t get_unique_id() const { return m_unique_id; }
	size_t get_magic_num() const { return m_magic_num; }

	const t_cell_data_vector& get_cell_data() const { return m_cells; }
	      t_cell_data_vector& get_cell_data()       { return m_cells; }

private:
	size_t m_unique_id;
	size_t m_magic_num;

	// holds (cell-idx, object-idx) pairs for all cells we are in
	// (reference the std::make_pair<> call in grid::add_object)
	t_cell_data_vector m_cells;
};



class t_object_grid {
public:
	typedef std::vector<t_grid_object*> t_object_vector;

	struct t_grid_cell {
	public:
		t_grid_cell(): m_used_idx(-1lu) {}

		size_t add_object(t_grid_object* object) {
			m_grid_objects.push_back(object);
			return (m_grid_objects.size() - 1);
		}

		void del_object(const t_grid_object::t_cell_data_pair& p) {
			assert(p.second < m_grid_objects.size());

			m_grid_objects[p.second] = m_grid_objects.back();
			m_grid_objects.pop_back();

			// repair the replaced object's (object-)index
			//
			// this is faster than deletion by pointer if an
			// object covers fewer cells on average than the
			// number of objects held per cell
			m_grid_objects[p.second]->update_cell_data(p);
		}

		#if 0
		// delete by pointer; does not repair
		void del_object(const t_grid_object* obj) {
			t_object_vector::iterator it = std::find(m_grid_objects.begin(), m_grid_objects.end(), obj);

			assert(it != m_grid_objects.end());

			*it = m_grid_objects.back();
			m_grid_objects.pop_back();
		}
		#endif

		void set_used_idx(size_t idx) { m_used_idx = idx; }
		size_t get_used_idx() const { return m_used_idx; }

		const t_object_vector& get_objects() const { return m_grid_objects; }

		const t_grid_object* get_object(size_t idx) const { return m_grid_objects[idx]; }
		      t_grid_object* get_object(size_t idx)       { return m_grid_objects[idx]; }

		size_t get_size() const { return (m_grid_objects.size()); }

		bool is_empty() const { return (m_grid_objects.empty()); }
		bool is_used() const { return (m_used_idx != -1lu); }

	private:
		// vector of proxy-objects present in this cell
		t_object_vector m_grid_objects;

		// index into grid::m_used_cells (which non-empty cell we are)
		size_t m_used_idx;
	};


	t_object_grid() { m_magic_num = 0;}
	t_object_grid(const t_tup4i& gdim, const t_vec4f& gmins, const t_vec4f& gmaxs): m_gdim(gdim), m_mins(gmins), m_maxs(gmaxs) {
		m_grid_cells.resize(m_gdim.x() * m_gdim.y() * m_gdim.z(), t_grid_cell());

		assert(m_gdim.x() > 0);
		assert(m_gdim.y() > 0);
		assert(m_gdim.z() > 0);

		m_cdim.x() = (m_maxs.x() - m_mins.x()) / m_gdim.x();
		m_cdim.y() = (m_maxs.y() - m_mins.y()) / m_gdim.y();
		m_cdim.z() = (m_maxs.z() - m_mins.z()) / m_gdim.z();

		m_magic_num = 0;
	}
	~t_object_grid() {
		m_grid_cells.clear();
		m_used_cells.clear();
	}



	size_t get_cells_on_ray(const t_ray& ray, std::vector<t_grid_cell>& cells) const {
		cells.reserve(16);

		const t_pos4f f_rpos = ray.pos();
		const t_vec4f f_rdir = ray.dir();

		// assert(f_rpos != t_pos4f::error_point());
		assert(f_rdir != t_vec4f::zero_vector());

		constexpr int64_t epsil    = 1000l;
		constexpr int64_t scale    = 1000l * 1000l;
		constexpr int64_t scale_sq = scale * scale;

		// handle diagonal corner-to-corner traces
		// (note: we need *2* iterations per cell!)
		const unsigned int max_iters = std::sqrt(square(m_gdim.x()) + square(m_gdim.y()) + square(m_gdim.z())) + 1;
		      unsigned int num_iters = 0;

		const t_vec4li i_cdim(m_cdim.x() * scale, m_cdim.y() * scale, m_cdim.z() * scale);
		const t_vec4li i_rdir(f_rdir.x() * scale, f_rdir.y() * scale, f_rdir.z() * scale);

		const t_vec4li i_mins(m_mins.x() * scale, m_mins.y() * scale, m_mins.z() * scale);
		const t_vec4li i_maxs(m_maxs.x() * scale, m_maxs.y() * scale, m_maxs.z() * scale);
		const t_vec4li i_size(i_maxs - i_mins);

		t_vec4li i_rvec;
		t_pos4li i_rpos(f_rpos.x() * scale, f_rpos.y() * scale, f_rpos.z() * scale);

		t_tup4i i_cidx(((i_rpos.x() - i_mins.x()) * m_gdim.x()) / i_size.x(), ((i_rpos.y() - i_mins.y()) * m_gdim.y()) / i_size.y(), ((i_rpos.z() - i_mins.z()) * m_gdim.z()) / i_size.z());
		t_tup4i j_cidx;

		// floating-point trace:
		//   dim    = (x = 1.0, y = 1.0)
		//   pos0   = (x = 3.1, y = 2.1)
		//   dir    = (x = 0.8, y = 0.2)
		//   dist0  = {x = 0.9, y = 0.9)
		//   step0  = (x = dist.x / dir.x, y = dist.y / dir.y) = (1.125, 4.5) (min=1.125)
		//   vec0   = dir * step0 = (0.8, 0.2) * 1.125 = (0.9, 0.225)
		//   pos1   = pos0 + vec0 = (3.1, 2.1) + (0.9, 0.225) = (4.0, 2.325)
		//
		// fixed-point trace (note the extra scaling of dist in step):
		//   scale  = 1000000
		//   dim    = (x = 1.0 * scale, y = 1.0 * scale) = (1000000, 1000000)
		//   pos0   = (x = 3.1 * scale, y = 2.1 * scale) = (3100000, 2100000)
		//   dir    = (x = 0.8 * scale, y = 0.2 * scale) = ( 800000,  200000)
		//   dist0  = (x = 0.9 * scale, y = 0.9 * scale) = ( 900000,  900000)
		//   step0  = (x = (dist.x * scale) / dir.x, y = (dist.y * scale) / dir.y) = (1125000, 4500000) (min=1125000)
		//   vec0   = (dir * step0) / scale = ((800000, 200000) * 1125000) / 1000000 = (900000, 225000)
		//   pos1   = pos0 + vec0 = (3100000, 2100000) + (900000, 225000) = (4000000, 2325000)
		//
		while ((num_iters++ < max_iters) && idx_in_bounds(i_cidx)) {
			if (i_cidx != j_cidx)
				cells.push_back(get_cell(i_cidx));

			const t_tup4li i_cmin = {i_mins.x() + (i_cidx.x()    ) * i_cdim.x(), i_mins.y() + (i_cidx.y()    ) * i_cdim.y(), i_mins.z() + (i_cidx.z()    ) * i_cdim.z()};
			const t_tup4li i_cmax = {i_mins.x() + (i_cidx.x() + 1) * i_cdim.x(), i_mins.y() + (i_cidx.y() + 1) * i_cdim.y(), i_mins.z() + (i_cidx.z() + 1) * i_cdim.z()};

			// note: min_dst values are negative, will flip sign after division by dir
			const t_vec4li i_dmin = {i_cmin.x() - i_rpos.x(), i_cmin.y() - i_rpos.y(), i_cmin.z() - i_rpos.z()};
			const t_vec4li i_dmax = {i_cmax.x() - i_rpos.x(), i_cmax.y() - i_rpos.y(), i_cmax.z() - i_rpos.z()};

			t_vec4li i_dist;
			t_vec4li i_step;

			for (unsigned int n = 0; n < 3; n++) {
				switch (signum(i_rdir[n])) {
					case -1: { i_dist[n] = i_dmin[n]; i_step[n] = (i_dist[n] * scale) / i_rdir[n]; } break;
					case +1: { i_dist[n] = i_dmax[n]; i_step[n] = (i_dist[n] * scale) / i_rdir[n]; } break;
					case  0: { i_dist[n] =  scale_sq; i_step[n] =                        scale_sq; } break;
					default: {} break;
				}
			}

			i_step.w()  = min3(i_step.x(), i_step.y(), i_step.z());
			i_dist.w()  = 1;
			i_dist.w() *= (((i_rdir.x() * i_step.w()) / scale) == 0);
			i_dist.w() *= (((i_rdir.y() * i_step.w()) / scale) == 0);
			i_dist.w() *= (((i_rdir.z() * i_step.w()) / scale) == 0);

			i_rvec  = i_rdir * i_step.w();
			i_rpos += (i_rvec / scale);

			// assert(i_rvec.x() != 0 || i_rvec.y() != 0 || i_rvec.z() != 0);

			// this (very slightly) translates the ray each cell where .w=1
			// and still allows getting stuck *on* a cell boundary (dists 0)
			// i_rpos += (i_dist * i_dist.w());

			// nudge the ray across instead to make sure it always advances
			// (epsilon should be made smaller if the ratio of dir to csize
			// is larger, otherwise we risk skipping cells)
			i_rpos += (((i_rdir * epsil) / scale) * i_dist.w());

			j_cidx = i_cidx;
			i_cidx = t_tup4i(((i_rpos.x() - i_mins.x()) * m_gdim.x()) / i_size.x(), ((i_rpos.y() - i_mins.y()) * m_gdim.y()) / i_size.y(), ((i_rpos.z() - i_mins.z()) * m_gdim.z()) / i_size.z());
		}

		return (cells.size());
	}

	// get all cells in the sphere defined by <position = wpos.xyz, radius = wpos.w>
	size_t get_cells_in_sphere(const t_tup4f& wpos, std::vector<t_grid_cell>& cells) const {
		const t_pos4f raw_pos = t_pos4f(wpos.x(), wpos.y(), wpos.z());

		const t_tup4i&  cell_idx = get_cell_idx(raw_pos);
		const t_vec4f&  cell_dim = get_cell_size();
		const t_tup4i  num_cells = t_tup4i((wpos.w() / cell_dim.x()) + 1, (wpos.w() / cell_dim.y()) + 1, (wpos.w() / cell_dim.z()) + 1);

		cells.reserve(16);

		for (int x = cell_idx.x() - num_cells.x(); x <= cell_idx.x() + num_cells.x(); x++) {
			for (int y = cell_idx.y() - num_cells.y(); y <= cell_idx.y() + num_cells.y(); y++) {
				for (int z = cell_idx.z() - num_cells.z(); z <= cell_idx.z() + num_cells.z(); z++) {
					const t_tup4i cell_idx = t_tup4i(x, y, z);

					if (!idx_in_bounds(cell_idx))
						continue;

					// get the point on cell that is closest to <position>; the
					// distance to this closest point must be less than <radius>
					const t_pos4f& cell_pos = clamp_pos_in_cell(raw_pos, cell_idx);

					// radius-check
					if ((cell_pos - raw_pos).sq_len() > (wpos.w() * wpos.w()))
						continue;

					cells.push_back(get_cell(cell_idx));
				}
			}
		}

		return (cells.size());
	}

	// get all *cells* in the cube defined by <wmins> and <wmaxs>
	size_t get_cells_in_cube(const t_pos4f& wmins, const t_pos4f& wmaxs, std::vector<t_grid_cell>& cells) const {
		t_pos4f mins;
		t_pos4f maxs;

		// ensure min-bounds are always smaller than max-bounds
		mins.x() = std::min(wmins.x(), wmaxs.x());
		mins.y() = std::min(wmins.y(), wmaxs.y());
		mins.z() = std::min(wmins.z(), wmaxs.z());
		maxs.x() = std::max(wmins.x(), wmaxs.x());
		maxs.y() = std::max(wmins.y(), wmaxs.y());
		maxs.z() = std::max(wmins.z(), wmaxs.z());

		// calculate indices *without* clamping; only that part
		// of the cube intersecting the grid will be considered
		const t_tup4i& min_idx = get_cell_idx(mins);
		const t_tup4i& max_idx = get_cell_idx(maxs);

		cells.clear();
		cells.reserve((max_idx.x() - min_idx.x()) * (max_idx.y() - min_idx.y()) * (max_idx.z() - min_idx.z()));

		for (int x = min_idx.x(); x <= max_idx.x(); x++) {
			for (int y = min_idx.y(); y <= max_idx.y(); y++) {
				for (int z = min_idx.z(); z <= max_idx.z(); z++) {
					const t_tup4i cell_idx = t_tup4i(x, y, z);

					if (!idx_in_bounds(cell_idx))
						continue;

					cells.push_back(get_cell(cell_idx));
				}
			}
		}

		return (cells.size());
	}

	// get all *objects* in the cube defined by <wmins> and <wmaxs>
	size_t get_objects_in_cube(const t_pos4f& wmins, const t_pos4f& wmaxs, t_object_vector& objs) {
		std::vector<t_grid_cell> cells;

		if (get_cells_in_cube(wmins, wmaxs, cells) == 0)
			return 0;

		objs.reserve(cells.size());

		// must always be increased before comparing objects against it
		m_magic_num += 1;

		for (size_t n = 0; n < cells.size(); n++) {
			const t_grid_cell& cell = cells[n];
			const t_object_vector& cell_objs = cell.get_objects();

			for (auto it = cell_objs.cbegin(); it != cell_objs.cend(); ++it) {
				t_grid_object* obj = *it;

				if (obj->get_magic_num() == m_magic_num)
					continue;

				// mark object as processed (filters duplicates)
				obj->set_magic_num(m_magic_num);
				objs.push_back(obj);
			}
		}

		return (objs.size());
	}

	// get all *objects* in the cube covering <wpos> to <wpos + wvel>
	// like get_cells_in_sphere this also includes the radius wpos.w
	size_t get_objects_in_cube_vel(const t_pos4f& wpos, const t_vec4f& wvel, t_object_vector& objs) {
		const t_vec4f cdim = get_cell_size();
		const t_vec4i ccnt = t_vec4i((wvel.x()    ) / cdim.x(), (wvel.y()    ) / cdim.y(), (wvel.z()    ) / cdim.z());
		const t_vec4f crng = t_vec4f((ccnt.x() + 1) * cdim.x(), (ccnt.y() + 1) * cdim.y(), (ccnt.z() + 1) * cdim.z());

		// grid ensure mins < maxs even when vel < 0
		const t_pos4f wmins = {
			lerp(wpos.x() - wpos.w(), wpos.x() - std::max(wpos.w(), crng.x()), float(wvel.x() <= 0.0f)),
			lerp(wpos.y() - wpos.w(), wpos.y() - std::max(wpos.w(), crng.y()), float(wvel.y() <= 0.0f)),
			lerp(wpos.z() - wpos.w(), wpos.z() - std::max(wpos.w(), crng.z()), float(wvel.z() <= 0.0f))
		};
		const t_pos4f wmaxs = {
			lerp(wpos.x() - wpos.w(), wpos.x() + std::max(wpos.w(), crng.x()), float(wvel.x() > 0.0f)),
			lerp(wpos.y() - wpos.w(), wpos.y() + std::max(wpos.w(), crng.y()), float(wvel.y() > 0.0f)),
			lerp(wpos.z() - wpos.w(), wpos.z() + std::max(wpos.w(), crng.z()), float(wvel.z() > 0.0f))
		};

		return (get_objects_in_cube(wmins, wmaxs, objs));
	}

	// get all *objects* in the cube of cells within <wdims> of <wpos>
	// (<wdims> should contain the *half*-sizes of the cube to search)
	size_t get_objects_in_cube(const t_pos4f& wpos, const t_vec4f& wdims, t_object_vector& objs) {
		const t_tup4i&  cell_idx = get_cell_idx(wpos);
		const t_vec4f&  cell_dim = get_cell_size();
		const t_tup4i  num_cells = t_tup4i((wdims.x() / cell_dim.x()) + 1, (wdims.y() / cell_dim.y()) + 1, (wdims.z() / cell_dim.z()) + 1);

		objs.reserve(num_cells.x() * num_cells.y() * num_cells.z());

		// must always be increased before comparing objects against it
		m_magic_num += 1;

		for (int x = cell_idx.x() - num_cells.x(); x <= cell_idx.x() + num_cells.x(); x++) {
			for (int y = cell_idx.y() - num_cells.y(); y <= cell_idx.y() + num_cells.y(); y++) {
				for (int z = cell_idx.z() - num_cells.z(); z <= cell_idx.z() + num_cells.z(); z++) {
					const t_tup4i cell_idx = t_tup4i(x, y, z);

					if (!idx_in_bounds(cell_idx))
						continue;

					const t_grid_cell& cell = get_cell(cell_idx);
					const t_object_vector& cell_objs = cell.get_objects();

					for (auto it = cell_objs.cbegin(); it != cell_objs.cend(); ++it) {
						t_grid_object* obj = *it;

						if (obj->get_magic_num() == m_magic_num)
							continue;

						// mark object as processed (filters duplicates)
						obj->set_magic_num(m_magic_num);
						objs.push_back(obj);
					}
				}
			}
		}

		return (objs.size());
	}



	size_t add_object(t_grid_object* object, const t_tup4f& wpos) {
		// note: wpos.xyz = position, wpos.w = (bounding) radius
		const t_pos4f obj_pos = t_pos4f(              wpos.x(),               wpos.y(),               wpos.z());
		const t_pos4f min_pos = t_pos4f(obj_pos.x() - wpos.w(), obj_pos.y() - wpos.w(), obj_pos.z() - wpos.w());
		const t_pos4f max_pos = t_pos4f(obj_pos.x() + wpos.w(), obj_pos.y() + wpos.w(), obj_pos.z() + wpos.w());

		// calculate indices *with* clamping; each object
		// must always be registered somewhere within the
		// grid
		const t_tup4i min_idx = get_cell_idx_clamped(min_pos);
		const t_tup4i max_idx = get_cell_idx_clamped(max_pos);

		std::vector<t_grid_object::t_cell_data_pair>& obj_cells = object->get_cell_data();

		obj_cells.clear();
		obj_cells.reserve(16);

		// add object to all overlapped cells
		for (int x = min_idx.x(); x <= max_idx.x(); x++) {
			for (int y = min_idx.y(); y <= max_idx.y(); y++) {
				for (int z = min_idx.z(); z <= max_idx.z(); z++) {
					const t_tup4i cell_idx = t_tup4i(x, y, z);
					const t_pos4f& cell_pos = clamp_pos_in_cell(obj_pos, cell_idx);

					// radius-check
					if ((cell_pos - obj_pos).sq_len() > (wpos.w() * wpos.w()))
						continue;

					assert(idx_in_bounds(cell_idx));
					t_grid_cell& grid_cell = get_cell(cell_idx);

					if (grid_cell.is_empty()) {
						// empty cell, mark as occupied (since object will be added)
						assert(!grid_cell.is_used());

						m_used_cells.push_back(&grid_cell);
						grid_cell.set_used_idx(m_used_cells.size() - 1);

						assert(grid_cell.is_used());
					}

					// keep a record of the cells into which this object was inserted
					// this speeds up deletion if we want to relocate the object later
					obj_cells.push_back(std::make_pair(z * (m_gdim.y() * m_gdim.z()) + y * (m_gdim.y()) + x, grid_cell.add_object(object)));
				}
			}
		}

		return (obj_cells.size());
	}

	void del_object(t_grid_object* object) {
		// delete object from every cell it overlaps
		std::vector<t_grid_object::t_cell_data_pair>& obj_cells = object->get_cell_data();

		for (auto it = obj_cells.cbegin(); it != obj_cells.cend(); ++it) {
			t_grid_cell&  grid_cell =  m_grid_cells[it->first];
			t_grid_cell** used_cell = &m_used_cells[grid_cell.get_used_idx()];

			assert(object == grid_cell.get_object(it->second));

			// NOTE:
			//   this will change ordering of objects within cell's object-vector
			//   therefore indices returned by cell::add_object for other objects
			//   in <grid_cell> will become invalid (if stored in a vector)
			grid_cell.del_object(*it);

			if (!grid_cell.is_empty())
				continue;

			// cell is no longer used, return it to the free pool
			assert(grid_cell.is_used());
			assert(&grid_cell == *used_cell);

			(*used_cell) = m_used_cells.back();
			(*used_cell)->set_used_idx(grid_cell.get_used_idx());

			m_used_cells.pop_back();
			grid_cell.set_used_idx(-1lu);

			assert(!grid_cell.is_used());
		}

		obj_cells.clear();
	}



	bool idx_in_bounds(const t_tup4i& cidx) const {
		unsigned int mask = 0;
		mask += (cidx.x() < 0 || cidx.x() >= m_gdim.x());
		mask += (cidx.y() < 0 || cidx.y() >= m_gdim.y());
		mask += (cidx.z() < 0 || cidx.z() >= m_gdim.z());
		return (mask == 0);
	}
	// check if a position is within the grid-bounds
	// (by converting it to cell-indices and testing
	// those)
	bool pos_in_bounds(const t_pos4f& wpos) const {
		return (idx_in_bounds(get_cell_idx(wpos)));
	}
	#if 0
	bool pos_in_bounds(const t_pos4f& wpos) const {
		unsigned int mask = 0;
		mask += (wpos.x() < m_mins.x() || wpos.x() >= m_maxs.x());
		mask += (wpos.y() < m_mins.y() || wpos.y() >= m_maxs.y());
		mask += (wpos.z() < m_mins.z() || wpos.z() >= m_maxs.z());
		return (mask == 0);
	}
	#endif


	// <cpos> is the world-position of the cell
	t_pos4f clamp_pos_in_cell(const t_pos4f& wpos, const t_pos4f& cpos) const {
		return (clamp_pos_in_cell(wpos, get_cell_idx(cpos)));
	}

	t_pos4f clamp_pos_in_cell(const t_pos4f& wpos, const t_tup4i& cidx) const {
		const t_pos4f cell_offs = t_pos4f(cidx.x() * m_cdim.x(), cidx.y() * m_cdim.y(), cidx.z() * m_cdim.z());

		const t_pos4f cell_mins = t_pos4f(m_mins.x() + cell_offs.x(), m_mins.y() + cell_offs.y(), m_mins.z() + cell_offs.z());
		const t_pos4f cell_maxs = t_pos4f(cell_mins.x() + m_cdim.x(), cell_mins.y() + m_cdim.y(), cell_mins.z() + m_cdim.z());

		t_pos4f clamp_pos;

		clamp_pos.x() = lib_math::clamp(wpos.x(), cell_mins.x(), cell_maxs.x());
		clamp_pos.y() = lib_math::clamp(wpos.y(), cell_mins.y(), cell_maxs.y());
		clamp_pos.z() = lib_math::clamp(wpos.z(), cell_mins.z(), cell_maxs.z());

		return clamp_pos;
	}


	// return the index of the cell that contains spatial position <pos>
	// note: lack of numerical accuracy can cause the indices to be OOB!
	t_tup4i get_cell_idx(const t_pos4f& wpos) const {
		t_tup4i idx;
		idx.x() = int(((wpos.x() - m_mins.x()) / (m_maxs.x() - m_mins.x())) * m_gdim.x());
		idx.y() = int(((wpos.y() - m_mins.y()) / (m_maxs.y() - m_mins.y())) * m_gdim.y());
		idx.z() = int(((wpos.z() - m_mins.z()) / (m_maxs.z() - m_mins.z())) * m_gdim.z());
		return idx;
	}

	t_tup4i get_cell_idx_clamped(const t_pos4f& wpos) const {
		t_tup4i idx = get_cell_idx(wpos);
		idx.x() = lib_math::clamp(idx.x(), 0, m_gdim.x() - 1);
		idx.y() = lib_math::clamp(idx.y(), 0, m_gdim.y() - 1);
		idx.z() = lib_math::clamp(idx.z(), 0, m_gdim.z() - 1);
		return idx;
	}


	// return the cell containing spatial position <pos>
	// (auto-clamps the indices corresponding to <pos>)
	const t_grid_cell& get_cell(const t_pos4f& wpos) const { return (get_cell(get_cell_idx_clamped(wpos))); }
	      t_grid_cell& get_cell(const t_pos4f& wpos)       { return (get_cell(get_cell_idx_clamped(wpos))); }

	// return the cell at grid index <idx>; uses
	// indexing scheme z * (W * H) + y * (W) + x
	// note: assumes index-vector is in-bounds
	const t_grid_cell& get_cell(const t_tup4i& cidx) const { return (m_grid_cells[get_array_idx(cidx)]); }
	      t_grid_cell& get_cell(const t_tup4i& cidx)       { return (m_grid_cells[get_array_idx(cidx)]); }


	const std::vector<t_grid_cell*>& get_used_cells() const { return m_used_cells; }

	const t_tup4i& get_grid_size() const { return m_gdim; }
	const t_vec4f& get_cell_size() const { return m_cdim; }

	size_t get_array_idx(const t_tup4i& cidx) const {
		return (cidx.z() * (m_gdim.y() * m_gdim.z()) + cidx.y() * m_gdim.y() + cidx.x());
	}

private:
	std::vector<t_grid_cell> m_grid_cells;
	// pointers to all currently non-empty cells
	std::vector<t_grid_cell*> m_used_cells;

	// number of cells (divisions) along each dimension
	t_tup4i m_gdim;
	// spatial size per dimension of each cell
	t_vec4f m_cdim;

	// spatial bounds of the grid
	t_vec4f m_mins;
	t_vec4f m_maxs;

	// number of calls made to get_objects_in_*()
	size_t m_magic_num;
};

#endif

