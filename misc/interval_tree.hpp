#ifndef NEWTONICS_INTERVAL_TREE_HDR
#define NEWTONICS_INTERVAL_TREE_HDR

#include <cassert>
#include <cstddef>

#include <set>

template<typename t_key_type, typename t_val_type> struct t_interval_node {
public:
	t_interval_node(
		t_key_type min_key = t_key_type(0),
		t_key_type max_key = t_key_type(0),
		t_val_type val_obj = t_val_type()
	) {
		assert(min_key <= max_key);

		m_min_key = min_key;
		m_max_key = max_key;
		m_val_obj = val_obj;
	}
	t_interval_node(const t_interval_node& i) {
		*this = i;
	}

	t_interval_node& operator = (const t_interval_node& i) {
		m_min_key = i.get_min_key();
		m_max_key = i.get_max_key();
		m_val_obj = i.get_val_obj();
		return *this;
	}

	#if 0
	bool operator < (const t_interval_node& i) const {
		// this treats the max-bound as exclusive, which
		// prevents us from inserting adjacent intervals
		//
		// e.g. when inserting [0, 16] and [16, 32], the
		// latter will be rejected, but it allows finding
		// a point-interval on boundaries
		return (m_max_key < i.get_min_key());
	}
	#else
	bool operator < (const t_interval_node& i) const {
		return (m_max_key <= i.get_min_key());
	}
	#endif
	bool operator == (const t_interval_node& i) const {
		return (m_min_key == i.get_min_key() && m_max_key == i.get_max_key());
	}

	bool overlaps(const t_interval_node& i) const {
		if (m_max_key <= i.get_min_key()) return false;
		if (m_min_key >= i.get_max_key()) return false;
		return true;
	}

	t_key_type get_min_key() const { return m_min_key; }
	t_key_type get_max_key() const { return m_max_key; }

	void set_val_obj(const t_val_type& m_val_obj) { m_val_obj = m_val_obj; }
	const t_val_type& get_val_obj() const { return m_val_obj; }

	static t_interval_node point_interval(t_key_type key) {
		// note: the "+1" allows points to lie on boundaries
		return (t_interval_node(key, key + t_key_type(1)));
	}

private:
	t_key_type m_min_key;
	t_key_type m_max_key;
	t_val_type m_val_obj;
};

template<typename t_interval_node_type> struct t_interval_tree {
public:
	typedef  typename std::set<t_interval_node_type>::iterator  t_interval_node_iter;
	typedef  typename std::set<t_interval_node_type>::const_iterator  t_interval_node_const_iter;
	typedef  typename std::set<t_interval_node_type>::size_type  t_interval_size_type;
	typedef  typename std::pair<t_interval_node_iter, bool>  t_interval_iter_bool;

	bool empty() const { return (m_tree.empty()); }
	size_t size() const { return (m_tree.size()); }
	void clear() { m_tree.clear(); }

	// t_interval_iter_bool insert(const t_interval_node_type& n) { return (m_tree.insert(n)); }
	t_interval_node_iter insert(const t_interval_node_type& n) { return (m_tree.insert(n)); }
	t_interval_size_type erase(const t_interval_node_type& n) { return (m_tree.erase(n)); }

	t_interval_node_iter find(const t_interval_node_type& n) const { return (m_tree.find(n)); }
	t_interval_node_iter begin() const { return (m_tree.begin()); }
	t_interval_node_iter end() const { return (m_tree.end()); }

	t_interval_node_const_iter cbegin() const { return (m_tree.cbegin()); }
	t_interval_node_const_iter cend() const { return (m_tree.cend()); }

private:
	// STL implements sets with ordered trees, all we need is a comparator
	std::multiset<t_interval_node_type> m_tree;
};

#endif

