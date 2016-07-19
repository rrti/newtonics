#include <cassert>
#include <cstring>

#include <algorithm>
#include <stack>
#include <unordered_map>

#include "../global_defs.hpp"
#include "scoped_timer.hpp"
#include "highres_clock.hpp"

#define STRING_TO_SIZE_T(s) reinterpret_cast<     size_t>(s)
#define SIZE_T_TO_STRING(s) reinterpret_cast<const char*>(s)

typedef std::stack<t_scoped_timer*> t_scoped_timer_stack;
typedef std::unordered_map<size_t, t_scoped_timer::t_timing_data> t_hashed_timing_table;
typedef std::pair<const char*, t_scoped_timer::t_timing_data> t_timing_entry;
typedef std::vector<t_timing_entry> t_sorted_timing_table;

static t_scoped_timer_stack* g_scoped_timer_stack_ref = nullptr;
static t_hashed_timing_table* g_hashed_timing_data_ref = nullptr;
static t_sorted_timing_table* g_sorted_timing_data_ref = nullptr;


static int extract_delim_index(const char* s, int i) {
	int k = i;

	while (k > 0) {
		if (s[k - 1] == ':' && s[k] == ':')
			break;
		if (s[k] == ' ')
			break;

		k--;
	}

	switch (s[k]) {
		case ':': { k--; } break;
		case ' ': { k++; } break;
	}

	return k;
}

static int compare_timing_data_pairs(const t_timing_entry& a, const t_timing_entry& b) {
	assert(a.first != nullptr);
	assert(b.first != nullptr);

	return ((a.second).stime > (b.second).stime);
}

static void print_timing_data_entry(const t_timing_entry& raw_entry) {
	const char* name = raw_entry.first;

	const t_scoped_timer::t_timing_data& raw_data = raw_entry.second;
	// const t_scoped_timer::t_timing_data& sum_data = sum_entry.second;

	const int64_t stime_ms = raw_data.stime / (1000 * 1000);
	const int64_t ctime_ms = raw_data.ctime / (1000 * 1000);
	const int64_t calls    = raw_data.calls;

	// name = "[[... ::]class]::func(...)"
	// look up the first '(', then scan backwards twice for
	// "::" or " " to get the function and class substrings
	const int i = std::max(0l, strstr(name, "(") - name);
	const int n = extract_delim_index(name, i);
	const int k = extract_delim_index(name, n);
	const int j = (name[k] == ':') * 2;

	char buf[512] = {0};

	if (i > 0) {
		memcpy(buf, name + (k + j), i - k - j);
	} else {
		memcpy(buf, name, strlen(name));
	}

	printf("[%s] :: time=%ldms (children=%ldms) calls=%ld depth=%ld\n", buf, stime_ms, ctime_ms, calls, raw_data.depth);
}



#if (ENABLE_SCOPED_TIMERS == 1)
t_scoped_timer::t_scoped_timer(const char* scope_name) {
	m_parent_timer = (!g_scoped_timer_stack_ref->empty())? g_scoped_timer_stack_ref->top(): nullptr;
	m_scope_name = scope_name;

	m_ctor_tick = t_hres_clock::get_ticks();
	m_dtor_tick = m_ctor_tick;

	g_scoped_timer_stack_ref->push(this);
}

t_scoped_timer::~t_scoped_timer() {
	m_dtor_tick = t_hres_clock::get_ticks();

	// create or update the entry for this scope
	// (reinterprets const char* as an integer so
	// we do not incur the cost of string-hashes)
	t_hashed_timing_table& table = *g_hashed_timing_data_ref;
	t_timing_data& data = table[ STRING_TO_SIZE_T(m_scope_name) ];

	// track cumulative time (in ns) and number of calls
	data.stime += (m_dtor_tick - m_ctor_tick);
	data.calls += 1;
	data.depth  = std::max(data.depth, int64_t(g_scoped_timer_stack_ref->size()));

	// add elapsed time to parent so it can calculate the
	// difference between its own scope and that spent in
	// child-scopes
	if (m_parent_timer != nullptr)
		table[ STRING_TO_SIZE_T(m_parent_timer->get_scope_name()) ].ctime += (m_dtor_tick - m_ctor_tick);

	assert(!g_scoped_timer_stack_ref->empty());
	g_scoped_timer_stack_ref->pop();
}

#else
t_scoped_timer::t_scoped_timer(const char*) {}
t_scoped_timer::~t_scoped_timer() {}
#endif



t_scoped_timer::t_timing_data t_scoped_timer::get_timing_data(const char* scope_name) {
	const auto it = g_hashed_timing_data_ref->find(STRING_TO_SIZE_T(scope_name));

	if (it != g_hashed_timing_data_ref->end())
		return it->second;

	return (t_timing_data());
}

void t_scoped_timer::init_timing_data() {
	static t_scoped_timer_stack g_scoped_timer_stack;
	static t_hashed_timing_table g_hashed_timing_data;
	static t_sorted_timing_table g_sorted_timing_data;

	g_scoped_timer_stack_ref = &g_scoped_timer_stack;
	g_hashed_timing_data_ref = &g_hashed_timing_data;
	g_sorted_timing_data_ref = &g_sorted_timing_data;
}

void t_scoped_timer::print_timing_data() {
	printf("[%s]\n", __PRETTY_FUNCTION__);

	if (g_hashed_timing_data_ref->empty())
		return;

	g_sorted_timing_data_ref->clear();

	// copy; turn timing-data keys back into const char*'s
	for (auto it = g_hashed_timing_data_ref->cbegin(); it != g_hashed_timing_data_ref->cend(); ++it) {
		g_sorted_timing_data_ref->push_back(std::pair<const char*, t_timing_data>(SIZE_T_TO_STRING(it->first), it->second));
	}

	if (!g_sorted_timing_data_ref->empty())
		std::sort(g_sorted_timing_data_ref->begin(), g_sorted_timing_data_ref->end(), compare_timing_data_pairs);

	for (const t_timing_entry& entry: *g_sorted_timing_data_ref) {
		print_timing_data_entry(entry);
	}
}

