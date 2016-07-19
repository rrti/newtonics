#ifndef NEWTONICS_ARG_PARSER_HDR
#define NEWTONICS_ARG_PARSER_HDR

#include <cassert>
#include <cstdio>

#include <string>
#include <unordered_map>

typedef         int t_int_type;
typedef       float t_flt_type;
typedef std::string t_str_type;


struct t_arg_obj {
public:
	t_arg_obj() {}
	t_arg_obj(const t_str_type& name, const t_str_type& value) {
		assert(!name.empty());
		assert(!value.empty());

		set_name(name);
		set_value(value);
	}

	const t_str_type& get_name() const { return m_name; }
	const t_str_type& get_value() const { return m_value; }

	t_int_type to_int() const { return std::atoi(m_value.c_str()); }
	t_flt_type to_flt() const { return std::atof(m_value.c_str()); }
	t_str_type to_str() const { return          (m_value        ); }

	void set_name(const t_str_type& n) { m_name = n; }
	void set_value(const t_str_type& v) { m_value = v; }

private:
	t_str_type m_name;
	t_str_type m_value;
};


struct t_arg_parser {
public:
	// insert one dummy arg
	t_arg_parser() { add_arg_obj(t_arg_obj("$dummy$", "$value$")); }

	t_int_type get_int_arg(const t_str_type& arg_name) const { return ((get_arg_obj(arg_name)).to_int()); }
	t_flt_type get_flt_arg(const t_str_type& arg_name) const { return ((get_arg_obj(arg_name)).to_flt()); }
	t_str_type get_str_arg(const t_str_type& arg_name) const { return ((get_arg_obj(arg_name)).to_str()); }

	void add_arg_obj(const t_arg_obj& arg_obj) {
		m_args.insert(std::make_pair(arg_obj.get_name(), arg_obj));
	}

	// make an alias map to the name of an arg; multiple aliases
	// can map to the same arg (but every alias string should be
	// unique)
	void add_arg_alias(const t_str_type& arg_name, const t_str_type& arg_alias) {
		m_arg_aliases.insert(std::make_pair(arg_alias, arg_name));
	}

	void set_arg_val(const t_str_type& arg_name, const t_str_type& arg_val) {
		assert(m_args.find(arg_name) != m_args.end());
		m_args[arg_name].set_value(arg_val);
	}

	const t_arg_obj& get_arg_obj(const t_str_type& arg_name, size_t depth = 0) const {
		const auto iter = m_args.find(arg_name);

		if (iter != m_args.end())
			return (iter->second);

		// do not recurse forever
		assert(depth == 0);

		return (get_arg_obj("$dummy$", depth + 1));
	}


	bool match_arg_name(const t_str_type& arg_name, t_arg_obj* ret) const {
		const auto iter = m_args.find(arg_name);

		if (iter != m_args.end()) {
			*ret = iter->second;
			return true;
		}

		return false;
	}

	bool match_arg(const char* arg_str, const char* arg_val) {
		t_arg_obj arg_obj;

		if (match_arg_name(arg_str, &arg_obj)) {
			set_arg_val(arg_obj.get_name(), arg_val);
			return true;
		}

		return false;
	}

	void parse_args(int argc, const char** argv) {
		if (argc <= 1)
			return;

		// skip program-name arg
		argv += (argc > 0);
		argc -= (argc > 0);

		// scan the CL args
		for (int i = 0; i < (argc - 1); i += 2) {
			const char* arg_str = argv[i    ];
			const char* arg_val = argv[i + 1];

			// check if argv[i] is an alias for a named arg
			const auto alias_it = m_arg_aliases.find(t_str_type(arg_str));

			// if so, substitute the name of the aliased arg
			if (alias_it != m_arg_aliases.end())
				arg_str = (alias_it->second).c_str();

			match_arg(arg_str, arg_val);
		}
	}

	void print_args() const {
		printf("[%s]\n", __func__);

		for (auto it = m_args.cbegin(); it != m_args.cend(); ++it) {
			printf("\t\"%s\"=%s\n", (it->first).c_str(), ((it->second).get_value()).c_str());
		}
	}

private:
	std::unordered_map<t_str_type, t_arg_obj> m_args;

	// maps keyed shortcuts to arg-names
	std::unordered_map<t_str_type, t_str_type> m_arg_aliases;
};

#endif

