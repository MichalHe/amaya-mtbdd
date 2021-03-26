#include <assert.h>
#include <cstring>
#include <iterator>
#include <sstream>
#include <sylvan.h>
#include <sylvan_bdd.h>
#include <sylvan_common.h>
#include <sylvan_gmp.h>
#include <sylvan_mt.h>
#include <sylvan_mtbdd.h>
#include <sylvan_int.h>
#include <sylvan_mtbdd_int.h>
#include <sylvan_stats.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <sstream>
#include <string>

#ifndef AMAYA_MTBDD_H
#define AMAYA_MTBDD_H
extern "C" {
	// Export constants (wrapped)
	extern const sylvan::MTBDD w_mtbdd_true = sylvan::mtbdd_true;
	extern const sylvan::MTBDD w_mtbdd_false = sylvan::mtbdd_false;


	// Functions
	sylvan::MTBDD amaya_unite_mtbdds(sylvan::MTBDD m1, sylvan::MTBDD m2);
	sylvan::MTBDD amaya_project_variables_away(
			sylvan::MTBDD m, uint32_t *variables, uint32_t var_count);

	void amaya_print_dot(sylvan::MTBDD m, int32_t fd);
	sylvan::MTBDD amaya_mtbdd_build_single_terminal(
		uint8_t *transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t transition_symbols_count,
		uint32_t variable_count,
		uint32_t *destination_set,
		uint32_t destination_set_size);

	void shutdown_machinery();
	void init_machinery();
}

#endif
