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
	int* amaya_mtbdd_get_transition_target(
			sylvan::MTBDD mtbdd, 
			uint8_t* cube, 
			uint32_t cube_size, 
			uint32_t* result_size);

	void amaya_print_dot(sylvan::MTBDD m, int32_t fd);
	sylvan::MTBDD amaya_mtbdd_build_single_terminal(
		uint8_t *transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t transition_symbols_count,
		uint32_t variable_count,
		uint32_t *destination_set,
		uint32_t destination_set_size);
	void amaya_mtbdd_rename_states(
			sylvan::MTBDD root, 
			int* names, // [(old, new), (old, new), (old, new)] 
			uint32_t name_count);

	int* amaya_mtbdd_get_leaves(
			sylvan::MTBDD root, 
			uint32_t** leaf_sizes,	// OUT, Array containing the sizes of leaves inside dest
			uint32_t* leaf_cnt);	// OUT, Number of leaves in the tree

	void amaya_do_free(void *ptr);

	void shutdown_machinery();
	void init_machinery();
}

std::set<int>* _get_transition_target(
		sylvan::MTBDD root, 
		uint32_t current_variable,
		uint8_t* variable_assigments, 
		uint32_t var_count);

void collect_mtbdd_leaves(sylvan::MTBDD root, std::set<sylvan::MTBDD>& dest);
#endif
