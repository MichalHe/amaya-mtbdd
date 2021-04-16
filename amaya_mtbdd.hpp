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

#define AMAYA_PAD_CLOSURE_OPERATION_ID 10

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

	/**
	 * Calculates the post set from given MTBDD.
	 * @param m The transitions MTBDD for the state for which we calculate the post set.
	 * @param post_size (out) The size of the returned array.
	 * @returns array containing the post set.
	 */
	int* amaya_mtbdd_get_state_post(
			sylvan::MTBDD m, 
			uint32_t *post_size);

	/**
	 * Apply the padding closure to a state with transition function `left`, with successor's function `right`.
	 * @param left 			The transitions of state into which we want to propagate finishing symbols.
	 * @param right 		The transitions of a state that has some transitions leading to the final state.
	 * @param final_states 		Pointer to an array containing the final states of the automaton.
	 * @param final_states_cnt 	Number of states in the final_state array.
	 * @returns Boolean indicating whether the left mtbdd was modified (new transitions to final state were added).
	 */
	bool amaya_mtbdd_do_pad_closure(
			sylvan::MTBDD left, 
			sylvan::MTBDD right, 
			int* final_states, 
			uint32_t final_states_cnt);

	/**
	 * Debug function that allows to retrieve all paths inside the MTBDD and corresponding destinations.
	 * @param root The MTBDD itself.
	 * @param vars Array containing the variables of interest
	 * @param varcount Size of vars array
	 * @param symbols_cnt (Out) Will contain the number of symbols returned.
	 * @param dest_states (Out) Will point to an array where destination states are stored.
	 * @param dest_states_cnt (Out) Array of sizes for each destination set stored in dest_states.
	 * @returns Array of symbols serialized. Each symbol has the length of var_count, and the array contains symbols_cnt of such symbols.
	 */
	uint8_t* amaya_mtbdd_get_transitions(
			sylvan::MTBDD root,
			uint32_t* vars,	
			uint32_t var_count,
			uint32_t* symbols_cnt,
			int** dest_states,
			uint32_t** dest_states_cnt
			);

	/**
	 * Calculate the intersection of two mtbdds (Operation on leaves is set intersection).
	 * @param a Some transition MTBDD
	 * @param b Other transition MTBDD
	 * @returns A MTBDD that contains only transitions that are present in both given MTBDDs.
	 */
	sylvan::MTBDD amaya_mtbdd_intersection(
			sylvan::MTBDD a, 
			sylvan::MTBDD b);
	void amaya_do_free(void *ptr);

	void shutdown_machinery();
	void init_machinery();
}

class Transition_Destination_Set {
public:
	uint32_t last_state_rename_id;
	std::set<int>* destination_set;
	Transition_Destination_Set();
	Transition_Destination_Set(const Transition_Destination_Set &other);
	Transition_Destination_Set(std::set<int>* destination_set);
	~Transition_Destination_Set();
	void print_dest_states();
};

Transition_Destination_Set* _get_transition_target(
		sylvan::MTBDD root, 
		uint32_t current_variable,
		uint8_t* variable_assigments, 
		uint32_t var_count);

typedef struct {
	bool had_effect;
	uint32_t final_states_cnt;
	int *final_states;
} pad_closure_info_t;

void collect_mtbdd_leaves(sylvan::MTBDD root, std::set<sylvan::MTBDD>& dest);
#endif
