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
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include <string>

#ifndef AMAYA_MTBDD_H
#define AMAYA_MTBDD_H

#define AMAYA_EXISTS_OPERATION_ID 		11
#define AMAYA_UNION_OP_ID 				12
#define AMAYA_INTERSECTION_OP_ID 		13

extern "C" {
	static uint32_t CUR_PADDING_CLOSURE_ID = 64;

	// Export constants (wrapped)
	extern const sylvan::MTBDD w_mtbdd_true = sylvan::mtbdd_true;
	extern const sylvan::MTBDD w_mtbdd_false = sylvan::mtbdd_false;


	// Functions
	sylvan::MTBDD amaya_unite_mtbdds(sylvan::MTBDD m1, sylvan::MTBDD m2, uint32_t automaton_id);
	sylvan::MTBDD amaya_project_variables_away(
			sylvan::MTBDD m, uint32_t *variables, uint32_t var_count);
	int* amaya_mtbdd_get_transition_target(
			sylvan::MTBDD mtbdd, 
			uint8_t* cube, 
			uint32_t cube_size, 
			uint32_t* result_size);

	void amaya_print_dot(sylvan::MTBDD m, int32_t fd);

	sylvan::MTBDD amaya_mtbdd_build_single_terminal(
		uint32_t  automaton_id,
		uint8_t*  transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t  transition_symbols_count,
		uint32_t  variable_count,
		int* 	  destination_set,
		uint32_t  destination_set_size);

	void amaya_mtbdd_rename_states(
			sylvan::MTBDD* mtbdd_roots, 
			uint32_t root_count,
			int* names, // [(old, new), (old, new), (old, new)] 
			uint32_t name_count);

	/**
	 * @param leaf_ptrs 	If not NULL will point to an array containing pointers 
	 * 						to transition destinations set of the leaves.
	 */
	int* amaya_mtbdd_get_leaves(
			sylvan::MTBDD root, 
			uint32_t** leaf_sizes,	// OUT, Array containing the sizes of leaves inside dest
			uint32_t*  leaf_cnt, 	// OUT, Number of leaves in the tree
			void***    leaf_ptrs);

	/**
	 * Given a valid pointer to a leaf transition destination set, replaces its contents.
	 * You can receive the pointers by calling amaya_mtbdd_get_leaves and passing leaves_ptr as
	 * non NULL.
	 * @param leaf_tds 		Transition destination set pointer.
	 * @param new_contents  Array containin the new contents.
	 * @param contents_size The length of the array.
	 */
	void amaya_replace_leaf_contents_with(
			void *leaf_tds, 
			int* new_contents, 
			uint32_t contents_size);

	/**
	 * Given an arary of MTBDD roots, collects their leaves and changes their automaton id
	 * to `new_id`. Should be called after a new automaton is a result of automaton union.
	 *
	 * @param roots 	The array of mtbdds roots.
	 * @param root_cnt 	The number of roots in the `roots` array.
	 * @param new_id 	New automaton identifier for the leaves.
	 */
	void amaya_mtbdd_change_automaton_id_for_leaves(
			sylvan::MTBDD* roots,
			uint32_t root_cnt,
			uint32_t new_id);

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
			int left_state,
			sylvan::MTBDD left, 
			int right_state,
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
	 * @param a 				Some transition MTBDD
	 * @param b 				Other transition MTBDD
	 * @param new_automaton_id 	Automaton ID for the created intersect leaves.
	 *
	 * @param discovered_states 		Flat array of discovered metastates and their generated int.
	 * 									[metastate_left, metastate_right, state_int ...]
	 * @param discovered_states_cnt 	Discovered_states count (size/3).
	 * @returns A MTBDD that contains only transitions that are present in both given MTBDDs.
	 */
	sylvan::MTBDD amaya_mtbdd_intersection(
		sylvan::MTBDD a, 
		sylvan::MTBDD b,
        uint32_t result_automaton_id, 
        int** discovered_states,         // OUT
        int*  discovered_states_cnt);    // OUT

	
	/**
	 * Marks the beginning of intersection. Setups global std::map that holds
	 * information about which intersection metastates (pairs) were mapped to which values.
	 */
	void amaya_begin_intersection();

	/**
	 * Updates the intersection state information by inserting provided metastates (flat 2d array)
	 * and their corresponding mappings into intersection state.
	 *
	 * @param metastates  			Flattened 2D array containing the intersection metastates (always a pair)
	 * @param renamed_metastates  	Renamed metastates.
	 * @param cnt 					Count of the metastates.
	 */
	void amaya_update_intersection_state(
			int* metastates, 
			int* renamed_metastates, 
			uint32_t cnt);
	/**
	 * Marks the end of intersection and frees up resources.
	 */
	void amaya_end_intersection();

	void amaya_set_debugging(bool debug);

	void amaya_do_free(void *ptr);

	/** 
	 * Renames the metastates that are contained withing the roots of mtbdds resulting 
	 * from the determinization procedure some automaton.
	 * @param roots 					The roots of the MTBDDs that were created during the determinization procedure.  
	 * @param root_cnt  				The number of given MTBDDs. 
	 * @param resulting_automaton_id  	The ID of the resulting automaton.
	 * @param out_metastates_sizes   	OUTPUT: The sized of the located metastates.
	 * @param out_metastates_cnt   		OUTPUT: The number of located metastates.
	 * @returns The array containing serialized metastates in the order they were located in the given MTBDDs.
	 */
	int* amaya_rename_metastates_to_int(
			sylvan::MTBDD* roots, 							// MTBDDs resulting from determinization
			uint32_t root_cnt,								// Root count
			int metastate_num_range_start,
			uint32_t resulting_automaton_id,
			uint32_t** out_metastates_sizes,
			uint32_t*  out_metastates_cnt
			);

	sylvan::MTBDD amaya_complete_mtbdd_with_trapstate(
			sylvan::MTBDD mtbdd,
			uint32_t automaton_id, 
			int trapstate,
			bool* had_effect
			);

	void shutdown_machinery();
	void init_machinery();
}

class Transition_Destination_Set {
public:
	std::set<int>* destination_set;
	uint32_t automaton_id;

	Transition_Destination_Set();
	Transition_Destination_Set(const Transition_Destination_Set &other);
	Transition_Destination_Set(uint32_t automaton_id, std::set<int>* destination_set);
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
	int left_state; 	// For debug purposes
	int right_state; 	// Actually used
	uint32_t final_states_cnt;
	int *final_states;
} pad_closure_info_t;

typedef struct {
	uint32_t automaton_id;
	int trapstate;
	bool had_effect;
} Complete_With_Trapstate_Op_Info;

typedef struct {
	uint32_t automaton_id;
	std::vector<int>* discoveries; // Flat array of [metastate_left, metastate_right, state, ...]
} intersect_info_t;

void collect_mtbdd_leaves(sylvan::MTBDD root, std::set<sylvan::MTBDD>& dest);
#endif
