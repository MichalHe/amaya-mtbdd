#ifndef AMAYA_WRAPPER_H
#define AMAYA_WRAPPER_H

#include "base.hpp"
#include <sylvan.h>
#include <inttypes.h>

struct Serialized_DFA {
    State*          states;
    uint64_t        state_count;
    State           initial_state;
    State*          final_states;
    uint64_t        final_state_count;
    sylvan::MTBDD*  mtbdds;
    uint64_t*       vars;
    uint64_t        var_count;
};

extern "C" {
	// Export constants (wrapped)
	const sylvan::MTBDD w_mtbdd_true  = sylvan::mtbdd_true; 
	const sylvan::MTBDD w_mtbdd_false = sylvan::mtbdd_false;

	// Functions
	sylvan::MTBDD amaya_unite_mtbdds(sylvan::MTBDD m1, sylvan::MTBDD m2);

	sylvan::MTBDD amaya_project_variables_away(
			sylvan::MTBDD m, 
			uint32_t *variables,
			uint32_t var_count);

	State* amaya_mtbdd_get_transition_target(
			sylvan::MTBDD mtbdd, 
			uint8_t* cube, 
			uint32_t cube_size, 
			uint32_t* result_size);

	void amaya_print_dot(sylvan::MTBDD m, int32_t fd);

	sylvan::MTBDD amaya_mtbdd_build_single_terminal(
		uint8_t*  transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t  transition_symbols_count,
		uint32_t  variable_count,
		State*    destination_set,
		uint32_t  destination_set_size);

	sylvan::MTBDD* amaya_mtbdd_rename_states(
			sylvan::MTBDD* mtbdd_roots, 
			uint32_t 	root_count,
			State* 		state_name_pairs,  // [(old, new), (old, new), (old, new)] 
			uint32_t 	state_name_pairs_cnt);

	/**
	 * @param leaf_ptrs 	If not NULL will point to an array containing pointers 
	 * 						to transition destinations set of the leaves.
	 * @returns 			Contents of leaves serialized into array: [leafA1, leafA2, ..., leafAN, leafB1, ...]
	 */
	State* amaya_mtbdd_get_leaves(
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
			void* 	 leaf_tds, 
			State* 	 new_leaf_contents, 
			uint32_t contents_size);

	/**
	 * Calculates the post set from given MTBDD.
	 * @param m The transitions MTBDD for the state for which we calculate the post set.
	 * @param post_size (out) The size of the returned array.
	 * @returns array containing the post set.
	 */
	State* amaya_mtbdd_get_state_post(
			sylvan::MTBDD  m, 
			uint32_t* 	   post_size);

	/**
	 * Apply the padding closure to a state with transition function `left`, with successor's function `right`.
	 * @param left 			The transitions of state into which we want to propagate finishing symbols.
	 * @param right 		The transitions of a state that has some transitions leading to the final state.
	 * @param final_states 		Pointer to an array containing the final states of the automaton.
	 * @param final_states_cnt 	Number of states in the final_state array.
	 * @returns The newly created MTBDD. If the MTBDD has the same value as the left_dd, 
     *          then the pad_closure had no effect.
	 */
	sylvan::MTBDD amaya_mtbdd_do_pad_closure(
			State 			left_state,
			sylvan::MTBDD 	left_dd, 
			State 			right_state,
			sylvan::MTBDD 	right_dd);

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
			uint32_t* 	vars,	
			uint32_t 	var_count,
			uint32_t* 	symbols_cnt,
			State** 	dest_states,
			uint32_t** 	dest_states_cnt);

	/**
	 * Calculate the intersection of two mtbdds (Operation on leaves is set intersection).
	 * @param a 				Some transition MTBDD
	 * @param b 				Other transition MTBDD
	 * @param discovered_states 		Flat array of discovered macrostates and their generated int.
	 * 									[macrostate_left, macrostate_right, state_int ...]
	 * @param discovered_states_cnt 	Discovered_states count (size/3).
	 * @returns A MTBDD that contains only transitions that are present in both given MTBDDs.
	 */
	sylvan::MTBDD amaya_mtbdd_intersection(
			sylvan::MTBDD a, 
			sylvan::MTBDD b,
			State**  	discovered_states,         // OUT
			uint32_t*  	discovered_states_cnt);    // OUT

	
	/**
	 * Marks the beginning of intersection. Setups global std::map that holds
	 * information about which intersection macrostates (pairs) were mapped to which values.
	 */
	void amaya_begin_intersection(
			bool 		should_do_early_prunining, 
			State* 		prune_final_states, 
			uint32_t 	final_states_cnt);

	/**
	 * Updates the intersection state information by inserting provided macrostates (flat 2d array)
	 * and their corresponding mappings into intersection state.
	 *
	 * @param macrostates  			Flattened 2D array containing the intersection macrostates (always a pair)
	 * @param renamed_macrostates  	Renamed macrostates.
	 * @param cnt 					Count of the macrostates.
	 */
	void amaya_update_intersection_state(
			State* 	 macrostates, 
			State* 	 renamed_macrostates, 
			uint32_t macrostates_cnt);
	/**
	 * Marks the end of intersection and frees up resources.
	 */
	void amaya_end_intersection();

	void amaya_set_debugging(bool debug);

	void amaya_do_free(void *ptr);

	/** 
	 * Renames the macrostates that are contained withing the roots of mtbdds resulting 
	 * from the determinization procedure some automaton.
	 * @param roots 					The roots of the MTBDDs that were created during the determinization procedure.  
	 * @param root_cnt  				The number of given MTBDDs. 
	 * @param out_macrostates_sizes   	OUTPUT: The sized of the located macrostates.
	 * @param out_macrostates_cnt   		OUTPUT: The number of located macrostates.
	 * @param out_serialized_macrostates OUTPUT: The macrostates located serialized one after another.
	 * @returns The transformed mtbdds.
	 */
	sylvan::MTBDD* amaya_rename_macrostates_to_int(
			sylvan::MTBDD* 		roots, 							// MTBDDs resulting from determinization
			uint32_t 			root_cnt,						// Root count
			State 				start_numbering_macrostates_from,
			State**				out_serialized_macrostates,
			uint64_t**			out_macrostates_sizes,
			uint64_t*			out_macrostates_cnt);


	sylvan::MTBDD amaya_complete_mtbdd_with_trapstate(
			sylvan::MTBDD mtbdd,
			State 		trapstate,
			bool* 		had_effect);
	
	/**
	 * Walks the MTBDD building a set of reachable states encoded within the MTBDD. For
	 * every located state also notes the transition symbol via which can the state be reached.
	 * Every reachable state is presented in the result array max. 1 times (it is an array representation
	 * of a set).
	 *	
	 * @param mtbdd 			The mtbdd for which the reachable states will be retrieved. 
	 * @param variables			An array containing the indices of used variables.
	 * @param variable_cnt		The number of variables in the variable array.
	 * @param out_symbols		(OUT) Will contain transition symbols corresponting the returned states.
	 * @param transition_cnt	(OUT) Will contain the number of states located.
	 * @returns 				The array containing the located states.
	 */
	State* amaya_get_state_post_with_some_transition(
			sylvan::MTBDD mtbdd,
			uint32_t* 	variables,
			uint32_t 	variable_cnt,
			uint8_t** 	out_symbols,
			uint32_t* 	transition_cnt);

	sylvan::MTBDD* amaya_remove_states_from_transitions(
			sylvan::MTBDD* transition_roots,
			uint32_t 	transition_cnt,
			State* 		states_to_remove,
			uint32_t 	states_to_remove_cnt);

	/**
	 * Retrives all leaves that are present in the leaves of given MTBDDs.
	 * Debug method.
	 *
	 * @param mtbdds 			An array of MTBDDs whose leaves will be searched for states.
	 * @param mtbdd_cnt  		The number of MTBDDs passed in.
	 * @param out_states_cnt  	Output parameter - the number of found states.
	 *
	 * @returns An array of states present in the leaves of given MTBDDs.
	 */
	State* amaya_get_states_in_mtbdd_leaves(
		sylvan::MTBDD* mtbdds,
		uint32_t mtbdd_cnt,
		uint32_t* out_state_cnt);

	void amaya_begin_pad_closure(
        State new_final_state,
		State *final_states,
		uint32_t final_states_cnt);

    struct Serialized_DFA* amaya_minimize_hopcroft(struct Serialized_DFA* serialized_dfa);

	void amaya_end_pad_closure();

	void amaya_mtbdd_ref(sylvan::MTBDD mtbdd);
	void amaya_mtbdd_deref(sylvan::MTBDD mtbdd);
	void amaya_sylvan_gc();
	void amaya_sylvan_try_performing_gc();
    
	void amaya_sylvan_clear_cache();

	void shutdown_machinery();
	void init_machinery();

}

Transition_Destination_Set* _get_transition_target(
		sylvan::MTBDD root, 
		uint32_t current_variable,
		uint8_t* variable_assigments, 
		uint32_t var_count);

void collect_mtbdd_leaves(sylvan::MTBDD root, std::set<sylvan::MTBDD>& dest);

#endif
