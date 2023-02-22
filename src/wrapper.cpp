#include "../include/wrapper.hpp"
#include "../include/base.hpp"
#include "../include/custom_leaf.hpp"
#include "../include/hopcroft_leaf.hpp"
#include "../include/operations.hpp"

#include <algorithm>
#include <sylvan.h>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <unistd.h>
#include <assert.h>
#include <unordered_set>
#include <lace.h>
#include <utility>
#include <iostream>
#include <cstring>

using namespace sylvan;

using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::map;
using std::pair;
using std::unordered_set;

static bool DEBUG_ON = false;

extern void* 		REMOVE_STATES_OP_PARAM;
extern uint64_t 	REMOVE_STATES_OP_COUNTER;
extern void*		ADD_TRAPSTATE_OP_PARAM;
extern uint64_t 	ADD_TRAPSTATE_OP_COUNTER;

extern Intersection_State* intersection_state;
extern uint64_t INTERSECTION_OP_COUNTER;

extern uint64_t 				STATE_RENAME_OP_COUNTER;
extern State_Rename_Op_Info 	*STATE_RENAME_OP_PARAM;

extern uint64_t    TRANSFORM_MACROSTATES_TO_INTS_COUNTER;
extern Transform_Macrostates_To_Ints_State *TRANSFORM_MACROSTATES_TO_INTS_STATE;

extern uint64_t 	PAD_CLOSURE_OP_COUNTER;
extern Pad_Closure_Info *PAD_CLOSURE_OP_STATE;



VOID_TASK_0(gc_start)
{
	fprintf(stderr, "Sylvan - Starting garbage collection.\n");
}

void init_machinery()
{
    int n_workers = 1;
    size_t dequeue_size = 10000000;

    lace_init(n_workers, dequeue_size);
    //lace_startup(program_stack_size, TASK(_main), NULL);

    // THIS SEEMS TO BE THE SECRET
    // When the TASK parameter (the middle one) is set to be NULL, it does not spawn a new thread
    // and instead uses the current thread for all tasks (makes it possible to call from python)
	const size_t stack_size = 1LL << 20;
    lace_startup(0, NULL, NULL);

	sylvan_set_limits(500LL*1024*1024, 3, 5); // Allocate 100MB
	//sylvan_set_sizes(1LL << 27, 1LL << 26, 1LL << 26, 1LL << 20);
	//sylvan_set_sizes(1LL << 24, 1LL << 28, 1LL << 24, 1LL << 28);
    sylvan_init_package();
    sylvan_init_mtbdd();

	//sylvan_gc_hook_pregc(TASK(gc_start));

    // Initialize leaf type for leaves containing sets - represents outgoing transition from a state
    mtbdd_leaf_type_set = sylvan_mt_create_type();
    sylvan_mt_set_hash(mtbdd_leaf_type_set, set_leaf_hash);
    sylvan_mt_set_equals(mtbdd_leaf_type_set, set_leaf_equals);
    sylvan_mt_set_create(mtbdd_leaf_type_set, mk_set_leaf);
    sylvan_mt_set_destroy(mtbdd_leaf_type_set, destroy_set_leaf);
    sylvan_mt_set_to_str(mtbdd_leaf_type_set, set_leaf_to_str);

    // Initialize leave type for efficient hopcroft minimization
    mtbdd_leaf_type_hopcroft = sylvan_mt_create_type();
    sylvan_mt_set_hash(mtbdd_leaf_type_hopcroft,    hopcroft_leaf_hash);
    sylvan_mt_set_equals(mtbdd_leaf_type_hopcroft,  hopcroft_leaf_equals);
    sylvan_mt_set_create(mtbdd_leaf_type_hopcroft,  hopcroft_leaf_create);
    sylvan_mt_set_destroy(mtbdd_leaf_type_hopcroft, hopcroft_leaf_destroy);
    sylvan_mt_set_to_str(mtbdd_leaf_type_hopcroft,  hopcroft_leaf_to_str);
}

void shutdown_machinery()
{
	sylvan_quit();
	lace_exit();
}

MTBDD amaya_mtbdd_build_single_terminal(
		uint8_t*  transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t  transition_symbols_count,
		uint32_t  variable_count,
		State* 	  destination_set,
		uint32_t  destination_set_size)
{
	// CALL(mtbdd_union_cube(mtbdd, variables, cube, terminal)
	// Cube (transitions symbols) have to be in format -- suppose the transition_symbols are ready
	// x_0 = 0  ===> low transition
	// x_0 = 1  ===> high transition
	// x_0 = 2  ===> any
	// Construct variables set containing all the variables.
	BDDSET variables = mtbdd_set_empty();
	for (uint32_t i=1; i <= variable_count; i++) {
		variables = mtbdd_set_add(variables, i); // Variables are numbered from 1
	}

	// Construct the destination set

	auto leaf_state_set = new std::set<State>();
	for (uint32_t i = 0; i < destination_set_size; i++) {
		leaf_state_set->insert(destination_set[i]);
	}

    Transition_Destination_Set* tds = new Transition_Destination_Set(leaf_state_set);
	MTBDD leaf = make_set_leaf(tds);

	// Construct the initial MTBDD, then add the rest of the symbols
	// signature: mtbdd_cube(MTBDD variables, uint8_t *cube, MTBDD terminal)
	// initial_cube = transition_symbols[0*variable_count + (0..variable_count)] (size: variable_count)
	MTBDD mtbdd = mtbdd_cube(variables, transition_symbols, leaf);

	LACE_ME;
	for (uint32_t i = 1; i < transition_symbols_count; i++) {
		// Cube for this iteration:
		//   i*variable_count + (0..variable_count)
		mtbdd = CALL(mtbdd_union_cube, mtbdd, variables, &transition_symbols[i*variable_count], leaf);
	}

	return mtbdd;
}

MTBDD amaya_complete_mtbdd_with_trapstate(
		MTBDD 	 dd,
		State 	 trapstate,
		bool* 	 had_effect)
{

	Complete_With_Trapstate_Op_Info op_info = {};
	op_info.trapstate = trapstate;
	op_info.had_effect = false;

	LACE_ME;
	ADD_TRAPSTATE_OP_PARAM = &op_info;
	MTBDD result = mtbdd_uapply(dd, TASK(complete_transition_with_trapstate_op), ADD_TRAPSTATE_OP_COUNTER);
	ADD_TRAPSTATE_OP_COUNTER += (1LL << 5);

	if (DEBUG_ON) printf("Had the complete with trapstate effect? %s\n", (op_info.had_effect ? "Yes": "No"));
	*had_effect = op_info.had_effect;
	return result;
}

Transition_Destination_Set* _get_transition_target(
		MTBDD root,
		uint32_t current_variable,
		uint8_t* variable_assigments,
		uint32_t var_count)
{
	int is_leaf = mtbdd_isleaf(root);
	if (var_count == 0) {
		if (is_leaf) {
			if (root == mtbdd_false) return NULL;
			auto tds = (Transition_Destination_Set*) mtbdd_getvalue(root);
			auto tds_copy = new Transition_Destination_Set(*tds);
			return tds_copy;
		} else {
			// We have exhausted all variables, yet we did not hit a leaf.
			// That means the cube did not contain all the variables.
			return NULL;  // Not enough information to decide.
		}
	}

	if (is_leaf) {
		if (root == mtbdd_false) return NULL;
		else {
			// Found a solution
			auto tds = (Transition_Destination_Set*) mtbdd_getvalue(root);
			auto tds_copy = new Transition_Destination_Set(*tds); // Return a copy (might modify elsewhere)
			return tds_copy;
		}
	}

	// If we got here then root is not a leaf and we still have some variables
	uint32_t root_var = mtbdd_getvar(root);

	if (root_var == current_variable) {
		if (*variable_assigments == 0) {
			// Go low
			return _get_transition_target(
					mtbdd_getlow(root),
					current_variable + 1,
					variable_assigments + 1,
					var_count - 1);
		}
		else if (*variable_assigments == 1) {
			// Go high
			return _get_transition_target(
					mtbdd_gethigh(root),
					current_variable + 1,
					variable_assigments + 1,
					var_count - 1);
		}
		else {
			// That means that the current assigment is don't care (2)
            // Then we need to explore both branches and return the results
			auto low_tds = _get_transition_target(
					mtbdd_getlow(root),
					current_variable + 1,
					variable_assigments + 1,
					var_count - 1);

			auto high_tds = _get_transition_target(
					mtbdd_gethigh(root),
					current_variable + 1,
					variable_assigments + 1,
					var_count - 1);

			// Now decide what to do with the results
			if (low_tds == NULL) return high_tds; // Might return NULL
			if (high_tds == NULL) return low_tds; // The low result is not NULL

            // Perform state merging.
			low_tds->destination_set->insert(
                    high_tds->destination_set->begin(),
                    high_tds->destination_set->end());
			delete high_tds; // Only low result will be propagated upwards
			return low_tds;
		}

	} else {
		// The current variable in the tree skipped a variable.
		// That means the current variable assignment did not matter.
		return _get_transition_target(
				root,
				current_variable + 1,   		// Move to the next variable
				variable_assigments + 1,		// --^
				var_count - 1);
	}
}

State* amaya_mtbdd_get_transition_target(
        MTBDD 		mtbdd,
        uint8_t* 	cube,
        uint32_t 	cube_size,
        uint32_t* 	result_size)
{

	auto search_result_tds = _get_transition_target(mtbdd, 1, cube, cube_size);
	if (search_result_tds == NULL)
	{
		*result_size = 0;
		return NULL;
	} else {
		const uint32_t rs = search_result_tds->destination_set->size();
		*result_size = rs;
		auto result_arr = (State*) malloc(sizeof(uint32_t) * rs);
		uint32_t i = 0;
		for (auto state : *search_result_tds->destination_set)
		{
			result_arr[i++] = state;
		}
		return result_arr;
	}
}

void amaya_do_free(void *ptr)
{
	free(ptr);
}


void collect_mtbdd_leaves(MTBDD root, std::set<MTBDD>& dest)
{
	LACE_ME;
	MTBDD support = mtbdd_support(root);
	uint32_t support_size = mtbdd_set_count(support);

	// Stores the tree path to the leaf
	uint8_t* arr = (uint8_t*) malloc(sizeof(uint8_t) * support_size);

	MTBDD leaf = mtbdd_enum_first(root, support, arr, NULL);

 	while (leaf != mtbdd_false)
	{
		dest.insert(leaf);
 	    leaf = mtbdd_enum_next(root, support, arr, NULL);
 	}

	free(arr);
}

MTBDD* amaya_mtbdd_rename_states(
		MTBDD* 		mtbdd_roots,
		uint32_t 	root_count,
		State* 		names, // [(old, new), (old, new), (old, new)]
		uint32_t 	name_count)
{
	State old_state, new_state;

	std::map<State, State> state_names_map;
	for (uint32_t i = 0; i < name_count; i++) {
		old_state = names[2*i];
		new_state = names[2*i + 1];
		state_names_map.insert(std::make_pair(old_state, new_state));
	}

	State_Rename_Op_Info op_info = {0};
	op_info.states_rename_map = &state_names_map;
	STATE_RENAME_OP_PARAM = &op_info;

	auto renamed_mtbdds = (MTBDD *) malloc(sizeof(MTBDD) * root_count);
	assert(renamed_mtbdds != NULL);

	LACE_ME;
	for (uint32_t i = 0; i < root_count; i++) {
		MTBDD result = mtbdd_uapply(mtbdd_roots[i], TASK(rename_states_op), STATE_RENAME_OP_COUNTER);
		mtbdd_ref(result);
		renamed_mtbdds[i] = result;
	}
	STATE_RENAME_OP_COUNTER += (1LL << 32);

	return renamed_mtbdds;
}


State* amaya_mtbdd_get_leaves(
		MTBDD root,
		uint32_t** leaf_sizes,	// OUT, Array containing the sizes of leaves inside dest
		uint32_t*  leaf_cnt,	// OUT, Number of leaves in the tree
        void***    leaf_ptrs)   // Pointer to an array of pointers to MTBDD leaves
{
	std::set<MTBDD> leaves {};
	collect_mtbdd_leaves(root, leaves);

	// This probably is not the most efficient way how to do it.
	// First compute the destination size (for malloc)
	uint32_t size_cnt = 0; // Counter for the total number of states.
	for (MTBDD leaf : leaves) {
		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		size_cnt += tds->destination_set->size();
	}

	// Do the allocations
	auto _leaf_sizes = (uint32_t*) malloc(sizeof(uint32_t) * leaves.size());  // One slot for each leaf
	auto _leaf_states = (State*) malloc(sizeof(State) * size_cnt);

	// Populate
	uint32_t state_i = 0;
	uint32_t size_i = 0;

	for (MTBDD leaf : leaves) {
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (int state : *tds->destination_set) {
			_leaf_states[state_i] = state;
			state_i++;
		}
		_leaf_sizes[size_i] = tds->destination_set->size();
		size_i++;
	}

    if (leaf_ptrs != NULL) {
        void** leaf_ptr_arr = (void**) malloc(sizeof(void*) * leaves.size());

        uint32_t i = 0;
        for (auto leaf : leaves) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
            leaf_ptr_arr[i] = (void*) tds;
            i++;
        }
        *leaf_ptrs = leaf_ptr_arr;
    }

	*leaf_sizes = _leaf_sizes;
	*leaf_cnt = (uint32_t) leaves.size();
	return _leaf_states;
}

void amaya_replace_leaf_contents_with(void *leaf_tds, State* new_contents, uint32_t contents_size)
{
	assert(false);
    auto tds = (Transition_Destination_Set*) leaf_tds;
    tds->destination_set->clear();
    for (uint32_t i = 0; i < contents_size; i++) {
        tds->destination_set->insert(new_contents[i]);
    }
}

MTBDD* amaya_rename_macrostates_to_int(
		MTBDD* 		roots, 							// MTBDDs resulting from determinization
		uint32_t 	root_cnt,						// Root count
		State 		start_numbering_macrostates_from,
		State**		out_serialized_macrostates,
		uint64_t**	out_macrostates_sizes,
		uint64_t*	out_macrostates_cnt)
{
	// Set up the tranform state
	TRANSFORM_MACROSTATES_TO_INTS_STATE = (Transform_Macrostates_To_Ints_State *) malloc(sizeof(Transform_Macrostates_To_Ints_State));
	TRANSFORM_MACROSTATES_TO_INTS_STATE->first_available_state_number = start_numbering_macrostates_from;
	TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_cnt = 0;
	TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_sizes = new vector<uint64_t>();
	TRANSFORM_MACROSTATES_TO_INTS_STATE->serialized_macrostates = new vector<State>();
	TRANSFORM_MACROSTATES_TO_INTS_STATE->alias_map = new std::map<std::set<State>, State>();

	MTBDD *transformed_mtbdds = (MTBDD *) malloc(sizeof(MTBDD) * root_cnt);
	assert(transformed_mtbdds);

	LACE_ME;
	for (uint32_t i = 0; i < root_cnt; i++) {
		MTBDD transformed_mtbdd = mtbdd_uapply(roots[i], TASK(transform_macrostates_to_ints_op), TRANSFORM_MACROSTATES_TO_INTS_COUNTER);
		transformed_mtbdds[i] = transformed_mtbdd;
		mtbdd_ref(transformed_mtbdd);
	}

	TRANSFORM_MACROSTATES_TO_INTS_COUNTER += 1;

	// Write values to the Python side.
	State* serialized_macrostates = (State *) malloc(sizeof(State) * TRANSFORM_MACROSTATES_TO_INTS_STATE->serialized_macrostates->size());
	assert(serialized_macrostates);
	for (uint64_t i = 0; i < TRANSFORM_MACROSTATES_TO_INTS_STATE->serialized_macrostates->size(); i++) {
		serialized_macrostates[i] = TRANSFORM_MACROSTATES_TO_INTS_STATE->serialized_macrostates->at(i);
	}
	*out_serialized_macrostates = serialized_macrostates;

	uint64_t* macrostate_sizes = (uint64_t*) malloc(sizeof(uint64_t) * TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_cnt);
	assert(macrostate_sizes);
	for (uint64_t i = 0; i < TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_sizes->size(); i++) {
		macrostate_sizes[i] = TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_sizes->at(i);
	}
	*out_macrostates_sizes = macrostate_sizes;

	*out_macrostates_cnt = TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_cnt;

	delete TRANSFORM_MACROSTATES_TO_INTS_STATE->macrostates_sizes;
	delete TRANSFORM_MACROSTATES_TO_INTS_STATE->serialized_macrostates;
	delete TRANSFORM_MACROSTATES_TO_INTS_STATE->alias_map;
	free(TRANSFORM_MACROSTATES_TO_INTS_STATE);
	TRANSFORM_MACROSTATES_TO_INTS_STATE = NULL;
	return transformed_mtbdds;
}


void amaya_print_dot(MTBDD m, int32_t fd)
{
	FILE* file = fdopen(fd, "w");
	mtbdd_fprintdot(file, m);
	fflush(file);
	// fclose(file);  The file will be closed in python
}

MTBDD amaya_project_variables_away(MTBDD m, uint32_t *variables, uint32_t var_count)
{
	// Construct the variables set.
	BDDSET var_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < var_count; i++) {
		var_set = mtbdd_set_add(var_set, variables[i]);
	}

	// Do the projection itself.
	LACE_ME;
	MTBDD result = mtbdd_abstract(m, var_set, TASK(project_variable_away_abstract_op));

	//mtbdd_fprintdot(stdout, result);
	return result;
}

State* amaya_mtbdd_get_state_post(MTBDD dd, uint32_t *post_size)
{
	set<MTBDD> leaves;
	collect_mtbdd_leaves(dd, leaves);

	set<State> state_post_set;

	for (auto leaf: leaves)	{
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		state_post_set.insert(
                tds->destination_set->begin(),
                tds->destination_set->end());
	}

	auto result = (State*) malloc(sizeof(State) * state_post_set.size());
	uint32_t i = 0;

	for (int state: state_post_set)	{
		result[i++] = state;
	}

	*post_size = state_post_set.size();
	return result;
}


void amaya_begin_pad_closure(
        State new_final_state,
		State *final_states,
		uint32_t final_states_cnt)
{
	PAD_CLOSURE_OP_STATE = (Pad_Closure_Info*) malloc(sizeof(Pad_Closure_Info));
	assert(PAD_CLOSURE_OP_STATE != NULL);

    // Make a copy of the final states, because they will be deallocated after
    // the Python function returns. The +1 is for the new_final_state.
	State* final_states_cpy = (State*) malloc(sizeof(State) * (final_states_cnt + 1));
	assert(final_states_cpy != NULL);

	std::memcpy(final_states_cpy, final_states, (sizeof(State) * final_states_cnt));

    // Treat new_final_state as if it was already present in the automaton, so that
    // no reallocs will need to happen during the pad closure
    final_states_cpy[final_states_cnt] = new_final_state;

    PAD_CLOSURE_OP_STATE->new_final_state  = new_final_state;
	PAD_CLOSURE_OP_STATE->final_states     = final_states_cpy;
	PAD_CLOSURE_OP_STATE->final_states_cnt = final_states_cnt + 1; // Include the new one aswell

    // TODO: Is the cache even operational?
	PAD_CLOSURE_OP_STATE->operation_id_cache = new std::unordered_map<State, std::pair<MTBDD, uint64_t>>();
	PAD_CLOSURE_OP_STATE->first_available_r_cache_id = (1LL << 33);
}

void amaya_end_pad_closure()
{
	free(PAD_CLOSURE_OP_STATE->final_states);
	delete PAD_CLOSURE_OP_STATE->operation_id_cache;
	free(PAD_CLOSURE_OP_STATE);
	PAD_CLOSURE_OP_STATE = NULL;

	PAD_CLOSURE_OP_COUNTER += (1LL << 38);
}


MTBDD amaya_mtbdd_do_pad_closure(
		State 		left_state,
		MTBDD 		left_dd,
		State 		right_state,
		MTBDD 		right_dd)
{
	// @Note: When performing pad closure the finality of states does not change, therefore, the final states are constant.
	// 		  This means that the operation really depends only on the provided mtbdds, and not so much on the param.
	// 		  The only *problem* is that the same mtbdds (or their parts) might be created in different automata, but the
	// 		  final states are different. ---> Provide transaction like interface.

	PAD_CLOSURE_OP_STATE->right_state = right_state;
	PAD_CLOSURE_OP_STATE->left_state  = left_state;

	// When we use the same R mtbdd **with the same r state** but with different L mtbbds we might reuse cache
	// This means that operation parameters which control the cache selection process can be chosen in such a way
	// that every time we encounter the same R mtbdd the result from previous applications sharing the same R will
	// be reused.
	auto r_state_op_cache_it = PAD_CLOSURE_OP_STATE->operation_id_cache->find(right_state);
	uint64_t r_op_cache_id;
	if (r_state_op_cache_it != PAD_CLOSURE_OP_STATE->operation_id_cache->end()) {
		// We have some cache for this r-state, check whether the r-mtbdd matches the stored one
		// IF not, that means the mtbdd for r-state was already modified during pad closure (at a different time)
		// And therefore new operation chache needs to be "assigned"
		auto r_cache_entry = r_state_op_cache_it->second;

		if (r_cache_entry.first != right_dd) {
			// The cache entry is not valid anymore
			r_op_cache_id = PAD_CLOSURE_OP_STATE->first_available_r_cache_id;

			PAD_CLOSURE_OP_STATE->operation_id_cache->insert(
					std::make_pair(
						right_state,
						std::make_pair(right_dd, r_op_cache_id)));

			PAD_CLOSURE_OP_STATE->first_available_r_cache_id += (1LL << 10);
		} else {
			r_op_cache_id = r_cache_entry.second;
		}
	} else {
			// There is no cache entry, generate new one.
			r_op_cache_id = PAD_CLOSURE_OP_STATE->first_available_r_cache_id;

			PAD_CLOSURE_OP_STATE->operation_id_cache->insert(
					std::make_pair(
						right_state,
						std::make_pair(right_dd, r_op_cache_id)));

			PAD_CLOSURE_OP_STATE->first_available_r_cache_id += (1LL << 10);
	}

	LACE_ME;
	//sylvan_clear_cache();

	PAD_CLOSURE_OP_COUNTER += (1LL << 38);

	MTBDD result =  mtbdd_applyp(
			left_dd,
			right_dd,
			1LL, 			         // Allows using caches for the same R-state + R-MTBDD results
			TASK(pad_closure_op),
			PAD_CLOSURE_OP_COUNTER); // Identifies the overall padding closure being performed,

	return result;
}

uint8_t* amaya_mtbdd_get_transitions(
		MTBDD 		root, 	// The mtbdd
		uint32_t* 	vars,
		uint32_t 	var_count,
		uint32_t* 	symbols_cnt,
		State** 	dest_states,
		uint32_t** 	dest_states_sizes)
{
	LACE_ME;

	// Construct MTBDD set from given vars
	MTBDD variable_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < var_count; i++) {
		variable_set= mtbdd_set_add(variable_set, vars[i]);
	}

	// Stores the tree path to the leaf
	auto arr = (uint8_t*) malloc(sizeof(uint8_t) * var_count);

	MTBDD leaf = mtbdd_enum_first(root, variable_set, arr, NULL);

	vector<State> dest_states_vec;
	vector<uint8_t> symbols;
	vector<uint32_t> state_sizes;
 	while (leaf != mtbdd_false)
	{
		// Add the destination states to the oveall vector
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (int state: *tds->destination_set) {
			dest_states_vec.push_back(state);
		}

		for (uint32_t i = 0; i < var_count; i++) {
			symbols.push_back(arr[i]);
		}

		state_sizes.push_back(tds->destination_set->size());

 	    leaf = mtbdd_enum_next(root, variable_set, arr, NULL);
 	}

	// Create output arrays
	auto symbols_out = (uint8_t*) malloc(sizeof(uint8_t) * symbols.size());
	for (uint32_t i = 0; i < symbols.size(); i++) symbols_out[i] = symbols.at(i);

	auto dest_states_sizes_out = (uint32_t*) malloc(sizeof(uint32_t) * state_sizes.size());
	for (uint32_t i = 0; i < state_sizes.size(); i++) dest_states_sizes_out[i] = state_sizes.at(i);

	auto  dest_states_out = (State*) malloc(sizeof(State) * dest_states_vec.size());
	for (uint32_t i = 0; i < dest_states_vec.size(); i++) dest_states_out [i] = dest_states_vec.at(i);

	// Do the return arrays assignment
	*dest_states = dest_states_out;
	*dest_states_sizes = dest_states_sizes_out;

	*symbols_cnt = state_sizes.size(); // For every destination size, there exists 1 symbol leading to it.
	return symbols_out;
}


MTBDD amaya_mtbdd_intersection(
        MTBDD 		a,
        MTBDD 		b,
        State** 	discovered_states,         // OUT
        uint32_t*  	discovered_states_cnt)     // OUT
{
    Intersection_Op_Info intersect_info = {0};
    intersect_info.discoveries = new vector<State>();

	LACE_ME;
	MTBDD result = mtbdd_applyp(a, b, (uint64_t) &intersect_info, TASK(transitions_intersection_op), INTERSECTION_OP_COUNTER);

    if (!intersect_info.discoveries->empty()) {
        auto discovered_states_arr = (State*) malloc(sizeof(State) * intersect_info.discoveries->size());
        for (uint32_t i = 0; i < intersect_info.discoveries->size(); i++) {
            discovered_states_arr[i] = intersect_info.discoveries->at(i);
        }

        // Send the new discoveries to the python side.
        *discovered_states = discovered_states_arr;
        *discovered_states_cnt = intersect_info.discoveries->size() / 3;
    } else {
        // Nothing new was discovered, do not malloc
        *discovered_states = NULL;
        *discovered_states_cnt = 0;
    }

    delete intersect_info.discoveries;

    return result;
}

MTBDD amaya_unite_mtbdds(MTBDD m1, MTBDD m2) {
	LACE_ME;
	MTBDD u = mtbdd_applyp(m1, m2, (uint64_t) 0, TASK(transitions_union_op), AMAYA_UNION_OP_ID);
	return u;
}

void amaya_begin_intersection(
		bool prune_state_pairs_with_one_final,
		State* prune_final_states,
		uint32_t final_states_cnt)
{
    intersection_state = (Intersection_State*) malloc(sizeof(Intersection_State));
	intersection_state->intersection_state_pairs_numbers = new map<pair<State, State>, State>();

	if (prune_state_pairs_with_one_final) {
		intersection_state->should_do_early_prunining = true;
		intersection_state->prune_final_states = new unordered_set<State>();
		for (uint32_t i = 0; i < final_states_cnt; i++) {
			intersection_state->prune_final_states->insert(prune_final_states[i]);
		}
	} else {
		intersection_state->should_do_early_prunining = false;
		intersection_state->prune_final_states = NULL;
	}
}

void amaya_update_intersection_state(State* macrostates, State* renamed_macrostates, uint32_t cnt)
{
    for (uint32_t i = 0; i < cnt; i++) {
        const auto macrostate = std::make_pair(macrostates[2*i], macrostates[2*i + 1]);
        intersection_state
			->intersection_state_pairs_numbers
			->insert(std::make_pair(macrostate, renamed_macrostates[i]));
    }
}

void amaya_end_intersection()
{
	delete intersection_state->intersection_state_pairs_numbers;
	if (intersection_state->should_do_early_prunining) {
		delete intersection_state->prune_final_states;
	}
    free(intersection_state);
    intersection_state = NULL;

	INTERSECTION_OP_COUNTER += (5LL << 31);
	//sylvan_gc();
}

void amaya_set_debugging(bool debug) {
	DEBUG_ON = debug;
}

State* amaya_get_state_post_with_some_transition(
		MTBDD mtbdd,
		uint32_t* variables,
		uint32_t variable_cnt,
		uint8_t** out_symbols,
		uint32_t* transition_cnt)
{
	LACE_ME;

	// Create a MTBDD Cube containing all requred variables.
	MTBDD variable_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < variable_cnt; i++) {
		variable_set = mtbdd_set_add(variable_set, variables[i]);
	}

	// Stores the tree path to the leaf
	uint8_t* arr = (uint8_t*) malloc(sizeof(uint8_t) * variable_cnt);

	MTBDD leaf = mtbdd_enum_first(mtbdd, variable_set, arr, NULL);

	vector<State> 	reachable_states;
	vector<uint8_t> transition_symbols;

 	while (leaf != mtbdd_false) {
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (auto state : *tds->destination_set) {

			// @Optimize: We use linear search over vector - there should be a relatively small number of states
			// 			  reachable from everystate, therefore it should be faster to use ordinary vector (chache lines)
			bool state_already_located = false;
			for (auto already_located_state : reachable_states) {
				if (already_located_state == state) {
					state_already_located = true;
					break;
				}
			}

			if (!state_already_located) {
				reachable_states.push_back(state);
				// Copy the transition symbol so that the Python side will have the information available
				// (this method is used in DFS)

				for (uint32_t i = 0; i < variable_cnt; i++) {
					transition_symbols.push_back(arr[i]);
				}
			}
		}

 	    leaf = mtbdd_enum_next(mtbdd, variable_set, arr, NULL);
 	}

	free(arr);

	auto _out_symbols = (uint8_t*) malloc(sizeof(uint8_t) * transition_symbols.size());
	for (uint32_t i = 0; i < transition_symbols.size(); i++) {
		_out_symbols[i] = transition_symbols.at(i);
	}

	*out_symbols = _out_symbols;
	*transition_cnt = (uint32_t) reachable_states.size();

	auto _reachable_states = (State*) malloc(sizeof(State) * reachable_states.size());
	for (uint32_t i = 0; i < reachable_states.size(); i++) _reachable_states[i] = reachable_states.at(i);

	return _reachable_states;
}


MTBDD* amaya_remove_states_from_transitions(
	MTBDD* 		transition_roots,
	uint32_t 	transition_cnt,
	State* 		states_to_remove,
	uint32_t 	states_to_remove_cnt)
{

	auto mtbdds_after_removal = (MTBDD*) malloc(sizeof(MTBDD) * transition_cnt);
	LACE_ME;

	set<State> states_to_remove_set;
	for (uint32_t i = 0; i < states_to_remove_cnt; i++) {
		states_to_remove_set.insert(states_to_remove[i]);
	}

	REMOVE_STATES_OP_PARAM = &states_to_remove_set;
	for (uint32_t i = 0; i < transition_cnt; i++) {
		MTBDD result_mtbdd = mtbdd_uapply(transition_roots[i], TASK(remove_states_op), REMOVE_STATES_OP_COUNTER);
		mtbdd_ref(result_mtbdd);
		mtbdds_after_removal[i] = result_mtbdd;
	}

	REMOVE_STATES_OP_COUNTER++;

	return mtbdds_after_removal;
}

State* amaya_get_states_in_mtbdd_leaves(
	MTBDD* mtbdds,
	uint32_t mtbdd_cnt,
	uint32_t* out_state_cnt)
{
	// Garther a set of all unique leaves present in the given MTBDDs
	std::set<MTBDD> mtbdd_leaves;
	for (uint32_t i = 0; i < mtbdd_cnt; i++) {
		collect_mtbdd_leaves(mtbdds[i], mtbdd_leaves);
	}

	// Gather unique states in the previously extraced leaves
	std::set<State> states;
	for (auto mtbdd_leaf : mtbdd_leaves) {
		auto leaf_tds = (Transition_Destination_Set*) mtbdd_getvalue(mtbdd_leaf);
		for (auto state : *leaf_tds->destination_set) {
			states.insert(state);
		}
	}

	State* out_states = (State*) malloc(states.size() * sizeof(State*));

	uint32_t i = 0;
	for (auto state : states) {
		out_states[i++] = state;
	}

	*out_state_cnt = states.size();

	return out_states;
}


void amaya_mtbdd_ref(MTBDD dd)
{
	mtbdd_ref(dd);
}

void amaya_mtbdd_deref(MTBDD dd)
{
	mtbdd_deref(dd);
}

void amaya_sylvan_gc()
{
	LACE_ME;
	sylvan_gc();
}

void amaya_sylvan_try_performing_gc()
{
	LACE_ME;
	sylvan_gc_test();
}

void amaya_sylvan_clear_cache()
{
	LACE_ME;
	sylvan_clear_cache();
}

struct Serialized_DFA* amaya_minimize_hopcroft(struct Serialized_DFA* serialized_dfa) 
{
    // Deserialize the DFA
    struct NFA dfa = {
        .states = {},
        .final_states = {},
        .initial_states = {serialized_dfa->initial_state},
        .transitions = {},
        .vars = sylvan::mtbdd_set_empty(),
        .var_count = serialized_dfa->var_count
    };

    for (uint64_t i = 0; i < serialized_dfa->state_count; i++){
        dfa.states.insert(serialized_dfa->states[i]);
        dfa.transitions[serialized_dfa->states[i]] = serialized_dfa->mtbdds[i];
    }

    for (uint64_t i = 0; i < serialized_dfa->final_state_count; i++) dfa.final_states.insert(serialized_dfa->final_states[i]);

    for (uint64_t i = 0; i < serialized_dfa->var_count; i++)
        dfa.vars = sylvan::mtbdd_set_add(dfa.vars, serialized_dfa->vars[i]);
    
    struct NFA minimized_dfa = minimize_hopcroft(dfa);
    
    struct Serialized_DFA* output_dfa = (struct Serialized_DFA*) malloc(sizeof(struct Serialized_DFA));
    assert(output_dfa != nullptr);

    set<State>::iterator state_it;
    uint64_t i = 0;
    
    output_dfa->states = (State*) malloc(sizeof(State) * minimized_dfa.states.size()); 
    assert(output_dfa->states != nullptr);
    for (state_it = minimized_dfa.states.begin(), i = 0;
         state_it != minimized_dfa.states.end();
         state_it++, i++) {
        output_dfa->states[i] = *state_it;
    }
    output_dfa->state_count = minimized_dfa.states.size(); 

    output_dfa->final_states = (State*) malloc(sizeof(State) * minimized_dfa.final_states.size()); 
    assert(output_dfa->final_states != nullptr);
    for (state_it = minimized_dfa.final_states.begin(), i = 0;
         state_it != minimized_dfa.final_states.end();
         state_it++, i++) {
        output_dfa->final_states[i] = *state_it;
    }
    output_dfa->final_state_count = minimized_dfa.final_states.size(); 

    output_dfa->initial_state = *minimized_dfa.initial_states.begin();
    
    output_dfa->mtbdds = (sylvan::MTBDD*) malloc(sizeof(sylvan::MTBDD) * minimized_dfa.states.size());

    i = 0;
    for (State state: minimized_dfa.states) {
        output_dfa->mtbdds[i] = minimized_dfa.transitions[state];
        i++;
    }

    output_dfa->vars = serialized_dfa->vars;
    output_dfa->var_count = serialized_dfa->var_count;

    return output_dfa;
}

