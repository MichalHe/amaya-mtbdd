#include "wrapper.hpp"

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

#include "base.hpp"
#include "custom_leaf.hpp"
#include "operations.hpp"

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
extern uint64_t mtbdd_leaf_type_set;

extern void* 		REMOVE_STATES_OP_PARAM;
extern uint64_t 	REMOVE_STATES_OP_COUNTER;
extern void*		ADD_TRAPSTATE_OP_PARAM;
extern uint64_t 	ADD_TRAPSTATE_OP_COUNTER;
extern uint32_t 	CUR_PADDING_CLOSURE_ID;
extern Intersection_State* intersection_state;

extern uint64_t 				STATE_RENAME_OP_COUNTER;
extern State_Rename_Op_Info 	*STATE_RENAME_OP_PARAM;

extern uint64_t    TRANSFORM_METASTATES_TO_INTS_COUNTER;
extern Transform_Metastates_To_Ints_State *TRANSFORM_METASTATES_TO_INTS_STATE;

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
	sylvan_set_sizes(1LL << 26, 1LL << 26, 1LL << 20, 1LL << 20);
	//sylvan_set_sizes(1LL << 24, 1LL << 28, 1LL << 24, 1LL << 28);
    sylvan_init_package();
    sylvan_init_mtbdd();
	
	//sylvan_gc_hook_pregc(TASK(gc_start));

  // Initialize my own sylvan type.

    mtbdd_leaf_type_set = sylvan_mt_create_type();
    sylvan_mt_set_hash(mtbdd_leaf_type_set, set_leaf_hash);
    sylvan_mt_set_equals(mtbdd_leaf_type_set, set_leaf_equals);
    sylvan_mt_set_create(mtbdd_leaf_type_set, mk_set_leaf);
    sylvan_mt_set_destroy(mtbdd_leaf_type_set, destroy_set_leaf);
    sylvan_mt_set_to_str(mtbdd_leaf_type_set, set_leaf_to_str);
}

void shutdown_machinery() 
{
	sylvan_quit();
	lace_exit();
}

MTBDD amaya_mtbdd_build_single_terminal(
		uint32_t  automaton_id,
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

    Transition_Destination_Set* tds = new Transition_Destination_Set(automaton_id, leaf_state_set);
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
		uint32_t automaton_id, 
		State 	 trapstate,
		bool* 	 had_effect) 
{

	Complete_With_Trapstate_Op_Info op_info = {};
	op_info.trapstate = trapstate;
	op_info.automaton_id = automaton_id;
	op_info.had_effect = false;
	
	LACE_ME;
	ADD_TRAPSTATE_OP_PARAM = &op_info;
	MTBDD result = mtbdd_uapply(dd, TASK(complete_transition_with_trapstate_op), ADD_TRAPSTATE_OP_COUNTER);
	ADD_TRAPSTATE_OP_COUNTER++;

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
	STATE_RENAME_OP_COUNTER += 1;

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

MTBDD* amaya_rename_metastates_to_int(
		MTBDD* 		roots, 							// MTBDDs resulting from determinization
		uint32_t 	root_cnt,						// Root count
		State 		start_numbering_metastates_from,
		uint32_t 	resulting_automaton_id,
		State**		out_serialized_metastates,
		uint64_t**	out_metastates_sizes,
		uint64_t*	out_metastates_cnt)
{
	// Set up the tranform state
	TRANSFORM_METASTATES_TO_INTS_STATE = (Transform_Metastates_To_Ints_State *) malloc(sizeof(Transform_Metastates_To_Ints_State));
	TRANSFORM_METASTATES_TO_INTS_STATE->first_available_state_number = start_numbering_metastates_from;
	TRANSFORM_METASTATES_TO_INTS_STATE->metastates_cnt = 0;
	TRANSFORM_METASTATES_TO_INTS_STATE->metastates_sizes = new vector<uint64_t>();
	TRANSFORM_METASTATES_TO_INTS_STATE->serialized_metastates = new vector<State>();
	
	MTBDD *transformed_mtbdds = (MTBDD *) malloc(sizeof(MTBDD) * root_cnt);
	assert(transformed_mtbdds);

	LACE_ME;
	for (uint32_t i = 0; i < root_cnt; i++) {
		MTBDD transformed_mtbdd = mtbdd_uapply(roots[i], TASK(transform_metastates_to_ints_op), TRANSFORM_METASTATES_TO_INTS_COUNTER);
		transformed_mtbdds[i] = transformed_mtbdd;
		mtbdd_ref(transformed_mtbdd);
	}

	TRANSFORM_METASTATES_TO_INTS_COUNTER += 1;

	// Write values to the Python side.	
	State* serialized_metastates = (State *) malloc(sizeof(State) * TRANSFORM_METASTATES_TO_INTS_STATE->serialized_metastates->size());
	assert(serialized_metastates);
	for (uint64_t i = 0; i < TRANSFORM_METASTATES_TO_INTS_STATE->serialized_metastates->size(); i++) {
		serialized_metastates[i] = TRANSFORM_METASTATES_TO_INTS_STATE->serialized_metastates->at(i);
	}
	*out_serialized_metastates = serialized_metastates;
	
	uint64_t* metastate_sizes = (uint64_t*) malloc(sizeof(uint64_t) * TRANSFORM_METASTATES_TO_INTS_STATE->metastates_cnt);
	assert(metastate_sizes);
	for (uint64_t i = 0; i < TRANSFORM_METASTATES_TO_INTS_STATE->metastates_sizes->size(); i++) {
		metastate_sizes[i] = TRANSFORM_METASTATES_TO_INTS_STATE->metastates_sizes->at(i);
	}
	*out_metastates_sizes = metastate_sizes;

	*out_metastates_cnt = TRANSFORM_METASTATES_TO_INTS_STATE->metastates_cnt;

	free(TRANSFORM_METASTATES_TO_INTS_STATE);
	TRANSFORM_METASTATES_TO_INTS_STATE = NULL;
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


bool amaya_mtbdd_do_pad_closure(
		State 		left_state, 
		MTBDD 		left_dd,
		State 		right_state, 
		MTBDD 		right_dd, 
		State* 		final_states, 
		uint32_t 	final_states_cnt)
{
	Pad_Closure_Info pci = {0};
	pci.final_states = final_states;
	pci.final_states_cnt = final_states_cnt;
	pci.had_effect = false;

	pci.right_state = right_state;
	pci.left_state = left_state;
	

	LACE_ME;
	//sylvan_clear_cache();
	mtbdd_applyp(left_dd, 
				 right_dd,
				 (uint64_t) &pci, 
				 TASK(pad_closure_op), 
				 CUR_PADDING_CLOSURE_ID);

	// Every padding closure that is performed (even repeated) should be treated
	// as a unique operation --> so that the caching problems will not occur
	CUR_PADDING_CLOSURE_ID += 4;

	return pci.had_effect;
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

void amaya_mtbdd_change_automaton_id_for_leaves(
        MTBDD* roots,
        uint32_t root_cnt,
        uint32_t new_id)
{
	set<MTBDD> leaves {};
    for (uint32_t i = 0; i < root_cnt; i++) {
        collect_mtbdd_leaves(roots[i], leaves);
    }

    for (auto leaf : leaves) {
        auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf); 
        tds->automaton_id = new_id;
    }
}

MTBDD amaya_mtbdd_intersection(
        MTBDD 		a, 
        MTBDD 		b,
        uint32_t 	result_automaton_id, 
        State** 	discovered_states,         // OUT
        uint32_t*  	discovered_states_cnt)     // OUT
{
    Intersection_Op_Info intersect_info = {0};
    intersect_info.automaton_id = result_automaton_id;
    intersect_info.discoveries = new vector<State>();

	LACE_ME;
	MTBDD result = mtbdd_applyp(a, b, (uint64_t) &intersect_info, TASK(transitions_intersection_op), AMAYA_INTERSECTION_OP_ID);
    
    
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

MTBDD amaya_unite_mtbdds(MTBDD m1, MTBDD m2, uint32_t automaton_id) {
	LACE_ME;
	MTBDD u = mtbdd_applyp(m1, m2, (uint64_t) automaton_id, TASK(transitions_union_op), AMAYA_UNION_OP_ID);
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

void amaya_update_intersection_state(State* metastates, State* renamed_metastates, uint32_t cnt)
{
    for (uint32_t i = 0; i < cnt; i++) {
        const auto metastate = std::make_pair(metastates[2*i], metastates[2*i + 1]);
        intersection_state
			->intersection_state_pairs_numbers
			->insert(std::make_pair(metastate, renamed_metastates[i])); 
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
	
	LACE_ME;
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

