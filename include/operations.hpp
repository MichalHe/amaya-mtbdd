#ifndef AMAYA_OPERATIONS_H
#define AMAYA_OPERATIONS_H

#include "base.hpp"
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <sylvan.h>

#define AMAYA_EXISTS_OPERATION_ID 0x2000000
#define AMAYA_UNION_OP_ID 0x4000000

TASK_DECL_3(sylvan::MTBDD, transitions_intersection_op, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);

/**
 * The *abstraction* F definition.
 * Task returns a MTBDD. Task accepts two MTBDDs - left and right child (subtree) of a
 * node for variable <v> specified in variable set passed to abstract_apply with this operator.
 * The final int is a number of variables that are missing when reaching a leaf in MTBDD, but there is
 * still <k> number of variables in the variable set passed to abstract_apply.
 */
TASK_DECL_3(sylvan::MTBDD, project_variable_away_abstract_op, sylvan::MTBDD, sylvan::MTBDD, int);
TASK_DECL_3(sylvan::MTBDD, pad_closure_op, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);
TASK_DECL_3(sylvan::MTBDD, transitions_union_op, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);

TASK_DECL_2(sylvan::MTBDD, remove_states_op, sylvan::MTBDD, uint64_t);
TASK_DECL_2(sylvan::MTBDD, complete_transition_with_trapstate_op, sylvan::MTBDD, uint64_t);
TASK_DECL_2(sylvan::MTBDD, rename_states_op, sylvan::MTBDD, uint64_t);
TASK_DECL_2(sylvan::MTBDD, transform_macrostates_to_ints_op, sylvan::MTBDD, uint64_t);

typedef struct
{
    State new_final_state;      // Final state to be added if the saturation property is broken
	State left_state;           // For debug purposes
	State right_state;          // Actually used

	std::unordered_map<State, std::pair<sylvan::MTBDD, uint64_t>> *operation_id_cache;
	uint64_t first_available_r_cache_id;
	State *final_states;
	uint32_t final_states_cnt;
} Pad_Closure_Info;

typedef struct
{
	bool had_effect;
	State trapstate;
} Complete_With_Trapstate_Op_Info;

typedef struct
{
	std::vector<State> *discoveries; // Flat array of [macrostate_left_part, macrostate_right_part, state, ...]
} Intersection_Op_Info;

typedef struct
{
	bool should_do_early_prunining;
	std::unordered_set<State> *prune_final_states;
	std::map<std::pair<State, State>, State> *intersection_state_pairs_numbers;
} Intersection_State;

typedef struct
{
	std::map<State, State> *states_rename_map;
} State_Rename_Op_Info;

typedef struct
{
	std::map<std::set<State>, State>* alias_map;
	std::vector<State>* serialized_macrostates;
	std::vector<uint64_t>* macrostates_sizes;
	uint64_t macrostates_cnt;
	State first_available_state_number;
} Transform_Macrostates_To_Ints_State;

struct NFA minimize_hopcroft(struct NFA& nfa);


// NFA operations
void nfa_print_transitions(struct NFA& nfa);

std::pair<std::set<State>, std::set<State>>
fragment_dest_states_using_partition(const std::set<State>& dest_states, const std::set<State>& partition);

template<typename C>
std::string states_to_str(const C& states);

#endif
