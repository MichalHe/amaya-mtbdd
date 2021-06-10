#ifndef AMAYA_OPERATIONS_H
#define AMAYA_OPERATIONS_H

#include "base.hpp"
#include <vector>
#include <utility>
#include <map>
#include <unordered_set>

#include <sylvan.h>

#define AMAYA_EXISTS_OPERATION_ID 0x2000000
#define AMAYA_UNION_OP_ID 0x4000000
#define AMAYA_INTERSECTION_OP_ID 0x6000000

TASK_DECL_3(sylvan::MTBDD, set_intersection_op, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);

/**
 * The *abstraction* F definition.
 * Task returns a MTBDD. Task accepts two MTBDDs - left and right child (subtree) of a
 * node for variable <v> specified in variable set passed to abstract_apply with this operator.
 * The final int is a number of variables that are missing when reaching a leaf in MTBDD, but there is 
 * still <k> number of variables in the variable set passed to abstract_apply.
 */
TASK_DECL_3(sylvan::MTBDD, my_abstract_exists_op, sylvan::MTBDD, sylvan::MTBDD, int);
TASK_DECL_3(sylvan::MTBDD, pad_closure_op, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);
TASK_DECL_3(sylvan::MTBDD, set_union, sylvan::MTBDD *, sylvan::MTBDD *, uint64_t);

TASK_DECL_2(sylvan::MTBDD, remove_states_op, sylvan::MTBDD, uint64_t);
TASK_DECL_2(sylvan::MTBDD, complete_mtbdd_with_trapstate_op, sylvan::MTBDD, uint64_t);

typedef struct
{
	bool had_effect;
	State left_state;  // For debug purposes
	State right_state; // Actually used
	State *final_states;
	uint32_t final_states_cnt;
} Pad_Closure_Info;

typedef struct
{
	bool had_effect;
	uint32_t automaton_id;
	State trapstate;
} Complete_With_Trapstate_Op_Info;

typedef struct
{
	uint32_t automaton_id;
	std::vector<State> *discoveries; // Flat array of [metastate_left_part, metastate_right_part, state, ...]
} Intersection_Op_Info;

typedef struct
{
	bool should_do_early_prunining;
	std::unordered_set<State> *prune_final_states;
	std::map<std::pair<State, State>, State> *intersection_state_pairs_numbers;
} Intersection_State;

#endif
