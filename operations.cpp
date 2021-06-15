#include "operations.hpp"
#include "custom_leaf.hpp"
#include "base.hpp"

#include <sylvan.h>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <sylvan_bdd.h>
#include <sylvan_mt.h>
#include <sylvan_int.h>
#include <sylvan_mtbdd_int.h>
#include <sylvan_stats.h>
#include <lace.h>

#include <unistd.h>
#include <assert.h>

#include <algorithm>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <vector>
#include <map>

using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::map;
using std::pair;
using std::unordered_set;

using sylvan::MTBDD;
using sylvan::mtbdd_makeleaf;
using sylvan::mtbdd_makenode;
using sylvan::mtbdd_false;
using sylvan::mtbdd_isleaf;
using sylvan::mtbdd_getvalue;
using sylvan::mtbdd_invalid;
using sylvan::mtbdd_applyp_CALL;

extern uint64_t mtbdd_leaf_type_set;

/**
 * Global variables:
 */
void*		REMOVE_STATES_OP_PARAM = NULL;
uint64_t 	REMOVE_STATES_OP_COUNTER = 0;
void*		ADD_TRAPSTATE_OP_PARAM = NULL;
uint32_t 	CUR_PADDING_CLOSURE_ID = 64;

Intersection_State* intersection_state = NULL;
State_Rename_Op_Info *STATE_RENAME_OP_PARAM = NULL;

uint64_t 	ADD_TRAPSTATE_OP_COUNTER = 	(1LL << 32);
uint64_t 	STATE_RENAME_OP_COUNTER = (1LL << 33);

/**
 * Performs an intersection of two given transitions. When performing an intersection the states tuples
 * that get created must be renamed to an integer right away (to be consistent with the overall design -
 * states are always integers). This requires that the information about mapping metastates to corresponding
 * integers will be persistent in between the intersection of individual transition MTBDDs. This is solved
 * via intersection_state - a global variable.
 */
TASK_IMPL_3(MTBDD, transitions_intersection_op, MTBDD *, pa, MTBDD *, pb, uint64_t, param)
{
    MTBDD a = *pa, b = *pb;
    // Intersection with an empty set (mtbdd_false) is empty set (mtbdd_false)
    if (a == mtbdd_false || b == mtbdd_false)
    {
        return mtbdd_false;
    }

    // If both are leaves calculate intersection
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b))
    {
        auto intersect_info = (Intersection_Op_Info *)param;

        auto &tds_a = *((Transition_Destination_Set *)mtbdd_getvalue(a));
        auto &tds_b = *((Transition_Destination_Set *)mtbdd_getvalue(b));

        auto left_states = tds_a.destination_set;
        auto right_states = tds_b.destination_set;

        if (left_states->empty() || right_states->empty())
        {
            return mtbdd_false;
        }

        // Calculate cross product
		std::pair<State, State> metastate;
        State state;

        auto intersection_leaf_states = new std::set<State>();

        auto already_discovered_intersection_states = intersection_state->intersection_state_pairs_numbers;

        for (auto left_state : *left_states)
        {
            for (auto right_state : *right_states)
            {
                metastate = std::make_pair(left_state, right_state);
                auto pos = already_discovered_intersection_states->find(metastate);
                bool contains_metastate = (pos != already_discovered_intersection_states->end());
                if (contains_metastate)
                {
                    state = pos->second;
                }
                else
                {
                    // We have discovered a new state.
                    // Check (if early pruning is on) whether the state should be pruned.
                    if (intersection_state->should_do_early_prunining)
                    {
                        const bool is_left_in_pruned = (intersection_state->prune_final_states->find(left_state) != intersection_state->prune_final_states->end());
                        const bool is_right_in_pruned = (intersection_state->prune_final_states->find(right_state) != intersection_state->prune_final_states->end());

                        // Pruning is performed only when exactly one of the states is in the pruned set
                        if (is_left_in_pruned != is_right_in_pruned)
                        {
                            // The state should be pruned
                            continue;
                        }
                    }

                    // Update the global intersection state, so in the future every such
                    // state will get the same integer
                    state = already_discovered_intersection_states->size();
                    already_discovered_intersection_states->insert(std::make_pair(metastate, state));

                    // Update the discoveries local for this intersection, so that the Python
                    // side knows that we have discovered some new metastate and what number we assigned
					// to it.
                    intersect_info->discoveries->push_back(left_state);
                    intersect_info->discoveries->push_back(right_state);
                    intersect_info->discoveries->push_back(state);
                }

                intersection_leaf_states->insert(state);
            }
        }

        auto intersection_tds = new Transition_Destination_Set(intersect_info->automaton_id, intersection_leaf_states);
        MTBDD intersection_leaf = make_set_leaf(intersection_tds);
        return intersection_leaf;
    }

    // TODO: Perform pointer swap here, so the cache would be utilized
    // (commutative).
    // TODO: Utilize operation cache.
    return mtbdd_invalid;
}

/**
 * Unites the two transition MTBDDs.
 */
TASK_IMPL_3(MTBDD, transitions_union_op, MTBDD *, pa, MTBDD *, pb, uint64_t, param)
{
    MTBDD a = *pa, b = *pb;
    if (a == mtbdd_false && b == mtbdd_false)
        return mtbdd_false;

    // When one leaf is empty set (false),
    // the algorithm should return all reachable states for the other one
    if (a == mtbdd_false)
    {
        return b;
    }

    if (b == mtbdd_false)
    {
        return a;
    }
    // If both are leaves, we calculate union
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b))
    {
        auto &tds_a = *((Transition_Destination_Set *)mtbdd_getvalue(a));
        auto &tds_b = *((Transition_Destination_Set *)mtbdd_getvalue(b));

        // If the passed parameter is missing (-1) that means that we should derive it
        // from the leaves (TDS). This can happen when padding closure is performed, as
        // the leaf automaton id is not modified.
        if (param == -1)
        {
            assert(tds_a.automaton_id == tds_b.automaton_id);
            param = (uint32_t)tds_a.automaton_id;
        }

        std::set<State> *union_set = new std::set<State>();
        std::set_union(
            tds_a.destination_set->begin(), tds_a.destination_set->end(),
            tds_b.destination_set->begin(), tds_b.destination_set->end(),
            std::inserter(*union_set, union_set->begin()));

        Transition_Destination_Set *union_tds = new Transition_Destination_Set((uint32_t)param, union_set);

        MTBDD union_leaf = make_set_leaf(union_tds); // Wrap the TDS with a MTBDD leaf.
        return union_leaf;
    }

	return sylvan::mtbdd_invalid;
}

/**
 * ABOUT ABSTRACTIONS:
 * Abstraction operation defines what should happen when a variable that should be projected
 * away should happen with the subtrees (one for high, one for low) that should have 
 * the previous variable as a parent.
 *
 * THIS OPERATOR:
 * When projecting away an alphabet variable we perform union of the two MTBDD trees underneath.
 *
 */
TASK_IMPL_3(MTBDD, project_variable_away_abstract_op, MTBDD, a, MTBDD, b, int, k)
{
    MTBDD u = mtbdd_applyp(a, b, (uint64_t)-1, TASK(transitions_union_op), AMAYA_EXISTS_OPERATION_ID);
    return u;
}

/**
 * NOTES:
 * The param is actually not used, since the (const) param would cause the sylvan cache to return
 * a result from some previous remove-states-op application.
 */
TASK_IMPL_2(MTBDD, remove_states_op, MTBDD, dd, uint64_t, param) {
	if (dd == mtbdd_false) return mtbdd_false;

	if (mtbdd_isleaf(dd)) {
		// Param is ignored, instead use the global variable REMOVE_STATES_OP_PARAM
		(void) param;

		auto states_to_remove = (set<State>*) (REMOVE_STATES_OP_PARAM);
		auto tds = (Transition_Destination_Set *) mtbdd_getvalue(dd);

		auto new_tds = new Transition_Destination_Set(*tds); // Make leaf value copy.
		
		for (auto state : *tds->destination_set) {
			bool should_be_removed = states_to_remove->find(state) != states_to_remove->end();	
			if (should_be_removed) new_tds->destination_set->erase(state);
		}

		if (new_tds->destination_set->empty()) {
			delete new_tds;
			return mtbdd_false;
		}
		
		return make_set_leaf(new_tds);
	}

	return mtbdd_invalid;
} 

/**
 * Completes the given transition MTBDD with trapstate.
 *
 * NOTE:
 * The param is not used - sylvan cachce problems, see remove_states_op. Instead the information
 * needed for the operation (like automaton ID and trapstate value) is passed via global variable
 * ADD_TRAPSTATE_OP_PARAM. This should be used with combination with ADD_TRAPSTATE_OP_COUNTER to 
 * achieve cache utilization only for the the current completition.
 */
TASK_IMPL_2(MTBDD, complete_transition_with_trapstate_op, MTBDD, dd, uint64_t, param)
{
	// param is used just to avoid sylvan miss-cachces
	if (dd == mtbdd_false) {
		(void) param;
		auto op_info = (Complete_With_Trapstate_Op_Info*) ADD_TRAPSTATE_OP_PARAM;

		Transition_Destination_Set* tds = new Transition_Destination_Set();
		tds->automaton_id = op_info->automaton_id;
		tds->destination_set = new set<State>();
		tds->destination_set->insert(op_info->trapstate);

		op_info->had_effect = true;
		return make_set_leaf(tds);
	} else if (mtbdd_isleaf(dd)) {
		// @TODO: Check that the complete with trapstate does have the same automaton ID as the node 
		// being returned.
		return dd;
	}
	return mtbdd_invalid;
}

/**
 * Performs pad closure on the two given transitions. Information about which states are final
 * and (out) information about whether the pad closure had any effect (the left transition was 
 * modified) are passed via `op_param` (a pointer to Pad_Closure_Info struct).
 */
TASK_IMPL_3(MTBDD, pad_closure_op, MTBDD *, p_left, MTBDD *, p_right, uint64_t, op_param) {
	MTBDD left = *p_left, right = *p_right;
    if (left == mtbdd_false) {
        return right; // Padding closure has no effect, when one mtbdd path leads to nothing
    }
    if (right == mtbdd_false) {
        return left; // Same here
    }
    // If both are non-false leaves, we can try doing the actual padding closure
    if (mtbdd_isleaf(left) && mtbdd_isleaf(right)) {
        auto left_tds  = (Transition_Destination_Set*) mtbdd_getvalue(left);
        auto right_tds = (Transition_Destination_Set*) mtbdd_getvalue(right);

	    auto pci = (Pad_Closure_Info*) op_param;

		// Check whether the transition destination state even leads the the right state (might not)
		if (left_tds->destination_set->find(pci->right_state) == left_tds->destination_set->end()) {
			// Does not lead to the right state, therefore cannot propagate padding.
			return left; 
		}

	    bool is_final;
	    for (State rs: *right_tds->destination_set)
	    {
			is_final = false;

			// Is the current right state final?
			for (uint32_t i = 0; i < pci->final_states_cnt; i++) {
				if (rs == pci->final_states[i]) {
					is_final = true;
					break; // We have the information we required, we don't need to iterate further.
				}
			}
			
			if (is_final) {
				const bool is_missing_from_left = left_tds->destination_set->find(rs) == left_tds->destination_set->end();
				if (is_missing_from_left) {
					// We've located a state that is final, and is not in left, pad closure adds it to the left states.
					left_tds->destination_set->insert(rs); // The actual pad closure
					pci->had_effect = true; // To utilize the python cache
				}
			}
	    }

		// The pad closure is in place modification 
		// (no new MTBDD is created, only leaf states are modified)
        return left;
    }

    return mtbdd_invalid;
}


TASK_IMPL_2(MTBDD, rename_states_op, MTBDD, dd, uint64_t, param) {
	if (dd == mtbdd_false) return mtbdd_false;

	if (mtbdd_isleaf(dd)) {
		(void) param;
		auto state_rename_info = STATE_RENAME_OP_PARAM;
		auto old_tds = (Transition_Destination_Set *) mtbdd_getvalue(dd);
		auto new_tds = new Transition_Destination_Set();

		new_tds->automaton_id = old_tds->automaton_id;
		auto renamed_leaf_contents = new set<State>();
		
		for (auto state : *old_tds->destination_set) {
			auto new_name_it = state_rename_info->states_rename_map->find(state);
			assert(new_name_it != state_rename_info->states_rename_map->end());

			auto new_state_name = new_name_it->second;
			renamed_leaf_contents->insert(new_state_name);
		}

		new_tds->destination_set = renamed_leaf_contents;

		return mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) new_tds);
	}

	return mtbdd_invalid;
}
