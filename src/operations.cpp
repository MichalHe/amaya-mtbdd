#include "../include/operations.hpp"
#include "../include/custom_leaf.hpp"
#include "../include/base.hpp"

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
#include <string>
#include <sstream>

using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::string;
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
using sylvan::mtbdd_support_CALL;
using sylvan::mtbdd_enum_first;
using sylvan::mtbdd_enum_next;
using sylvan::mtbdd_set_count;

extern uint64_t mtbdd_leaf_type_set;

/**
 * Global variables:
 */
void*		REMOVE_STATES_OP_PARAM = NULL;
uint64_t 	REMOVE_STATES_OP_COUNTER = 0;
void*		ADD_TRAPSTATE_OP_PARAM = NULL;

Intersection_State* intersection_state = NULL;
uint64_t INTERSECTION_OP_COUNTER = (1LL << 35);

State_Rename_Op_Info *STATE_RENAME_OP_PARAM = NULL;
Transform_Metastates_To_Ints_State *TRANSFORM_METASTATES_TO_INTS_STATE = NULL;

Pad_Closure_Info *PAD_CLOSURE_OP_STATE = NULL;
uint64_t 	PAD_CLOSURE_OP_COUNTER = 64;

uint64_t 	ADD_TRAPSTATE_OP_COUNTER = (1LL << 32);
uint64_t 	STATE_RENAME_OP_COUNTER = (1LL << 33);
uint64_t    TRANSFORM_METASTATES_TO_INTS_COUNTER = (1LL << 34);

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

        auto intersection_tds = new Transition_Destination_Set(intersection_leaf_states);
        MTBDD intersection_leaf = make_set_leaf(intersection_tds);
		delete intersection_tds;
        return intersection_leaf;
    }

    // TODO: Perform pointer swap here, so the cache would be utilized (operation is commutative).
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

		// @RefactorMe
        //if (param == -1)
        //{
            //param = (uint32_t)tds_a.automaton_id;
        //}

        std::set<State> *union_set = new std::set<State>();
        std::set_union(
            tds_a.destination_set->begin(), tds_a.destination_set->end(),
            tds_b.destination_set->begin(), tds_b.destination_set->end(),
            std::inserter(*union_set, union_set->begin()));

        Transition_Destination_Set *union_tds = new Transition_Destination_Set(union_set);

        MTBDD union_leaf = make_set_leaf(union_tds); // Wrap the TDS with a MTBDD leaf.
		delete union_tds;
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

		MTBDD leaf = make_set_leaf(new_tds);
		delete new_tds;
		return leaf;
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
		// RefactorMe
		//tds->automaton_id = op_info->automaton_id;
		tds->destination_set = new set<State>();
		tds->destination_set->insert(op_info->trapstate);

		op_info->had_effect = true;
		MTBDD leaf =  make_set_leaf(tds);
		delete tds;
		return leaf;
	} else if (mtbdd_isleaf(dd)) {
		// @TODO: Check that the complete with trapstate does have the same automaton ID as the node
		// being returned.
		return dd;
	}
	return mtbdd_invalid;
}


inline bool contains_final_state(Pad_Closure_Info* pci, Transition_Destination_Set* tds)
{
    for (State current_post_state: *tds->destination_set) {
        for (uint32_t i = 0; i < pci->final_states_cnt; i++) {
            if (current_post_state == pci->final_states[i])
                return true;
        }
    }
    return false;
}


/**
 * Performs pad closure on the two given transitions. Information about which states are final
 * and (out) information about whether the pad closure had any effect (the left transition was
 * modified) are passed via `op_param` (a pointer to Pad_Closure_Info struct).
 */
TASK_IMPL_3(MTBDD, pad_closure_op, MTBDD *, p_left, MTBDD *, p_right, uint64_t, op_param)
{
	MTBDD left = *p_left, right = *p_right;

    if (left == mtbdd_false) {
        // The pre-state (left) leads to nothing, nothing can be propagated
        return mtbdd_false;
    }

    if (right == mtbdd_false) {
        // The state (right) leads to nothing, nothing can be propagated
        // Return the unmodified left leaf (no change)
        return left;
    }

    if (mtbdd_isleaf(left) && mtbdd_isleaf(right)) {
        // Both are leaves, do the actual padding closure
        auto left_tds  = (Transition_Destination_Set*) mtbdd_getvalue(left);
        auto right_tds = (Transition_Destination_Set*) mtbdd_getvalue(right);

	    auto pci = PAD_CLOSURE_OP_STATE;

		// Check whether the MTBDD of the pre-state (left) even leads to the current state (right)
		if (left_tds->destination_set->find(pci->right_state) == left_tds->destination_set->end()) {
			// Does not lead to the current state, therefore, the saturation property was not broken (here)
			return left;
		}

        // Check whether any final state is present in the post of the current state
        bool is_final_reachable_from_current = contains_final_state(pci, right_tds);

        if (!is_final_reachable_from_current) {
            // If there is no final state in the post of the current state
            // the saturation cannot be broken here
            return left; // No propagation
        }

        // Some final state is reachable from the current state via the current symbol
        // Check that some final state is reachable from the pre-state aswell, otherwise
        // the saturation property is broken
        bool is_final_reachable_from_pre = contains_final_state(pci, left_tds);

        if (is_final_reachable_from_pre) {
            // Some final state is reachable from both the pre-state, and the current state
            // the saturation property is satisfied.
            return left;  // Nothing to propagate
        }

        // The saturation property is broken, fix it by adding the new final
        // state to the pre-state (left) TDS
		Transition_Destination_Set* new_leaf_contents = new Transition_Destination_Set();

		// @Refactoring
		// new_leaf_contents->automaton_id = left_tds->automaton_id;

        // Make new state set from the original leaf + add the new final state
		auto new_leaf_destination_set = new std::set<State>();
		for (auto original_state: *left_tds->destination_set) new_leaf_destination_set->insert(original_state);
        new_leaf_destination_set->insert(pci->new_final_state);

		new_leaf_contents->destination_set = new_leaf_destination_set;

		MTBDD leaf = make_set_leaf(new_leaf_contents);
		delete new_leaf_contents;

		return leaf;
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

		// @Refactoring(codeboy): automaton_id is not used anymore, that is why it is commented out
		// new_tds->automaton_id = old_tds->automaton_id;
		auto renamed_leaf_contents = new set<State>();

		for (auto state : *old_tds->destination_set) {
			auto new_name_it = state_rename_info->states_rename_map->find(state);

			if (new_name_it == state_rename_info->states_rename_map->end()) {
				printf("We have found a leaf state with no mapping for it??");
				printf("State name: %lu\n", state);
				printf("Available mappings:");
				for (auto mapping : *state_rename_info->states_rename_map) {
					printf("(%lu, %lu)", mapping.first, mapping.second);
				}

				assert(false);
			}

			auto new_state_name = new_name_it->second;
			renamed_leaf_contents->insert(new_state_name);
		}

		new_tds->destination_set = renamed_leaf_contents;

		return mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) new_tds);
	}

	return mtbdd_invalid;
}


TASK_IMPL_2(MTBDD, transform_metastates_to_ints_op, MTBDD, dd, uint64_t, param) {
	if (dd == mtbdd_false) return mtbdd_false;

	if (mtbdd_isleaf(dd)) {
		(void) param;

		auto transform_state = TRANSFORM_METASTATES_TO_INTS_STATE;
		auto old_tds = (Transition_Destination_Set *) mtbdd_getvalue(dd);
		auto new_tds = new Transition_Destination_Set();

		// @Refactoring: this is commented because we do not use automaton_ids anymore
		// new_tds->automaton_id = old_tds->automaton_id;

		State metastate_state_number;
		bool is_cache_miss = false;
		auto iterator = transform_state->alias_map->find(*old_tds->destination_set);
		if (iterator == transform_state->alias_map->end()) {
			metastate_state_number = transform_state->first_available_state_number++;
		} else {
			// Cache entry for this leaf must have gotten evicted,
			// we need to return the previously returned leaf with
			// the same alias number.
			is_cache_miss = true;
			metastate_state_number = iterator->second;
		}

		// @Warn: This relies on the fact that the state sets are represented in a canoical fashion - the std::set
		// 		  keeps them sorted. That means that two metastates e.g {1, 2, 3} and {3, 2, 1} will get always hashed to the
		// 		  same value --- Otherwise the same metastates would get more than 1 ID which would cause troubles.

		auto transformed_leaf_contents = new set<State>();
		transformed_leaf_contents->insert(metastate_state_number);
		new_tds->destination_set = transformed_leaf_contents;

		if (!is_cache_miss) {
			// Serialize the current metastate, so that the python side will get notified about the created mapping.
			for (auto state : *old_tds->destination_set) {
				transform_state->serialized_metastates->push_back(state);
			}

			transform_state->metastates_sizes->push_back(old_tds->destination_set->size());
			transform_state->metastates_cnt += 1;

			transform_state->alias_map->insert(std::make_pair(*old_tds->destination_set, metastate_state_number));
		}

		MTBDD leaf = make_set_leaf(new_tds);
		delete new_tds;
		return leaf;
	}

	return mtbdd_invalid;
}


template<typename C>
std::string states_to_str(C const& states) {
    std::stringstream partition_str;

    partition_str << "{";
    if (!states.empty()) {
        auto state_it = states.begin();
        partition_str << *(state_it++);
        for (; state_it != states.end(); ++state_it) {
            partition_str << ", " << *state_it;
        }
    }
    partition_str << "}";

    return partition_str.str();
}


std::string nfa_to_str(struct NFA& nfa) {
    std::stringstream nfa_str;
    nfa_str << "NFA:" << std::endl
            << "> States:         " << states_to_str(nfa.states) << std::endl
            << "> Final states:   " << states_to_str(nfa.final_states) << std::endl
            << "> Initial states: " << states_to_str(nfa.initial_states) << std::endl;
    return nfa_str.str();
}

std::pair<std::set<State>, std::set<State>>
fragment_dest_states_using_partition(const std::set<State>& dest_states, const std::set<State>& partition)
{
    std::set<State> intersection, difference;

    std::set_intersection(dest_states.begin(), dest_states.end(),
                          partition.begin(), partition.end(),
                          std::inserter(intersection, intersection.begin()));

    // @Optimize: Do not compute the difference if we know that intersection is the same as dest_states
    std::set_difference(partition.begin(), partition.end(),
                        intersection.begin(), intersection.end(),
                        std::inserter(difference, difference.begin()));

    return {intersection, difference};
}

struct NFA minimize_hopcroft(struct NFA& nfa)
{
    LACE_ME;

    // Create initial partitions
    std::set<State> nonfinal_states_ordered;
    std::set_difference(nfa.states.begin(), nfa.states.end(),
                        nfa.final_states.begin(), nfa.final_states.end(),
                        std::inserter(nonfinal_states_ordered, nonfinal_states_ordered.begin()));
#if DEBUG
    std::cout << "Minimizing automaton " << nfa_to_str(nfa) << std::endl;
#endif

    std::vector<std::set<State>> partitions_to_check {nfa.final_states, nonfinal_states_ordered};

    std::set<std::set<State>> existing_partitions {nfa.final_states, nonfinal_states_ordered};


    uint8_t path_in_mtbdd_to_leaf[nfa.var_count];

    // Refine partitions untill fixpoint
    while (!partitions_to_check.empty()) {
        auto current_partition = partitions_to_check.back();
        partitions_to_check.pop_back();

#if DEBUG
        std::cout << "Processing partition: " << states_to_str(current_partition) << std::endl;
#endif

        // Compute MTBDD encoding all outgoing transitions of the current partition
        MTBDD partition_mtbdd = mtbdd_false;
        for (auto state : current_partition){
            MTBDD mtbdd_for_current_partition_component = nfa.transitions[state];

            partition_mtbdd = mtbdd_applyp(partition_mtbdd,
                                           mtbdd_for_current_partition_component,
                                           (uint64_t) 0,
                                           TASK(transitions_union_op),
                                           AMAYA_UNION_OP_ID);
        }

        // Iterate over all leaves of the created MTBDD. Every such a leaf represents a single post set over a symbol
        // given implicitly by the path in the MTBDD.
        MTBDD    support      = mtbdd_support(partition_mtbdd);
        uint32_t support_size = mtbdd_set_count(support);

        MTBDD leaf = mtbdd_enum_first(partition_mtbdd, support, path_in_mtbdd_to_leaf, NULL);

        while (leaf != mtbdd_false) {
            // See whether the states reachable via this symbol belong to one equivalence class
            Transition_Destination_Set* leaf_contents = (Transition_Destination_Set*) mtbdd_getvalue(leaf);

            // If the entire partition leads to one state it cannot be fragmented using this single state
            if (leaf_contents->destination_set->size() == 1) {
                leaf = mtbdd_enum_next(partition_mtbdd, support, path_in_mtbdd_to_leaf, NULL);
                continue;
            }

            // Check if the entire destination set belongs to same equivalence class, or we need to fragment it
            for (auto existing_partition: existing_partitions) {
                // @Optimize: Pass in references to sets - reuse them, and avoid needless allocations
                auto fragment = fragment_dest_states_using_partition(*leaf_contents->destination_set, existing_partition);
                if (fragment.first.size() == leaf_contents->destination_set->size() || fragment.first.size() == 0) {
#if DEBUG
                    std::cout << "Unable to fragment " << states_to_str(*leaf_contents->destination_set)
                              << " using partition " << states_to_str(existing_partition)
                              << std::endl;
#endif

                } else {
#if DEBUG
                    std::cout << "Fragmenting " << states_to_str(*leaf_contents->destination_set)
                              << " with partition " << states_to_str(existing_partition)
                              << std::endl;
#endif
                    // @Optimize: For now, we compute fragments of the current partition iteratively. Instead, we should keep the information
                    //            about what state led to a component of the destination set in the MTBDD

                    std::set<State> fragment_leading_to_partition, fragment_not_leading_to_partition;

                    // Compute the fragments
                    for (auto state: current_partition) {
                        MTBDD state_mtbdd = nfa.transitions[state];
                        uint64_t path_index = 0;
                        while (!mtbdd_isleaf(state_mtbdd)) {
                            assert(path_in_mtbdd_to_leaf[path_index] != 2);
                            if (path_in_mtbdd_to_leaf[path_index]) state_mtbdd = sylvan::mtbdd_gethigh(state_mtbdd);
                            else state_mtbdd = sylvan::mtbdd_getlow(state_mtbdd);
                        }
                        Transition_Destination_Set* tds = (Transition_Destination_Set*) sylvan::mtbdd_getvalue(state_mtbdd);
                        State dest_state = *(tds->destination_set->begin());

                        if (existing_partition.find(dest_state) != existing_partition.end()) {
                            fragment_leading_to_partition.insert(state);
                        } else {
                            fragment_not_leading_to_partition.insert(state);
                        }
                    }
#if DEBUG
                    std::cout << "Fragmented " << states_to_str(current_partition) << " into " << states_to_str(fragment_leading_to_partition)
                              << " (states leading to current partition) and " << states_to_str(fragment_not_leading_to_partition)
                              << " (states not leading to current partition)" << std::endl;
#endif

                    // Update the overall partitions with the fragments
                    existing_partitions.erase(current_partition);
                    existing_partitions.insert(fragment_leading_to_partition);
                    existing_partitions.insert(fragment_not_leading_to_partition);

                    // We cannot check only smaller of the two fragments - imagine a situation where we produce a smaller fragment with only one state - that means we never fragment bigger the partition
                    auto smaller_fragment = fragment_leading_to_partition.size() < fragment_not_leading_to_partition.size() ? fragment_leading_to_partition : fragment_not_leading_to_partition;
                    if (std::find(partitions_to_check.begin(), partitions_to_check.end(), fragment_leading_to_partition) == partitions_to_check.end()) {
                        partitions_to_check.push_back(fragment_leading_to_partition);
                    }

                    if (std::find(partitions_to_check.begin(), partitions_to_check.end(), fragment_not_leading_to_partition) == partitions_to_check.end()) {
                        partitions_to_check.push_back(fragment_not_leading_to_partition);
                    }

                    break;
                }
            }

            leaf = mtbdd_enum_next(partition_mtbdd, support, path_in_mtbdd_to_leaf, NULL);
        }
    }

#if DEBUG
    std::cout << "Fixpoint found. Partitions:" << std::endl;
    for (auto partition: existing_partitions) {
        std::cout << " - " << states_to_str(partition) << std::endl;
    }
#endif

    // Construct the NFA with states resulting from condensation according to computed eq. partitions
    NFA result_nfa;
    result_nfa.vars      = nfa.vars;
    result_nfa.var_count = nfa.var_count;

    std::map<State, const std::set<State>*> state_to_its_parition;
    for (auto& partition: existing_partitions) {
        for (auto state: partition) {
            state_to_its_parition[state] = &partition;
        }
    }

    State partition_index = 0;
    State next_available_partition_index = 1;
    State initial_state = *nfa.initial_states.begin();  // It is an DFA, therefore there is only one final state
    auto initial_partition_ptr = state_to_its_parition[initial_state];

    // We cannot populate the partition_to_state_index right now, since some partitions
    // might not be reachable and we would like to avoid assigning them a state number
    // However, we know that the initial partition will have its state number = 0
    std::map<const std::set<State>*, State> partition_ptr_to_state_index {{initial_partition_ptr, 0}};
    result_nfa.initial_states.insert(0);

    // Start with the partition containing the initial state and discover reachable partitions from there
    std::set<const std::set<State>*> partitions_to_add_to_minimized_dfa{initial_partition_ptr};

    while (!partitions_to_add_to_minimized_dfa.empty()) {
        auto current_partition_ptr = *partitions_to_add_to_minimized_dfa.begin();
        partitions_to_add_to_minimized_dfa.erase(current_partition_ptr);

        auto partition_index_ptr = &partition_ptr_to_state_index[current_partition_ptr];
        partition_index = *partition_index_ptr; // Partition index will not have its uppermost bit 1
        *partition_index_ptr |= (1ul << 63);    // Set the topmost bit to 1 to indicate that we already did explore the partition

#if DEBUG
        printf("Processing state=%lu (value=%lu)\n", (partition_index & ~(1ul << 63)), partition_index);
#endif
        result_nfa.states.insert(partition_index);

        // It is sufficient to check only one state, as there cannot be a partition with containing a final and a nonfinal state
        State some_state = *current_partition_ptr->begin();
        if (nfa.final_states.find(some_state) != nfa.final_states.end())
            result_nfa.final_states.insert(partition_index);

        // Compute transitions - create MTBDD of this partition
        MTBDD mtbdd_for_current_partition_index = sylvan::mtbdd_false;
        MTBDD mtbdd_for_some_state              = nfa.transitions[some_state];
        MTBDD leaf                              = mtbdd_enum_first(mtbdd_for_some_state, nfa.vars, path_in_mtbdd_to_leaf, NULL);

        while (leaf != mtbdd_false) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
            State tds_state = *tds->destination_set->begin();
            auto tds_partition_ptr = state_to_its_parition[tds_state];

            // Check if we have already seen the destination partition before and it already has a number
            State dest_partition_index;
            bool was_dest_partition_explored = false;
            auto partition_ptr_to_state_index_iter = partition_ptr_to_state_index.find(tds_partition_ptr);
            if (partition_ptr_to_state_index_iter != partition_ptr_to_state_index.end()) {
                dest_partition_index = (*partition_ptr_to_state_index_iter).second;
                was_dest_partition_explored = (1ul << 63) & dest_partition_index; // Get the topmost bit - was explorated info
                dest_partition_index &= ~(1ul << 63); // Zero out the topmost bit
            } else {
                dest_partition_index = next_available_partition_index;
                partition_ptr_to_state_index[tds_partition_ptr] = next_available_partition_index++;
            }

            // Create a MTBDD encoding only the transition to the state resulting from condensing the partition in which was the state located
            Transition_Destination_Set destination_tds;
            destination_tds.destination_set = new std::set<State>{dest_partition_index};

            MTBDD dest_leaf  = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) &destination_tds);
            MTBDD dest_mtbdd = sylvan::mtbdd_cube(nfa.vars, path_in_mtbdd_to_leaf, dest_leaf);

            // Add the created MTBDD to the other transitions originating in the current state
            mtbdd_for_current_partition_index = mtbdd_applyp(mtbdd_for_current_partition_index, dest_mtbdd, (uint64_t) 0, TASK(transitions_union_op), AMAYA_UNION_OP_ID);

            // Check whether we have already explored the destination partition
            if (!was_dest_partition_explored) partitions_to_add_to_minimized_dfa.insert(tds_partition_ptr);

            leaf = mtbdd_enum_next(mtbdd_for_some_state, nfa.vars, path_in_mtbdd_to_leaf, NULL);
        }
        result_nfa.transitions[partition_index] = mtbdd_for_current_partition_index;
    }
    std::cout << "Result has #states=" << result_nfa.states.size() << std::endl;
    return result_nfa;
}
