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
#include <stdio.h>

#include <algorithm>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::string;
using std::map;
using std::pair;

using sylvan::MTBDD;
using sylvan::mtbdd_makeleaf;
using sylvan::mtbdd_false;
using sylvan::mtbdd_isleaf;
using sylvan::mtbdd_getvalue;
using sylvan::mtbdd_invalid;
using sylvan::mtbdd_applyp_CALL;
using sylvan::mtbdd_enum_first;
using sylvan::mtbdd_enum_next;

extern uint64_t mtbdd_leaf_type_set;

/**
 * Global variables:
 */
void*       REMOVE_STATES_OP_PARAM = NULL;
uint64_t    REMOVE_STATES_OP_COUNTER = 0;
void*       ADD_TRAPSTATE_OP_PARAM = NULL;

State_Rename_Op_Info *STATE_RENAME_OP_PARAM = NULL;
Transform_Macrostates_To_Ints_State *TRANSFORM_MACROSTATES_TO_INTS_STATE = NULL;

Pad_Closure_Info *PAD_CLOSURE_OP_STATE = NULL;
uint64_t    PAD_CLOSURE_OP_COUNTER = 64;

uint64_t  ADD_TRAPSTATE_OP_COUNTER = (1LL << 32);
uint64_t  STATE_RENAME_OP_COUNTER = (1LL << 33);
uint64_t  TRANSFORM_MACROSTATES_TO_INTS_COUNTER = (1LL << 34);

/**
 * Unites the two transition MTBDDs.
 */
TASK_IMPL_3(MTBDD, transitions_union_op, MTBDD *, left_op_ptr, MTBDD *, right_op_ptr, uint64_t, param)
{
    MTBDD left_mtbdd = *left_op_ptr, right_mtbdd = *right_op_ptr;

    if (left_mtbdd == mtbdd_false)  return right_mtbdd;
    if (right_mtbdd == mtbdd_false) return left_mtbdd;

    if (mtbdd_isleaf(left_mtbdd) && mtbdd_isleaf(right_mtbdd)) {

        auto left_contents  = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(left_mtbdd));
        auto right_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(right_mtbdd));

        Transition_Destination_Set leaf_contents;

        std::set_union(left_contents->destination_set.begin(), left_contents->destination_set.end(),
                       right_contents->destination_set.begin(), right_contents->destination_set.end(),
                       std::inserter(leaf_contents.destination_set, leaf_contents.destination_set.begin()));

        MTBDD union_leaf = make_set_leaf(&leaf_contents);
        return union_leaf;
    }

    return sylvan::mtbdd_invalid;
}

/**
 * 'Abstraction' operation defines what happens to MTBDD subtrees when a variable is projected away.
 * We perform an union of the subtrees trees.
 */
TASK_IMPL_3(MTBDD, project_variable_away_abstract_op, MTBDD, left_mtbdd, MTBDD, right_mtbdd, int, k)
{
    MTBDD u = mtbdd_applyp(left_mtbdd, right_mtbdd, (uint64_t)-1, TASK(transitions_union_op), AMAYA_EXISTS_OPERATION_ID);
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

        for (auto state : tds->destination_set) {
            bool should_be_removed = states_to_remove->find(state) != states_to_remove->end();
            if (should_be_removed) new_tds->destination_set.erase(state);
        }

        if (new_tds->destination_set.empty()) {
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
 * Completes the given transition MTBDD by adding transitions along missing symbols leading to a trapstate.
 *
 * NOTE:
 * The param is not used - sylvan cachce problems, see remove_states_op. Instead the information
 * needed for the operation (like automaton ID and trapstate value) is passed via global variable
 * ADD_TRAPSTATE_OP_PARAM. This should be used with combination with ADD_TRAPSTATE_OP_COUNTER to
 * achieve cache utilization only for the the current completition.
 */
TASK_IMPL_2(MTBDD, complete_transition_with_trapstate_op, MTBDD, dd, uint64_t, param)
{
    // param is used just to avoid sylvan cache-miss
    if (dd == mtbdd_false) {
        (void) param;
        auto op_info = reinterpret_cast<Complete_With_Trapstate_Op_Info*>(ADD_TRAPSTATE_OP_PARAM);

        Transition_Destination_Set leaf_contents;
        leaf_contents.destination_set.insert(op_info->trapstate);

        op_info->had_effect = true;
        MTBDD leaf = make_set_leaf(&leaf_contents);
        return leaf;
    } else if (mtbdd_isleaf(dd)) {
        return dd;
    }
    return mtbdd_invalid;
}


inline bool contains_final_state(Pad_Closure_Info* pci, Transition_Destination_Set* tds)
{
    for (State current_post_state: tds->destination_set) {
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
        auto left_tds  = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(left));
        auto right_tds = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(right));

        auto pci = PAD_CLOSURE_OP_STATE;

        // Check whether the MTBDD of the pre-state (left) even leads to the current state (right)
        if (left_tds->destination_set.find(pci->right_state) == left_tds->destination_set.end()) {
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

        // The saturation property is broken, fix it by adding the new final state to the pre-state (left) TDS
        Transition_Destination_Set new_leaf_contents(*left_tds);
        new_leaf_contents.destination_set.insert(pci->new_final_state);
        MTBDD leaf = make_set_leaf(&new_leaf_contents);

        return leaf;
    }

    return mtbdd_invalid;
}


TASK_IMPL_2(MTBDD, rename_states_op, MTBDD, dd, uint64_t, param) {
    if (dd == mtbdd_false) return mtbdd_false;

    if (mtbdd_isleaf(dd)) {
        (void) param;
        auto state_rename_info = STATE_RENAME_OP_PARAM;
        auto old_tds = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(dd));

        // Do a heap allocation here, so that we can call mtbdd_makeleaf directly, avoiding copying state contents.
        auto new_leaf_contents = new Transition_Destination_Set();

        for (auto state : old_tds->destination_set) {
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
            new_leaf_contents->destination_set.insert(new_state_name);
        }

        return mtbdd_makeleaf(mtbdd_leaf_type_set, reinterpret_cast<uint64_t>(new_leaf_contents));
    }

    return mtbdd_invalid;
}


TASK_IMPL_2(MTBDD, transform_macrostates_to_ints_op, MTBDD, dd, uint64_t, param) {
    if (dd == mtbdd_false) return mtbdd_false;

    if (mtbdd_isleaf(dd)) {
        (void) param;

        auto transform_state = TRANSFORM_MACROSTATES_TO_INTS_STATE;
        auto old_tds = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(dd));

        State macrostate_state_number;
        bool is_cache_miss = false;
        auto iterator = transform_state->alias_map->find(old_tds->destination_set);
        if (iterator == transform_state->alias_map->end()) {
            macrostate_state_number = transform_state->first_available_state_number++;
        } else {
            // Cache entry for this leaf must have gotten evicted, we need to
            // return the previously returned leaf with the same alias number.
            is_cache_miss = true;
            macrostate_state_number = iterator->second;
        }

        // @Warn: This relies on the fact that the state sets are represented in a canoical fashion - the std::set
        //        keeps them sorted. That means that two macrostates e.g {1, 2, 3} and {3, 2, 1} will get always hashed to the
        //        same value --- Otherwise the same macrostates would get more than 1 ID which would cause troubles.

        Transition_Destination_Set new_leaf_contents;
        new_leaf_contents.destination_set.insert(macrostate_state_number);

        if (!is_cache_miss) {
            // Serialize the current macrostate, so that the python side will get notified about the created mapping.
            for (auto state : old_tds->destination_set) {
                transform_state->serialized_macrostates->push_back(state);
            }

            transform_state->macrostates_sizes->push_back(old_tds->destination_set.size());
            transform_state->macrostates_cnt += 1;

            transform_state->alias_map->emplace(old_tds->destination_set, macrostate_state_number);
        }

        MTBDD leaf = make_set_leaf(&new_leaf_contents);
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


void write_mtbdd_dot_to_tmp_file(const std::string& filename, sylvan::MTBDD mtbdd)
{
    std::string output_path = "/tmp/" + filename;
    auto output_file_handle = fopen(output_path.c_str(), "w");
    sylvan::mtbdd_fprintdot(output_file_handle, mtbdd);
    fclose(output_file_handle);
}

template<typename T>
std::string array_to_str(T* array, const uint64_t array_size)
{
    stringstream string_stream;
    string_stream << "[";
    if (!array_size) {
        string_stream << "]";
        return string_stream.str();
    }

    string_stream << (uint64_t) array[0];
    for (uint64_t i = 1; i < array_size; i++) {
        string_stream << "," << (uint64_t) array[i];
    }
    string_stream << "]";
    return string_stream.str();
}


std::pair<std::set<State>, std::set<State>>
fragment_dest_states_using_partition(const std::set<State>& dest_states, const std::set<State>& partition)
{
    std::set<State> intersection, difference;

    std::set_intersection(dest_states.begin(), dest_states.end(),
                          partition.begin(), partition.end(),
                          std::inserter(intersection, intersection.begin()));

    // @Optimize: Do not compute the difference if we know that intersection is the same as dest_states
    std::set_difference(dest_states.begin(), dest_states.end(),
                        intersection.begin(), intersection.end(),
                        std::inserter(difference, difference.begin()));

    return {intersection, difference};
}


inline bool should_explore_partition(std::vector<std::set<State>>& to_explore,
                                     std::set<std::set<State>>& known_partitions,
                                     std::set<State>& new_partition)
{
    bool is_partition_in_to_explore = std::find(to_explore.begin(), to_explore.end(), new_partition) != to_explore.end();
    // bool is_partition_known = known_partitions.find(new_partition) != known_partitions.end();
    return !is_partition_in_to_explore; // !is_partition_known;
}

inline sylvan::MTBDD traverse_mtbdd_along_variables(sylvan::MTBDD mtbdd, uint32_t* vars, uint64_t var_count, uint8_t* path)
{
    // We have to do a more careful traversal, as we generate path in the overall MTBDD for the partition
    // (union of all MTBDDs for the states in the partition), therefore, the MTBDD `mtbdd` for state
    // or the overall MTBDD might have different variables as don't care.
    auto top_mtbdd_var = sylvan::mtbdd_getvar(mtbdd);
    uint64_t var_index = 0;
    while (var_index < var_count && !sylvan::mtbdd_isleaf(mtbdd)) {
        if (top_mtbdd_var == vars[var_index]) {
            switch (path[var_index]) {
                case 0:
                    mtbdd = sylvan::mtbdd_getlow(mtbdd);
                    break;
                case 1:
                    mtbdd = sylvan::mtbdd_gethigh(mtbdd);
                    break;
                default:
                    // Overall MTBDD is don't-care for this variable, but the state MTBDD cares. As the overall
                    // MTBDD is a union of all state MTBDDs, if it is a don't care for some variable, then it means
                    // that the same subtrees have been the result of the union operation. Therefore, does not matter
                    // which way we go in the state MTBDD, as long as we have a consistent decision strategy in all
                    // state MTBDDs. Using a consistent strategy means that in an individual state MTBDD we might reach
                    // a different leaf as when using some other strategy, however, as the subtrees in the overall MTBDD
                    // are the same, there must exist other state(s) from the partition, that will contribute with the missing
                    // destination states.
                    mtbdd = sylvan::mtbdd_getlow(mtbdd);
                    break;
            }
            top_mtbdd_var = sylvan::mtbdd_getvar(mtbdd);
        }
        var_index++;
    }

    assert(mtbdd_isleaf(mtbdd));
    return mtbdd;
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
    std::cout << "Minimizing automaton " << nfa << std::endl;
#endif

    std::vector<std::set<State>> partitions_to_check {nfa.final_states, nonfinal_states_ordered};
    std::set<std::set<State>> existing_partitions {nfa.final_states, nonfinal_states_ordered};

    uint8_t path_in_mtbdd_to_leaf[nfa.var_count];
    uint32_t vars_to_check[nfa.var_count];
    sylvan::mtbdd_set_to_array(nfa.vars, vars_to_check);

    // Refine partitions untill fixpoint
    while (!partitions_to_check.empty()) {
        auto current_partition = partitions_to_check.back();
        bool was_current_partition_fragmented = false;
        partitions_to_check.pop_back();

#if DEBUG
        std::cout << "Processing partition: " << states_to_str(current_partition) << " (size=" << current_partition.size() << ")" << std::endl;
#endif

        // Compute MTBDD encoding all outgoing transitions of the current partition
        MTBDD partition_mtbdd = mtbdd_false;
        sylvan::mtbdd_refs_push(partition_mtbdd);
        for (auto state : current_partition){
            MTBDD mtbdd_for_current_partition_component = nfa.transitions[state];

            MTBDD new_partition_mtbdd = mtbdd_applyp(partition_mtbdd,
                                                     mtbdd_for_current_partition_component,
                                                     (uint64_t) 0,
                                                     TASK(transitions_union_op),
                                                     AMAYA_UNION_OP_ID);
            sylvan::mtbdd_refs_pop(1);
            sylvan::mtbdd_refs_push(new_partition_mtbdd);
            partition_mtbdd = new_partition_mtbdd;
        }

        // Iterate over all leaves of the created MTBDD. Every such a leaf represents a single post set over a symbol
        // given implicitly by the path in the MTBDD.
        MTBDD leaf = mtbdd_enum_first(partition_mtbdd, nfa.vars, path_in_mtbdd_to_leaf, NULL);

        while (leaf != mtbdd_false && !was_current_partition_fragmented) {
            // See whether the states reachable via this symbol belong to the same equivalence class
            Transition_Destination_Set* leaf_contents = (Transition_Destination_Set*) mtbdd_getvalue(leaf);

            // If the entire partition leads to one state it cannot be fragmented using this single state
            if (leaf_contents->destination_set.size() == 1) {
                leaf = mtbdd_enum_next(partition_mtbdd, nfa.vars, path_in_mtbdd_to_leaf, NULL);
                continue;
            }

            // Check if the entire destination set belongs to same equivalence class, or we need to fragment it
            for (auto existing_partition: existing_partitions) {
                // @Optimize: Pass in references to sets - reuse them, and avoid needless allocations
                auto fragment = fragment_dest_states_using_partition(leaf_contents->destination_set, existing_partition);
                auto dest_states_from_existing_partition = fragment.first;
                auto dest_states_not_from_existing_partition = fragment.second;

                if (dest_states_from_existing_partition.size() && dest_states_not_from_existing_partition.size()) {
                    // @Optimize: For now, we compute fragments of the current partition iteratively. Instead, we should keep the information
                    //            about what state led to a component of the destination set in the MTBDD

                    std::set<State> fragment_leading_to_partition, fragment_not_leading_to_partition;

                    // Iterate over the states in current partition and divide them into those that lead to the existing partition, and those that do not
                    for (auto state: current_partition) {
                        MTBDD state_mtbdd = nfa.transitions[state];

                        state_mtbdd = traverse_mtbdd_along_variables(state_mtbdd, vars_to_check, nfa.var_count, path_in_mtbdd_to_leaf);

                        Transition_Destination_Set* tds = (Transition_Destination_Set*) sylvan::mtbdd_getvalue(state_mtbdd);


                        assert(tds->destination_set.size() == 1);  // We should be dealing with complete DFAs
                        State dest_state = *(tds->destination_set.begin());

                        if (existing_partition.find(dest_state) != existing_partition.end()) {
                            fragment_leading_to_partition.insert(state);
                        } else {
                            fragment_not_leading_to_partition.insert(state);
                        }
                    }


#if DEBUG
                    assert(fragment_not_leading_to_partition.size());
                    assert(fragment_leading_to_partition.size());
#endif

                    // Update the overall partitions with the fragments
                    existing_partitions.erase(current_partition);
                    existing_partitions.insert(fragment_leading_to_partition);
                    existing_partitions.insert(fragment_not_leading_to_partition);

                    if (should_explore_partition(partitions_to_check, existing_partitions, fragment_leading_to_partition)) {
                        partitions_to_check.push_back(fragment_leading_to_partition);
                    }

                    if (should_explore_partition(partitions_to_check, existing_partitions, fragment_not_leading_to_partition)) {
                        partitions_to_check.push_back(fragment_not_leading_to_partition);
                    }

                    was_current_partition_fragmented = true;
                    break;
                }
            }

            leaf = mtbdd_enum_next(partition_mtbdd, nfa.vars, path_in_mtbdd_to_leaf, NULL);
        }
        sylvan::mtbdd_refs_pop(1);
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
        sylvan::mtbdd_refs_push(mtbdd_for_current_partition_index);

        MTBDD mtbdd_for_some_state              = nfa.transitions[some_state];
        MTBDD leaf                              = mtbdd_enum_first(mtbdd_for_some_state, nfa.vars, path_in_mtbdd_to_leaf, NULL);

        while (leaf != mtbdd_false) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
            State tds_state = *tds->destination_set.begin();
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
            destination_tds.destination_set = {dest_partition_index};

            MTBDD dest_leaf  = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) &destination_tds);
            MTBDD dest_mtbdd = sylvan::mtbdd_cube(nfa.vars, path_in_mtbdd_to_leaf, dest_leaf);

            // Add the created MTBDD to the other transitions originating in the current state
            MTBDD new_mtbdd_for_current_part_idx = mtbdd_applyp(mtbdd_for_current_partition_index, dest_mtbdd, (uint64_t) 0, TASK(transitions_union_op), AMAYA_UNION_OP_ID);
            sylvan::mtbdd_refs_pop(1);
            sylvan::mtbdd_refs_push(new_mtbdd_for_current_part_idx);
            mtbdd_for_current_partition_index = new_mtbdd_for_current_part_idx;

            // Check whether we have already explored the destination partition
            if (!was_dest_partition_explored) {
                partitions_to_add_to_minimized_dfa.insert(tds_partition_ptr);
            } else {
#if DEBUG
                std::cout << "Reachable partition " << states_to_str(*tds_partition_ptr) << " aka " << tds_state << " was already explored." << std::endl;
#endif
            }

            leaf = mtbdd_enum_next(mtbdd_for_some_state, nfa.vars, path_in_mtbdd_to_leaf, NULL);
        }

        // @Note: We cannot ref all transition mtbdds at the end becaouse GC might trigger during the execution
        //        of the contruction of the result
        sylvan::mtbdd_ref(mtbdd_for_current_partition_index);
        sylvan::mtbdd_refs_pop(1);
        result_nfa.transitions[partition_index] = mtbdd_for_current_partition_index;
    }
#if DEBUG
    std::cout << "Result has #states=" << result_nfa.states.size() << std::endl;
#endif

    return result_nfa;
}

TASK_IMPL_3(MTBDD, build_pad_closure_fronier_op, MTBDD *, p_extension, MTBDD *, p_frontier, u64, raw_extension_origin_state)
{
    State extension_origin_state = static_cast<State>(raw_extension_origin_state);
    MTBDD extension = *p_extension, frontier = *p_frontier;

    if (extension == mtbdd_false) {
        return frontier;
    }

    if (frontier == mtbdd_false) {
        return mtbdd_false;
    }

    if (mtbdd_isleaf(extension) && mtbdd_isleaf(frontier)) {
        auto extension_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(extension));
        auto frontier_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(frontier));


        if (frontier_contents->destination_set.contains(extension_origin_state)) {
            return frontier;
        }

        if (!is_set_intersection_empty(extension_contents->destination_set, frontier_contents->destination_set)) {
            // It is possible to reach a state from extension_origin_state by a single transition along the current
            // symbol from which it is possible to reach a final state along the same symbol by reading finitely many
            // such symbols
            Transition_Destination_Set new_frontier_contents(*frontier_contents);
            new_frontier_contents.destination_set.insert(extension_origin_state);
            MTBDD new_frontier = make_set_leaf(&new_frontier_contents);
            return new_frontier;
        }

        return frontier;
    }

    return mtbdd_invalid;
}

TASK_IMPL_3(MTBDD, add_pad_transitions_op, MTBDD *, p_transitions, MTBDD *, p_frontier, u64, raw_pad_closure_info)
{
    Pad_Closure_Info2* pad_closure_info = reinterpret_cast<Pad_Closure_Info2*>(raw_pad_closure_info);
    MTBDD transitions = *p_transitions, frontier = *p_frontier;

    if (transitions == mtbdd_false) {
        return mtbdd_false;
    }

    if (frontier == mtbdd_false) {
        return transitions;
    }

    if (mtbdd_isleaf(transitions) && mtbdd_isleaf(frontier)) {
        auto transition_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(transitions));
        auto frontier_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(frontier));

        // It must be possible to reach a final state from the current origin
        if (!frontier_contents->destination_set.contains(pad_closure_info->origin_state))
            return transitions;

        // Make sure that also its post is contained in the frontier as states that are final are contained in the frontier
        // but they do not necesserily can have a path to another final state via the symbol
        // @FixMe: Maybe we should fix the frontier instead to not contain the final states (exclude epsilon reachablility)
        if (is_set_intersection_empty(transition_contents->destination_set, frontier_contents->destination_set)) {
            return transitions;
        }

        if (is_set_intersection_empty(transition_contents->destination_set, *pad_closure_info->final_states)) {
            Transition_Destination_Set new_transition_contents(*transition_contents);
            new_transition_contents.destination_set.insert(pad_closure_info->new_final_state);
            MTBDD new_destinations = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, reinterpret_cast<u64>(&new_transition_contents));
            return new_destinations;
        }

        return transitions;
    }

    return mtbdd_invalid;
}

inline
bool skip_product_state_due_to_no_post(Intersection_Info2& info, std::pair<State, State>& product) {
#ifdef INTERSECTION_DETECT_STATES_WITH_NO_POST
    auto left_pos_in_postless = std::find(info.left_states_without_post.begin(), info.left_states_without_post.end(), product.first);
    auto right_pos_in_postless = std::find(info.right_states_without_post.begin(), info.right_states_without_post.end(), product.second);

    bool is_left_postless = (left_pos_in_postless != info.left_states_without_post.end());
    bool is_right_postless = (right_pos_in_postless != info.right_states_without_post.end());

    if (is_left_postless || is_right_postless) {
        bool product_has_lang = (info.left_final_states->contains(product.first) && info.right_final_states->contains(product.second));
        return !product_has_lang;
    }
    return false;
#else
    return false;
#endif
}



TASK_IMPL_3(MTBDD, transitions_intersection2_op, MTBDD *, pa, MTBDD *, pb, uint64_t, param) {
    MTBDD a = *pa, b = *pb;

    if (a == mtbdd_false || b == mtbdd_false) {
        return mtbdd_false; // The result is an empty set when one of the operands is an empty set (mtbdd_false)
    }

    if (!mtbdd_isleaf(a) || !mtbdd_isleaf(b)) {
        // @Note: According to the GMP implementation in Sylvan source code, swapping pointers should increase the cache
        // performance. However, doing a pointer swap causes some of the tests to fail, not sure what is up with it.
        return mtbdd_invalid;
    }

    auto intersect_info = reinterpret_cast<Intersection_Info2*>(param);

    auto left_leaf_contents  = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(a));
    auto right_leaf_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(b));

    if (left_leaf_contents->destination_set.empty() || right_leaf_contents->destination_set.empty()) {
        return mtbdd_false;  // Equivalent to a leaf with an empty set
    }

    Transition_Destination_Set leaf_contents;

    for (auto left_state : left_leaf_contents->destination_set) {
        for (auto right_state : right_leaf_contents->destination_set) {
            std::pair<State, State> product_state = std::make_pair(left_state, right_state);

            if (skip_product_state_due_to_no_post(*intersect_info, product_state)) {
                ++intersect_info->skipped_states_with_no_post;
                continue;
            }

            State product_handle = intersect_info->seen_products.size();

            auto [existing_entry_it, did_insert] = intersect_info->seen_products.emplace(product_state, product_handle);

            Intersection_Discovery discovery = {.left = left_state, .right = right_state, .handle = product_handle};
            if (did_insert) {
                intersect_info->work_queue.push_back(discovery);
            } else {
                product_handle = existing_entry_it->second;
            }
            leaf_contents.destination_set.insert(product_handle);
        }
    }

    MTBDD intersection_leaf = make_set_leaf(&leaf_contents);
    return intersection_leaf;
}

TASK_IMPL_2(MTBDD, replace_macrostates_with_handles_op, MTBDD, dd, uint64_t, param) {
    if (dd == mtbdd_false) {
        auto ctx = reinterpret_cast<Determinization_Context*>(param);
        ctx->is_trapstate_needed = true;

        Transition_Destination_Set new_leaf_contents;
        new_leaf_contents.destination_set.insert(ctx->trapstate_handle);
        MTBDD new_leaf = make_set_leaf(&new_leaf_contents);

        return new_leaf;
    }

    if (!sylvan::mtbdd_isleaf(dd)) return sylvan::mtbdd_invalid;

    auto raw_leaf_contents = sylvan::mtbdd_getvalue(dd);
    auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(raw_leaf_contents);
    Macrostate macrostate(leaf_contents->destination_set.begin(), leaf_contents->destination_set.end());

    auto ctx = reinterpret_cast<Determinization_Context*>(param);
    auto [known_macrostates_elem, was_inserted] = ctx->known_macrostates.emplace(macrostate, ctx->known_macrostates.size());

    if (was_inserted) ctx->work_queue.push_back(&(*known_macrostates_elem));

    State macrostate_handle = known_macrostates_elem->second;
    Transition_Destination_Set new_leaf_contents;
    new_leaf_contents.destination_set.insert(macrostate_handle);

    MTBDD new_leaf = make_set_leaf(&new_leaf_contents);
    return new_leaf;
}

TASK_IMPL_2(MTBDD, remove_states2_op, MTBDD, dd, uint64_t, param) {
    if (dd == mtbdd_false) return mtbdd_false;

    if (mtbdd_isleaf(dd)) {
        auto states_to_remove = reinterpret_cast<std::set<State>*>(param);
        auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(mtbdd_getvalue(dd));

        Transition_Destination_Set new_leaf_contents;
        std::set_difference(leaf_contents->destination_set.begin(), leaf_contents->destination_set.end(),
                            states_to_remove->begin(), states_to_remove->end(),
                            std::inserter(new_leaf_contents.destination_set, new_leaf_contents.destination_set.begin()));

        MTBDD leaf = make_set_leaf(&new_leaf_contents);
        return leaf;
    }

    return mtbdd_invalid;
}
