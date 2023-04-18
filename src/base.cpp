#include "../include/base.hpp"
#include "../include/custom_leaf.hpp"
#include "../include/operations.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include <sstream>
#include <string>

using sylvan::MTBDD;


Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other) {
    // NOTE: Copy is called when a new leaf is created
    this->destination_set = std::set<State>(other.destination_set);
}

Transition_Destination_Set::Transition_Destination_Set(std::set<State>& destination_set) {
    // @Cleanup: Check whether is is even called as it does a copy of the destination set,
    //           and thus it is essentially the same as the copy constructor
    this->destination_set = destination_set;
}

void Transition_Destination_Set::print_dest_states() {
    int cnt = 1;
    std::cout << "{";
    for (auto state : this->destination_set)
    {
        std::cout << state;
        if (cnt < this->destination_set.size())
        {
            std::cout << ", ";
        }
        ++cnt;
    }
    std::cout << "}" << std::endl;
}


void unpack_dont_care_bits_in_mtbdd_path(
        uint8_t* mtbdd_path,
        const uint64_t path_size,
        std::vector<struct Transition>& unpacked_transitions,
        uint64_t symbol_index_to_start_at,
        const State origin,
        const State destination)
{
    // Find first don't-care symbol
    uint64_t first_dont_care_symbol_index = 0;
    bool dont_care_symbol_found = false;
    for (uint64_t i = symbol_index_to_start_at; i < path_size; ++i) {
        if (mtbdd_path[i] == 2) {
            dont_care_symbol_found = true;
            first_dont_care_symbol_index = i;
        }
    }

    if (dont_care_symbol_found) {
        uint8_t path_with_dont_care_bit_interpreted[path_size];
        for (uint64_t i = 0; i < path_size; i++) path_with_dont_care_bit_interpreted[i] = mtbdd_path[i];

        // Case split - interpret the don't care bit by setting first bit = 0, and then bit = 1
        path_with_dont_care_bit_interpreted[first_dont_care_symbol_index] = 0;
        unpack_dont_care_bits_in_mtbdd_path(path_with_dont_care_bit_interpreted, path_size, unpacked_transitions, first_dont_care_symbol_index+1, origin, destination);

        path_with_dont_care_bit_interpreted[first_dont_care_symbol_index] = 1;
        unpack_dont_care_bits_in_mtbdd_path(path_with_dont_care_bit_interpreted, path_size, unpacked_transitions, first_dont_care_symbol_index+1, origin, destination);
    } else {
        // The mtbdd_path does not contain any don't care bits, we can emit a transition
        std::vector<u8> symbol(path_size);
        for (uint64_t i = 0; i < path_size; i++) symbol[i] = mtbdd_path[i];

        Transition transition = {.from = origin, .to = destination, .symbol = symbol};
        unpacked_transitions.push_back(transition);
    }
}


std::vector<Transition> NFA::get_symbolic_transitions_for_state(State state) const {
    std::vector<Transition> symbolic_transitions;

    auto state_mtbdd_it = transitions.find(state);
    if (state_mtbdd_it == transitions.end()) {
        return {};
    }

    auto mtbdd = state_mtbdd_it->second;

    u8 raw_symbol[this->var_count];

    MTBDD leaf = sylvan::mtbdd_enum_first(mtbdd, this->vars, raw_symbol, NULL);
    while (leaf != sylvan::mtbdd_false) {
        auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(sylvan::mtbdd_getvalue(leaf));

        for (auto& dest_state: leaf_contents->destination_set) {

            std::vector<u8> symbol(this->var_count);
            for (u64 i = 0; i < this->var_count; i++) symbol[i] = raw_symbol[i];

            Transition transition = {.from = state, .to = dest_state, .symbol = symbol};
            symbolic_transitions.push_back(transition);
        }

        leaf = sylvan::mtbdd_enum_next(mtbdd, this->vars, raw_symbol, NULL);
    }

    return symbolic_transitions;
}

void NFA::add_transition(State from, State to, const u64 symbol, const u64 quantified_bits_mask) {
    u8 rich_symbol[var_count];

    for (u64 bit_i = 0u; bit_i < var_count; bit_i++) {
        u64 current_bit = (1u << bit_i);
        bool dont_care = (current_bit & quantified_bits_mask) > 0;
        u8 care_val = (symbol & current_bit) > 0;

        rich_symbol[bit_i] = dont_care ? 2 : care_val;
    }

    add_transition(from, to, rich_symbol);
}

void NFA::add_transition(State from, State to, std::vector<u8>&& symbol) {
    add_transition(from, to, symbol.data());
}

void NFA::add_transition(State from, State to, u8* rich_symbol) {
    Transition_Destination_Set leaf_contents({to});
    sylvan::MTBDD leaf = make_set_leaf(&leaf_contents);

    sylvan::MTBDD transition_mtbdd = sylvan::mtbdd_cube(this->vars, rich_symbol, leaf);

    LACE_ME;
    auto [present_key_val_it, was_inserted] = transitions.emplace(from, transition_mtbdd);
    if (was_inserted) {
        sylvan::mtbdd_ref(transition_mtbdd);
    } else {
        auto existing_transitions = present_key_val_it->second;

        using namespace sylvan; // Pull the entire sylvan namespace because mtbdd_applyp is a macro
        sylvan::MTBDD updated_mtbdd = mtbdd_applyp(existing_transitions, transition_mtbdd, 0u, TASK(transitions_union_op), AMAYA_UNION_OP_ID);

        sylvan::mtbdd_deref(existing_transitions);
        sylvan::mtbdd_ref(updated_mtbdd);

        transitions[from] = updated_mtbdd;
    }
}


std::vector<Transition> nfa_unpack_transitions(struct NFA& nfa) {
    LACE_ME;

    std::vector<struct Transition> transitions;

    uint8_t* path_in_mtbdd_to_leaf = (uint8_t*) malloc(sizeof(uint8_t) * nfa.var_count);
    assert(path_in_mtbdd_to_leaf != NULL);

    for (auto origin_state: nfa.states) {
        MTBDD state_transitions = nfa.transitions[origin_state];
        MTBDD leaf = sylvan::mtbdd_enum_first(state_transitions, nfa.vars, path_in_mtbdd_to_leaf, NULL);

        while (leaf != sylvan::mtbdd_false) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) sylvan::mtbdd_getvalue(leaf);
            for (auto dest_state: tds->destination_set) {
                unpack_dont_care_bits_in_mtbdd_path(path_in_mtbdd_to_leaf, nfa.var_count, transitions, 0, origin_state, dest_state);
            }

            leaf = sylvan::mtbdd_enum_next(state_transitions, nfa.vars, path_in_mtbdd_to_leaf, NULL);
        }
    }

    return transitions;
}

std::string transition_to_str(const struct Transition& transition) {
    std::stringstream str_stream;
    str_stream << "Transition{.origin = " << transition.from;
    str_stream << ", .symbols = (";
    if (!transition.symbol.empty()) {
        auto bit_iter = transition.symbol.begin();
        str_stream << (int) *bit_iter++;
        for (; bit_iter != transition.symbol.end(); bit_iter++) str_stream << ", " << (int) *bit_iter;
    }
    str_stream <<  "), .destination = " << transition.to << "}";
    return str_stream.str();
}

bool Transition::operator==(const Transition& other) const {
    return from == other.from && to == other.to && symbol == other.symbol;
}

/*
The MTBDD pad closure works in two steps:
1. Build a 'frontier' - an MTBDD mapping alphabet symbols to all states that can reach
   an accepting state via the corresponding symbol
2. Use the frontier to augment all transition MTBDDs in case that no final state
   can be reached directly, but it can be reached via a series of transitions along
   certain alhabet symbol (decided using the frontier)
 */
void NFA::perform_pad_closure(u64 leaf_type_id) {
    if (states.empty()) return;

    LACE_ME;
    using namespace sylvan; // @Cleanup: Remove this once we migrate to the new Sylvan version

    Transition_Destination_Set frontier_init;
    frontier_init.destination_set = std::set<State>(this->final_states);
    MTBDD frontier = sylvan::mtbdd_makeleaf(leaf_type_id, reinterpret_cast<u64>(&frontier_init));

    bool was_frontier_modified = true;
    while (was_frontier_modified) {
        MTBDD new_frontier = frontier;
        for (auto& [origin, state_mtbdd]: transitions) {
            new_frontier = mtbdd_applyp(state_mtbdd, new_frontier, origin, TASK(build_pad_closure_fronier_op), AMAYA_EXTEND_FRONTIER_OP_ID );
        }
        was_frontier_modified = (new_frontier != frontier);
        frontier = new_frontier;
    }


    State new_final_state = *states.rbegin() + 1;
    Pad_Closure_Info2 pad_closure_info = {.new_final_state = new_final_state, .final_states = &final_states};

    bool was_any_transition_added = false;
    for (auto& state_transition_pair: transitions) {
        pad_closure_info.origin_state = state_transition_pair.first;
        MTBDD old_mtbdd = state_transition_pair.second;
        state_transition_pair.second = mtbdd_applyp(state_transition_pair.second,
                                                    frontier,
                                                    reinterpret_cast<u64>(&pad_closure_info),
                                                    TASK(add_pad_transitions_op), AMAYA_ADD_PAD_TRANSITIONS_OP_ID);
        if (old_mtbdd != state_transition_pair.second) was_any_transition_added = true;
    }

    if (was_any_transition_added) {
        states.insert(new_final_state);
        final_states.insert(new_final_state);
    }
}

void collect_reachable_states_from_mtbdd_leaves(MTBDD root, std::set<State>& output) {
    // @Todo: This should be implemented as a proper sylvan (standalong) task utilizing cache
    if (root == sylvan::mtbdd_false) {
        return;
    }

    if (sylvan::mtbdd_isleaf(root)) {
        auto raw_leaf_value = sylvan::mtbdd_getvalue(root);
        auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(raw_leaf_value);

        output.insert(leaf_contents->destination_set.begin(), leaf_contents->destination_set.end());
    } else {
        collect_reachable_states_from_mtbdd_leaves(sylvan::mtbdd_getlow(root), output);
        collect_reachable_states_from_mtbdd_leaves(sylvan::mtbdd_gethigh(root), output);
    }
}

std::set<State> NFA::get_state_post(State state) {
    std::set<State> reachable_states;

#if 0
    MTBDD transition_mtbdd = transitions[state];

    BDDSET relevant_vars = this->vars;
    u8 symbol[this->var_count];

    MTBDD leaf = sylvan::mtbdd_enum_first(transition_mtbdd, relevant_vars, symbol, NULL);
    while (leaf != sylvan::mtbdd_false) {
        auto raw_leaf_value = sylvan::mtbdd_getvalue(leaf);
        auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(raw_leaf_value);

        reachable_states.insert(leaf_contents->destination_set.begin(), leaf_contents->destination_set.end());

        leaf = sylvan::mtbdd_enum_next(transition_mtbdd, relevant_vars, symbol, NULL);
    }
#else
    collect_reachable_states_from_mtbdd_leaves(transitions[state], reachable_states);
#endif

    return reachable_states;
}

void NFA::remove_states(std::set<State>& states_to_remove) {
    LACE_ME;
    for (auto transitions_it = transitions.begin(); transitions_it != transitions.end(); ) {
        auto [state, state_transitions] = *transitions_it;

        if (states_to_remove.contains(state)) {
            transitions_it = transitions.erase(transitions_it);
        } else {
            using namespace sylvan;
            MTBDD new_state_transitions = mtbdd_uapply(state_transitions, TASK(remove_states2_op), reinterpret_cast<u64>(&states_to_remove));
            sylvan::mtbdd_deref(state_transitions);
            sylvan::mtbdd_ref(new_state_transitions);
            transitions[state] = new_state_transitions;
            ++transitions_it;
        }
    }

    in_place_set_difference(states, states_to_remove);
    in_place_set_difference(final_states, states_to_remove);
    in_place_set_difference(initial_states, states_to_remove);
}

MTBDD compute_states_reaching_set_by_repeated_symbol(NFA& nfa, std::set<State>& states_to_reach) {
    LACE_ME;

    Transition_Destination_Set frontier_init;
    frontier_init.destination_set = std::set<State>(states_to_reach);

    MTBDD frontier = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, reinterpret_cast<u64>(&frontier_init));
    sylvan::mtbdd_refs_push(frontier);

    bool was_frontier_modified = true;
    while (was_frontier_modified) {
        MTBDD new_frontier = frontier;
        sylvan::mtbdd_refs_push(new_frontier);

        for (auto& [origin, state_mtbdd]: nfa.transitions) {
            using namespace sylvan;
            new_frontier = mtbdd_applyp(state_mtbdd, new_frontier, origin, TASK(build_pad_closure_fronier_op), AMAYA_EXTEND_FRONTIER_OP_ID);
            sylvan::mtbdd_refs_pop(1);
            sylvan::mtbdd_refs_push(new_frontier);
        }

        was_frontier_modified = (new_frontier != frontier);

        frontier = new_frontier;
        sylvan::mtbdd_refs_pop(1);
        sylvan::mtbdd_refs_push(frontier);
    }

    return frontier;
}


std::set<State> compute_states_reaching_set(NFA nfa, std::set<State>& states_to_reach) {
    std::set<State> reaching_states(states_to_reach);
    std::set<State> potential_states;
    std::set_difference(nfa.states.begin(), nfa.states.end(),
                        states_to_reach.begin(), states_to_reach.end(),
                        std::inserter(potential_states, potential_states.end()));

    // @Optimize: A hash table would not be needed here, if we are sure that nfa states always span integers in [0, states.size()-1]
    std::unordered_map<State, std::set<State>> state_posts;
    for (State state: potential_states) {
        state_posts[state] = nfa.get_state_post(state);;
    }

    bool was_fixed_point_found = false;
    u64 reaching_set_size = reaching_states.size();  // It is sufficient to just check whether the state partition is still growing

    u64 total_post_size = 0;
    for (auto [state, post]: state_posts) {
        total_post_size += post.size();
    }

    // std::cout << "Average post size: " << static_cast<double>(total_post_size) / static_cast<double>(state_posts.size()) << std::endl;

    while (!was_fixed_point_found) {
        for (auto potential_states_it = potential_states.begin(); potential_states_it != potential_states.end(); ) {
            auto state = *potential_states_it;
            auto& state_post = state_posts[state];

            if (!is_set_intersection_empty(state_post.begin(), state_post.rbegin(), state_post.end(), reaching_states.begin(), reaching_states.rbegin(), reaching_states.end())) {
                reaching_states.insert(state);
                potential_states_it = potential_states.erase(potential_states_it);
            } else {
                ++potential_states_it;
            }
        }
        was_fixed_point_found = reaching_set_size == reaching_states.size();
        reaching_set_size = reaching_states.size();
    }

    return reaching_states;
}


void remove_nonfinishing_states(NFA& nfa) {
    if (nfa.states.empty()) return;

    auto states_reaching_final = compute_states_reaching_set(nfa, nfa.final_states);

    std::set<State> states_to_remove;
    std::set_difference(nfa.states.begin(), nfa.states.end(),
                        states_reaching_final.begin(), states_reaching_final.end(),
                        std::inserter(states_to_remove, states_to_remove.begin()));

    nfa.remove_states(states_to_remove);

    if (nfa.states.empty()) {
        nfa.states.insert(1);
        nfa.initial_states.insert(1);

        u8 symbol[nfa.var_count];
        for (u64 i = 0; i < nfa.var_count; i++) {
            symbol[i] = 2;
        }
        nfa.add_transition(1, 1, symbol);
    }
}


void setup_intersection_info_for_postless_pruning(Intersection_Info2& info, NFA& left, NFA& right) {
#if INTERSECTION_DETECT_STATES_WITH_NO_POST
    for (auto& final_state: left.final_states) {
        if (left.transitions[final_state] == sylvan::mtbdd_false) {
            info.left_states_without_post.push_back(final_state);
        }
    }
    for (auto& final_state: right.final_states) {
        if (right.transitions[final_state] == sylvan::mtbdd_false) {
            info.right_states_without_post.push_back(final_state);
        }
    }

    info.left_final_states  = &left.final_states;
    info.right_final_states = &right.final_states;
#endif
}


NFA compute_nfa_intersection(NFA& left, NFA& right) {
    LACE_ME;
    typedef std::pair<State, State> Product_State;

    std::vector<Intersection_Discovery> work_queue;
    Intersection_Info2 intersection_info = {.seen_products = {}, .work_queue = work_queue};

    setup_intersection_info_for_postless_pruning(intersection_info, left, right);

    sylvan::BDDSET intersection_vars = left.vars;

    const u64 right_vars_cnt = sylvan::mtbdd_set_count(right.vars);
    u32 right_vars[right_vars_cnt];
    sylvan::mtbdd_set_to_array(right.vars, right_vars);
    for (u64 i = 0; i < right_vars_cnt; i++) {
        intersection_vars = sylvan::mtbdd_set_add(intersection_vars, right_vars[i]);
    }

    NFA intersection_nfa = {.vars = intersection_vars, .var_count = sylvan::mtbdd_set_count(intersection_vars)};

    work_queue.reserve(left.initial_states.size() * right.initial_states.size());

    for (auto& left_init_state: left.initial_states) {
        for (auto& right_init_state: right.initial_states) {
            auto product_state = std::make_pair(left_init_state, right_init_state);
            auto handle = static_cast<State>(intersection_info.seen_products.size());
            intersection_info.seen_products.emplace(product_state, handle);
            intersection_nfa.initial_states.emplace(handle);

            work_queue.push_back({.left = product_state.first, .right = product_state.second, .handle = handle});
        }
    }

    while (!work_queue.empty()) {
        auto explored_product = work_queue.back();
        work_queue.pop_back();

        intersection_nfa.states.insert(explored_product.handle);

        MTBDD left_mtbdd = left.transitions[explored_product.left];
        MTBDD right_mtbdd = right.transitions[explored_product.right];

        if (left.final_states.contains(explored_product.left) && right.final_states.contains(explored_product.right)) {
            intersection_nfa.final_states.insert(explored_product.handle);
        }

        using namespace sylvan;
        sylvan::MTBDD product_mtbdd = mtbdd_applyp(left_mtbdd, right_mtbdd, reinterpret_cast<u64>(&intersection_info),
                                                   TASK(transitions_intersection2_op), AMAYA_TRANSITIONS_INTERSECT_OP_ID);
        sylvan::mtbdd_ref(product_mtbdd);
        intersection_nfa.transitions[explored_product.handle] = product_mtbdd;
    }

#if INTERSECTION_REMOVE_NONFINISHING_STATES
    remove_nonfinishing_states(intersection_nfa);
#endif

    return intersection_nfa;
}

NFA determinize_nfa(NFA& nfa) {
    LACE_ME;
    NFA result = {.vars = nfa.vars, .var_count = nfa.var_count};

    std::vector<std::pair<const Macrostate, State>*> work_queue;
    Determinization_Context ctx = {.known_macrostates = {}, .work_queue = work_queue};
    {
        Macrostate initial_state(nfa.initial_states.begin(), nfa.initial_states.end());
        State initial_handle = 0;
        auto emplace_status = ctx.known_macrostates.emplace(initial_state, 0);
        work_queue.push_back(&(*emplace_status.first));
    }

    u8 path_in_transitions_mtbdd[nfa.var_count];

    while (!work_queue.empty()) {
        auto current_macrostate_entry = work_queue.back();
        work_queue.pop_back();

        auto& macrostate = current_macrostate_entry->first;
        auto& handle = current_macrostate_entry->second;

        result.states.insert(handle);
        if (!is_set_intersection_empty(macrostate.begin(), macrostate.rbegin(), macrostate.end(),
                                       nfa.final_states.begin(), nfa.final_states.rbegin(), nfa.final_states.end())) {
            nfa.final_states.insert(handle);
        }

        sylvan::MTBDD macrostate_transitions = sylvan::mtbdd_false;
        sylvan::mtbdd_refs_push(macrostate_transitions);
        for (auto& macrostate_member: macrostate) {
            auto& member_mtbdd = nfa.transitions[macrostate_member];

            using namespace sylvan;
            sylvan::MTBDD new_macrostate_transitions = mtbdd_applyp(macrostate_transitions, member_mtbdd, 0u,
                                                       TASK(transitions_union_op), AMAYA_UNION_OP_ID);
            sylvan::mtbdd_refs_pop(1);
            sylvan::mtbdd_refs_push(new_macrostate_transitions);
            macrostate_transitions = new_macrostate_transitions;
        }

        using namespace sylvan;
        MTBDD transition_mtbdd = mtbdd_uapply(macrostate_transitions,
                                             TASK(replace_macrostates_with_handles_op),
                                             reinterpret_cast<u64>(&(ctx)));
        sylvan::mtbdd_ref(transition_mtbdd);
        sylvan::mtbdd_refs_pop(1);
        result.transitions.emplace(handle, transition_mtbdd);
    }

    return nfa;
}

template<typename T>
std::ostream& operator<<(std::ostream& output, const std::set<T>& set) {
    output << "{";
    if (!set.empty()) {
        auto set_it = set.begin();
        output << *set_it;
        ++set_it;
        for (; set_it != set.end(); ++set_it) output << ", " << *set_it;
    }
    output << "}";
    return output;
}

std::ostream& operator<<(std::ostream& output, const NFA& nfa) {
    output << "NFA{" << std::endl;
    output << "  states: "  << nfa.states << std::endl;
    output << "  final_states: "  << nfa.final_states << std::endl;
    output << "  initial_states: "  << nfa.initial_states << std::endl;
    output << "  var_count: "  << nfa.var_count << std::endl;
    output << "  vars: "  << nfa.vars << std::endl;
    output << "}";
    return output;
}
