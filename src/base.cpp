#include "../include/base.hpp"
#include "../include/custom_leaf.hpp"
#include "../include/operations.hpp"

#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include <sstream>
#include <string>

using sylvan::MTBDD;
using sylvan::mtbdd_support_CALL;

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
        for (uint64_t i = 0; i < path_size; i++) symbol.push_back(mtbdd_path[i]);

        Transition transition = {.from = origin, .to = destination, .symbol = symbol};
        unpacked_transitions.push_back(transition);
    }
}


std::vector<Transition> NFA::get_symbolic_transitions_for_state(State state) const {
    std::vector<Transition> symbolic_transitions;

    auto state_mtbdd_it = transitions.find(state);

    // @Robustness: Keeping the assert as we should access the transitions only when iterating trough states
    assert(state_mtbdd_it != transitions.end());

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

void NFA::add_transition(State from, State to, u8* rich_symbol) {
    Transition_Destination_Set leaf_contents({to});
    sylvan::MTBDD leaf = make_set_leaf(&leaf_contents);

    sylvan::MTBDD transition_mtbdd = sylvan::mtbdd_cube(this->vars, rich_symbol, leaf);

    LACE_ME;
    auto [present_key_val_it, was_inserted] = transitions.emplace(from, transition_mtbdd);
    if (!was_inserted) {
        auto existing_transitions = present_key_val_it->second;

        using namespace sylvan; // Pull the entire sylvan namespace because mtbdd_applyp is a macro
        sylvan::MTBDD updated_mtbdd = mtbdd_applyp(existing_transitions, transition_mtbdd, 0u, TASK(transitions_union_op), AMAYA_UNION_OP_ID);
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
