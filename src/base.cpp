#include "../include/base.hpp"

#include <iostream>
#include <set>
#include <vector>
#include <sstream>
#include <string>

using sylvan::MTBDD;
using sylvan::mtbdd_support_CALL;

Transition_Destination_Set::Transition_Destination_Set()
{
    this->destination_set = NULL;
}

Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other)
{
    // NOTE: Copy is called when a new leaf is created
    this->destination_set = new std::set<State>(*other.destination_set);
}

Transition_Destination_Set::Transition_Destination_Set(std::set<State> *destination_set)
{
    this->destination_set = destination_set;
}

Transition_Destination_Set::~Transition_Destination_Set()
{
    if (this->destination_set != NULL)
    {
        delete this->destination_set;
    }
}

void Transition_Destination_Set::print_dest_states()
{
    int cnt = 1;
    std::cout << "{";
    for (auto state : *this->destination_set)
    {
        std::cout << state;
        if (cnt < this->destination_set->size())
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
        std::vector<uint8_t> symbols;
        symbols.reserve(path_size);
        for (uint64_t i = 0; i < path_size; i++) symbols.push_back(mtbdd_path[i]);
        struct Transition transition = {
            .origin = origin,
            .destination = destination,
            .symbols = symbols
        };
        unpacked_transitions.push_back(transition);
    }
}


// @Unoptimized
std::vector<struct Transition> nfa_unpack_transitions(struct NFA& nfa)
{
    LACE_ME;

    std::vector<struct Transition> transitions;

    uint8_t* path_in_mtbdd_to_leaf = (uint8_t*) malloc(sizeof(uint8_t) * nfa.var_count);
    assert(path_in_mtbdd_to_leaf != NULL);

    for (auto origin_state: nfa.states) {
        MTBDD    state_transitions = nfa.transitions[origin_state];
        MTBDD    leaf         = sylvan::mtbdd_enum_first(state_transitions, nfa.vars, path_in_mtbdd_to_leaf, NULL);

        while (leaf != sylvan::mtbdd_false) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) sylvan::mtbdd_getvalue(leaf);
            for (auto dest_state: *tds->destination_set) {
                unpack_dont_care_bits_in_mtbdd_path(path_in_mtbdd_to_leaf, nfa.var_count, transitions, 0, origin_state, dest_state);
            }

            leaf = sylvan::mtbdd_enum_next(state_transitions, nfa.vars, path_in_mtbdd_to_leaf, NULL);
        }
    }

    return transitions;
}

std::string transition_to_str(const struct Transition& transition) {
    std::stringstream str_stream;
    str_stream << "Transition{.origin = " << transition.origin;
    str_stream << ", .symbols = (";
    if (!transition.symbols.empty()) {
        auto bit_iter = transition.symbols.begin();
        str_stream << (int) *bit_iter++;
        for (; bit_iter != transition.symbols.end(); bit_iter++) str_stream << ", " << (int) *bit_iter;
    }
    str_stream <<  "), .destination = " << transition.destination << "}";
    return str_stream.str();
}

bool transition_is_same_as(const struct Transition& transition_a, const struct Transition& transition_b)
{
    return transition_a.origin == transition_b.origin &&
           transition_a.destination == transition_b.destination &&
           transition_a.symbols == transition_b.symbols;
}
