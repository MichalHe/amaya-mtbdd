#include "../include/base.hpp"
#include "../include/hopcroft_leaf.hpp"
#include "../include/custom_leaf.hpp"

#include <functional>
#include <sstream>
#include <cstring>
#include <iostream>

#include <sylvan.h>

uint64_t mtbdd_leaf_type_hopcroft;

Hopcroft_Leaf_Contents::Hopcroft_Leaf_Contents()
{
    this->destination_to_origin_states = std::map<State, std::set<State>>();
}


Hopcroft_Leaf_Contents::Hopcroft_Leaf_Contents(Hopcroft_Leaf_Contents& orignal)
{
    this->destination_to_origin_states = std::map<State, std::set<State>>(orignal.destination_to_origin_states);
}

Hopcroft_Leaf_Contents::~Hopcroft_Leaf_Contents() {}

std::string Hopcroft_Leaf_Contents::to_str()
{
    std::stringstream string_stream;

    uint64_t destination_count = 0;
    string_stream << "[";
    for (auto dest_to_origin_states_pair : this->destination_to_origin_states) {
        string_stream << dest_to_origin_states_pair.first
                      << " <- {";
        uint64_t state_count = 0;
        for (auto origin_state : dest_to_origin_states_pair.second) {
            string_stream << origin_state;
            if (state_count < dest_to_origin_states_pair.second.size() - 1)
                string_stream << ", ";
            state_count++;
        }
        string_stream << "}";

        if (destination_count < this->destination_to_origin_states.size() - 1)
            string_stream << ", ";
        destination_count++;
    }
    string_stream << "]";

    std::string contents_str = string_stream.str();
    return contents_str;
}

void hopcroft_leaf_create(uint64_t* hopcroft_leaf_contents_ptr_param)
{
    auto hopcroft_leaf_contents_ptr = (Hopcroft_Leaf_Contents**)(hopcroft_leaf_contents_ptr_param);
    auto original_leaf_contents = *hopcroft_leaf_contents_ptr;
    auto leaf_contents_copy = new Hopcroft_Leaf_Contents(*original_leaf_contents);

    *hopcroft_leaf_contents_ptr = leaf_contents_copy;
}

uint64_t hopcroft_leaf_hash(const uint64_t contents_ptr, const uint64_t seed)
{
    uint64_t leaf_hash = seed;
    auto hopcroft_leaf_contents = (Hopcroft_Leaf_Contents*) contents_ptr;

    std::hash<State> hash_fn;

    for (auto pair : hopcroft_leaf_contents->destination_to_origin_states) {
        auto dest_state = pair.first;
        auto origin_state_set = pair.second;

        leaf_hash ^= hash_fn(dest_state) + (seed<<6) + 0x9e3779b9 + (seed>>2);
        for (auto origin_state: origin_state_set) {
            // TODO: Add citation (boost hash combine)
            leaf_hash ^= hash_fn(origin_state) + (seed<<6) + 0x9e3779b9 + (seed>>2);
        }
    }

    return leaf_hash;
}


char* hopcroft_leaf_to_str(int is_leaf_complemented, uint64_t leaf_contents_untyped_ptr, char* buffer, size_t buffer_size)
{
    (void) is_leaf_complemented;
    auto hopcroft_leaf_contents_ptr = (Hopcroft_Leaf_Contents*) leaf_contents_untyped_ptr;

    auto hopcroft_leaf_string = hopcroft_leaf_contents_ptr->to_str();

    char* output_buffer;

    if (buffer_size <= hopcroft_leaf_string.size() + 1) { // +1 for the nullbyte
        output_buffer = (char*) malloc(sizeof(char) * (hopcroft_leaf_string.size() + 1));
    } else {
        output_buffer = buffer; // Use the preallocated buffer as it has sufficient capacity
    }
    std::memcpy(output_buffer, hopcroft_leaf_string.c_str(), hopcroft_leaf_string.size());
    output_buffer[hopcroft_leaf_string.size()] = '\0';
    return output_buffer;
}

void hopcroft_leaf_destroy(uint64_t hl_untyped_ptr)
{
    auto hl_ptr = (Hopcroft_Leaf_Contents*) hl_untyped_ptr;
    delete hl_ptr;
}

int hopcroft_leaf_equals(uint64_t hl_a_untyped_ptr, uint64_t hl_b_untyped_ptr)
{
    auto hl_a_ptr = (Hopcroft_Leaf_Contents*) hl_a_untyped_ptr;
    auto hl_b_ptr = (Hopcroft_Leaf_Contents*) hl_b_untyped_ptr;

    return hl_a_ptr->destination_to_origin_states == hl_b_ptr->destination_to_origin_states;
}


TASK_IMPL_2(sylvan::MTBDD, hopcroft_leaf_from_normal, sylvan::MTBDD, mtbdd, uint64_t, origin_state_untyped)
{
    if (mtbdd == sylvan::mtbdd_false) {
        return sylvan::mtbdd_false;
    }

    if (!sylvan::mtbdd_isleaf(mtbdd)) {
        return sylvan::mtbdd_invalid;
    }

    assert(sylvan::mtbdd_gettype(mtbdd) == mtbdd_leaf_type_set);
    Transition_Destination_Set* tds = (Transition_Destination_Set*) sylvan::mtbdd_getvalue(mtbdd);

    State origin_state = (State) origin_state_untyped;
    Hopcroft_Leaf_Contents leaf_contents;

    std::map<State, std::set<State>> dest_state_to_origin_states;
    for (State dest_state: tds->destination_set) {
         dest_state_to_origin_states.insert({dest_state, {origin_state}});
    }

    leaf_contents.destination_to_origin_states = dest_state_to_origin_states;

    return sylvan::mtbdd_makeleaf(mtbdd_leaf_type_hopcroft, (uint64_t) &leaf_contents);
}

TASK_IMPL_3(sylvan::MTBDD,
            hopcroft_leaf_union,
            sylvan::MTBDD*, l_mtbdd_ptr,
            sylvan::MTBDD*, r_mtbdd_ptr,
            uint64_t, origin_state_untyped)
{
    sylvan::MTBDD l_mtbdd = *l_mtbdd_ptr;
    sylvan::MTBDD r_mtbdd = *r_mtbdd_ptr;

    if (l_mtbdd == sylvan::mtbdd_false) {
        return r_mtbdd;
    }

    if (r_mtbdd == sylvan::mtbdd_false) {
        return l_mtbdd; // Left is not False
    }

    if (!(sylvan::mtbdd_isleaf(l_mtbdd) && sylvan::mtbdd_isleaf(r_mtbdd))) {
        return sylvan::mtbdd_invalid;
    }

    // Both mtbdds are leaves
    Hopcroft_Leaf_Contents* l_hl_contents = (Hopcroft_Leaf_Contents*) sylvan::mtbdd_getvalue(l_mtbdd);
    Hopcroft_Leaf_Contents* r_hl_contents = (Hopcroft_Leaf_Contents*) sylvan::mtbdd_getvalue(r_mtbdd);

    Hopcroft_Leaf_Contents new_leaf_contents;
    new_leaf_contents.destination_to_origin_states = std::map<State, std::set<State>>(l_hl_contents->destination_to_origin_states);
    auto& new_dto = new_leaf_contents.destination_to_origin_states;  // For brewity

    // Merge the two leaf contents
    for (auto dest_to_origin_pair: r_hl_contents->destination_to_origin_states) {
        const State& dest_state = dest_to_origin_pair.first;
        std::set<State>& origin_states = dest_to_origin_pair.second;

        // Check whether the leaf contents already contains the destination, if so merge it
        if (new_dto.find(dest_to_origin_pair.first) != new_dto.end()) {
            // Just unite the two sets
            new_dto[dest_state].insert(origin_states.begin(), origin_states.end());
        } else {
            new_dto[dest_state] = std::set<State>(origin_states);  // Make copy
        }
    }

    return sylvan::mtbdd_makeleaf(mtbdd_leaf_type_hopcroft, (uint64_t) &new_leaf_contents);
}
