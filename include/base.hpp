#ifndef AMAYA_BASE_H
#define AMAYA_BASE_H

#define DEBUG 1

#include <vector>
#include <set>
#include <unordered_map>
#include <string>

#include <inttypes.h>

#include <sylvan.h>

using std::set;
using std::unordered_map;

using sylvan::MTBDD;
using sylvan::BDDSET;

typedef int64_t State;

class Transition_Destination_Set {
public:
	std::set<State>* destination_set;

	Transition_Destination_Set();
	Transition_Destination_Set(const Transition_Destination_Set &other);
	Transition_Destination_Set(std::set<State>* destination_set);
	~Transition_Destination_Set();
	void print_dest_states();
};

struct NFA {
    set<State> states;
    set<State> final_states;
    set<State> initial_states;
    unordered_map<State, sylvan::MTBDD> transitions;

    sylvan::BDDSET vars;
    uint64_t var_count;
};

struct Transition {
    State origin;
    State destination;
    std::vector<uint8_t> symbols;
};


std::vector<struct Transition> nfa_unpack_transitions(struct NFA& nfa);
std::string transition_to_str(const struct Transition& transition);
bool transition_is_same_as(const struct Transition& transition_a, const struct Transition& transition_b);

#endif
