#ifndef AMAYA_BASE_H
#define AMAYA_BASE_H

#include <set>
#include <inttypes.h>

typedef int64_t State;

class Transition_Destination_Set {
public:
	std::set<State>* destination_set;
	uint32_t automaton_id;

	Transition_Destination_Set();
	Transition_Destination_Set(const Transition_Destination_Set &other);
	Transition_Destination_Set(uint32_t automaton_id, std::set<State>* destination_set);
	~Transition_Destination_Set();
	void print_dest_states();
};

#endif
