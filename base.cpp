#include "base.hpp"

#include <iostream>
#include <set>

Transition_Destination_Set::Transition_Destination_Set()
{
    this->destination_set = NULL;
    this->automaton_id = 0;
}

Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other)
{
    this->destination_set = new std::set<State>(*other.destination_set); // Do copy
    this->automaton_id = other.automaton_id;
}

Transition_Destination_Set::Transition_Destination_Set(uint32_t automaton_id, std::set<State> *destination_set)
{
    this->automaton_id = automaton_id;
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
