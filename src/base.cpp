#include "../include/base.hpp"

#include <iostream>
#include <set>

Transition_Destination_Set::Transition_Destination_Set()
{
    this->destination_set = NULL;
}

Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other)
{
    this->destination_set = new std::set<State>(*other.destination_set); // Do copy
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
