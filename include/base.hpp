#ifndef AMAYA_BASE_H
#define AMAYA_BASE_H

#define DEBUG 0

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

typedef uint64_t u64;
typedef uint8_t  u8;
typedef int64_t  s64;

struct Transition {
    State from;
    State to;
    std::vector<uint8_t> symbol;

    bool operator==(const Transition& other) const;
};

// @Todo: Make this an ordinary struct
class Transition_Destination_Set {
public:
    set<State> destination_set;

    Transition_Destination_Set() {};
    Transition_Destination_Set(const Transition_Destination_Set &other);
    Transition_Destination_Set(set<State>& destination_set);
    Transition_Destination_Set(set<State>&& destination_set) : destination_set(destination_set) {};
    void print_dest_states();
};

struct NFA {
    set<State> states;
    set<State> final_states;
    set<State> initial_states;
    unordered_map<State, sylvan::MTBDD> transitions;

    sylvan::BDDSET vars;
    uint64_t var_count;

    // @Cleanup: Factor out the sylvan global configuration into a context struct
    void perform_pad_closure(u64 leaf_type_id);

    void add_transition(State from, State to, u64 symbol, u64 quantified_bits_mask);
    void add_transition(State from, State to, u8* symbol);

    std::string show_transitions() const;
    std::vector<Transition> get_symbolic_transitions_for_state(State state) const;
};


std::vector<struct Transition> nfa_unpack_transitions(struct NFA& nfa);
std::string transition_to_str(const struct Transition& transition);
bool transition_is_same_as(const struct Transition& transition_a, const struct Transition& transition_b);

#endif
