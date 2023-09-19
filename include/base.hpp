#ifndef AMAYA_BASE_H
#define AMAYA_BASE_H

#define DEBUG 1

#define INTERSECTION_DETECT_STATES_WITH_NO_POST 0
#define INTERSECTION_REMOVE_NONFINISHING_STATES 0

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>

#include <inttypes.h>

#include <sylvan.h>

#if DEBUG
#define PRINT_DEBUG(it) do { std::cerr << it << std::endl; } while (0)
#define PRINTF_DEBUG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define PRINT_DEBUG(it)
#define PRINTF_DEBUG(...)
#endif


using std::set;
using std::unordered_map;

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef int64_t  s64;
typedef int32_t  s32;

typedef int64_t State;

typedef std::vector<State> Macrostate;

struct Transition {
    State from;
    State to;
    std::vector<uint8_t> symbol;

    bool operator==(const Transition& other) const;
};
std::ostream& operator<<(std::ostream& output, const Transition& transition);

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
    u64 var_count;

    NFA(sylvan::BDDSET vars = sylvan::mtbdd_set_empty(),
        u64 var_count = 0,
        std::set<State>&& states={},
        std::set<State>&& final_states={},
        std::set<State>&& initial_states={}) : vars(vars), var_count(var_count), states(states), final_states(final_states), initial_states(initial_states), transitions({}) {
        sylvan::mtbdd_ref(vars);
    };
    ~NFA() {
        sylvan::mtbdd_deref(vars);
    }
    NFA(const NFA& other) : vars(other.vars), var_count(other.var_count), states(other.states), final_states(other.final_states), initial_states(other.initial_states), transitions(other.transitions)  {
        sylvan::mtbdd_ref(vars);
    }


    // @Cleanup: Factor out the sylvan global configuration into a context struct
    void perform_pad_closure();

    void add_transition(State from, State to, u64 symbol, u64 quantified_bits_mask);
    void add_transition(State from, State to, u8* symbol);
    void add_transition(State from, State to, std::vector<u8>&& symbol);
    void add_universal_transition(State from, State to);

    std::string show_transitions() const;
    std::vector<Transition> get_symbolic_transitions_for_state(State state) const;

    std::set<State> get_state_post(State state);

    void remove_states(std::set<State>& states_to_remove);
};


NFA compute_nfa_intersection(NFA& left, NFA& right);
void remove_nonfinishing_states(NFA& nfa);
std::set<State> compute_states_reaching_set(NFA nfa, std::set<State>& states_to_reach);
NFA determinize_nfa(NFA& nfa);
NFA minimize_hopcroft(NFA& nfa);


std::vector<struct Transition> nfa_unpack_transitions(struct NFA& nfa);
std::string transition_to_str(const struct Transition& transition);
bool transition_is_same_as(const struct Transition& transition_a, const struct Transition& transition_b);

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

template<typename T>
std::ostream& operator<<(std::ostream& output, const std::vector<T>& arr) {
    output << "[";
    if (!arr.empty()) {
        auto arr_it = arr.begin();
        output << *arr_it;
        ++arr_it;
        for (; arr_it != arr.end(); ++arr_it) output << ", " << *arr_it;
    }
    output << "]";
    return output;
}

std::ostream& operator<<(std::ostream& output, const NFA& nfa);

template <>
struct std::hash<Macrostate> {
    std::size_t operator() (const Macrostate& macrostate) const {
        std::size_t hash = 0;
        for (auto state: macrostate) {
            std::size_t atom_val_hash = std::hash<State>{}(state);
            hash = hash + 0x9e3779b9 + (atom_val_hash << 6) + (atom_val_hash >> 2);
        }
        return hash;
    }
};

template <typename InputIterator1, typename ReverseInputIterator1, typename InputIterator2, typename ReverseInputIterator2>
bool is_set_intersection_empty(InputIterator1 left_begin, ReverseInputIterator1 left_rbegin, InputIterator1 left_end,
                               InputIterator2 right_begin, ReverseInputIterator2 right_rbegin, InputIterator2 right_end) {

    if (left_begin == left_end || right_begin == right_end) return true;
    if (*left_begin > *right_rbegin || *right_begin > *left_rbegin) return true;

    while (left_begin != left_end && right_begin != right_end) {
        if (*left_begin == *right_begin) return false;
        else if (*left_begin < *right_begin) ++left_begin;
        else ++right_begin;
    }

    return true;
}

template <typename T>
bool is_set_intersection_empty(const std::set<T>& left, const std::set<T>& right) {
    return is_set_intersection_empty(left.begin(), left.rbegin(), left.end(), right.begin(), right.rbegin(), right.end());
}

template <typename T>
void in_place_set_difference(std::set<T>& left, const std::set<T>& right) {
    auto left_it = left.begin();
    auto right_it = right.begin();

    while (left_it != left.end() && right_it != right.end()) {
        if (*left_it < *right_it) {
            ++left_it;
        } else if (*right_it < *left_it) {
            ++right_it;
        } else {
            left_it = left.erase(left_it);
            ++right_it;
        }
    }
}

#endif
