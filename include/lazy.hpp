#ifndef AMAYA_LAZY_H
#define AMAYA_LAZY_H

#include "base.hpp"
#include "custom_leaf.hpp"
#include "operations.hpp"

#include <bitset>
#include <list>
#include <cassert>
#include <optional>
#include <cstring>
#include <iostream>
#include <string>
#include <map>
#include <unordered_set>

using std::optional;
using std::vector;
using std::unordered_set;
using std::map;
using std::list;

template <typename T>
struct Sized_Array {
    T* items;
    u64 size;

    struct Iterator {
        const Sized_Array<T>* arr;
        u64                   idx;

        using value_type = T;
        using pointer    = T*;
        using reference  = T&;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;

        Iterator():                    arr(nullptr), idx(0) {}   
        Iterator(const Sized_Array<T>* arr, u64 start_idx): arr(arr), idx(start_idx) {}

        reference  operator*() {
            return (arr->items)[idx];
        }
        pointer operator->() {
            return &(arr->items[idx]);
        }
        const pointer operator->() const {
            return &(arr->items[idx]);
        }
        reference  operator[](int offset) {
            return arr->items[idx + offset];
        }
        // const reference operator*() const {
        //     return arr->items[idx];
        // }
        // const reference operator[](int offset) const {
        //     return arr->items[idx + offset];
        // }

        Iterator& operator++() {
            ++idx;
            return *this;
        }
        Iterator operator++(int) {
            Iterator r(*this);
            ++idx;
            return r;
        }

        Iterator& operator--() {
            --idx;
            return *this;
        }
        Iterator operator--(int) {
            Iterator r(*this);
            --idx;
            return r;
        }

        Iterator& operator+=(int offset) {
            idx += offset;
            return *this;
        }
        Iterator& operator-=(int offset) {
            idx -= offset;
            return *this;
        }

        Iterator operator+(int offset) const {
            Iterator r(*this);
            return r += offset;
        }
        Iterator operator-(int offset) const {
            Iterator r(*this);
            return r -= offset;
        }

        difference_type operator-(Iterator const& r) const {
            return idx - r.idx;
        }

        bool operator<(Iterator const& r)  const {
            return idx <  r.idx;
        }
        bool operator<=(Iterator const& r) const {
            return idx <= r.idx;
        }
        bool operator>(Iterator const& r)  const {
            return idx >  r.idx;
        }
        bool operator>=(Iterator const& r) const {
            return idx >= r.idx;
        }
        bool operator!=(const Iterator &r) const {
            return idx != r.idx;
        }
        bool operator==(const Iterator &r) const {
            return idx == r.idx;
        }
    };

    Iterator begin() const {
        return Iterator(this, 0);
    };

    Iterator end() const {
        return Iterator(this, this->size);
    };
};

std::size_t hash_array(Sized_Array<s64>& arr);

inline std::size_t hash_combine(std::size_t hash1, std::size_t hash2) {
    return hash1 + 0x9e3779b9 + (hash2 << 6) + (hash2 >> 2);
}

struct Congruence {
    Sized_Array<s64>* coefs;
    s64 modulus_odd;
    s64 modulus_2pow;
};

struct Inequation {
    Sized_Array<s64>* coefs;
};

struct Equation {
    Sized_Array<s64>* coefs;
};

enum Presburger_Atom_Type {
    PR_ATOM_INVALID = 0,
    PR_ATOM_INEQ = 1,
    PR_ATOM_EQ = 2,
    PR_ATOM_CONGRUENCE = 3
};

struct Presburger_Atom {
    Presburger_Atom_Type type;
    vector<s64>          coefs;
    s64                  modulus;

    Presburger_Atom() : type(Presburger_Atom_Type::PR_ATOM_INVALID) {
        coefs = {};
    }

    Presburger_Atom(Presburger_Atom_Type atom_type, const vector<s64>& coefficients, s64 modulus=0) : type(atom_type), coefs(coefficients), modulus(modulus) {}

    optional<s64> compute_post_along_sym(s64 constant, u64 symbol_bits) const;
    bool accepts_last_symbol(s64 state, u64 symbol_bits) const;

    std::string fmt_with_rhs(s64 rhs) const;
    bool operator==(const Presburger_Atom& other) const;
};

template <>
struct std::hash<Presburger_Atom> {
    std::size_t operator() (const Presburger_Atom& state) const {
        std::size_t hash = 0;
        for (auto coef: state.coefs) {
            std::size_t coef_hash = std::hash<s64>{}(coef);
            hash = hash + 0x9e3779b9 + (coef_hash << 6) + (coef_hash >> 2);
        }
        return hash;
    }
};

struct Sparse_Presburger_Atom {
    Presburger_Atom_Type type;
    u64* variables;
    s64* coefs;
    u64 var_count;
};

enum class Preferred_Var_Value_Type : unsigned {
    LOW  = 0x01,
    HIGH = 0x02,
};
constexpr Preferred_Var_Value_Type operator&(Preferred_Var_Value_Type left, Preferred_Var_Value_Type right);
constexpr Preferred_Var_Value_Type operator|(Preferred_Var_Value_Type left, Preferred_Var_Value_Type right);

struct Conjunction_State;

enum class Linear_Node_Type : u8 {
    EQ   = 0x01,
    INEQ = 0x02,
};

struct Linear_Node {
    Linear_Node_Type type;
    vector<s64>      coefs; // Cannot share coefs with original atom as it might be modified in-place during simplification
    vector<u64>      vars;
    bool             is_satisfied;
};

struct Congruence_Node {
    vector<s64> coefs;
    vector<u64> vars;
    s64         modulus_2pow;
    s64         modulus_odd;
    bool        is_satisfied;
};

struct Hard_Bound {
    bool is_present;
    u64  atom_i;
};

struct Var_Node {
    vector<u64> affected_free_vars; // Free vars affected via congruence
    Hard_Bound  hard_upper_bound;
    Hard_Bound  hard_lower_bound;
    vector<u64> upper_bounds;
    vector<u64> lower_bounds;
    vector<u64> congruences;

    bool is_hard_lower_bound(u64 atom_i) {
        return (hard_lower_bound.is_present && hard_lower_bound.atom_i == atom_i);
    }

    bool is_hard_upper_bound(u64 atom_i) {
        return (hard_upper_bound.is_present && hard_upper_bound.atom_i == atom_i);
    }

};

struct Dep_Graph {
    vector<u64>             potential_vars;
    vector<u64>             quantified_vars;
    vector<Var_Node>        var_nodes;
    vector<Linear_Node>     linear_nodes;
    vector<Congruence_Node> congruence_nodes;
};

struct Quantified_Atom_Conjunction {
    Sized_Array<Equation>*   equations;
    Sized_Array<Congruence>* congruences;
    Sized_Array<Inequation>* inequations;
    vector<u64> bound_vars;
    u64         var_count;
    Dep_Graph   dep_graph;
    bool        bottom = false;

    Quantified_Atom_Conjunction() {}
    Quantified_Atom_Conjunction(bool top) : bottom(!top) {}
    Quantified_Atom_Conjunction(const vector<Equation>& eqs, const vector<Congruence>& congruences, const vector<Inequation> ineqs, const vector<u64>& bound_vars, u64 var_count);

    bool operator==(const Quantified_Atom_Conjunction& other) const;

    u64 atom_count() const {
        return congruences->size + equations->size + inequations->size;
    }
    
    bool has_atoms() const {
        return this->atom_count() > 0;
    }

    bool is_top() const {
        return !this->has_atoms() && !bottom;
    }

    bool is_bottom() const {
        return !this->has_atoms() && bottom;
    }

    std::string fmt_with_state(Conjunction_State& state) const;
};

template <>
struct std::hash<Quantified_Atom_Conjunction> {
    std::size_t operator() (const Quantified_Atom_Conjunction& formula) const {
        std::size_t hash = formula.bottom ? 0 : 33;
        // @Optimize: The equations and inequations do not mutate during execution so we can just hash a pointer there. Only congruences need real hashing.
        for (auto& eq: *formula.equations) {
            hash = hash_combine(hash, hash_array(*eq.coefs));
        }

        for (auto& cg: *formula.congruences) {
            hash = hash_combine(hash, hash_array(*cg.coefs));
            hash = hash_combine(hash, cg.modulus_odd);
            hash = hash_combine(hash, cg.modulus_2pow);
        }
        for (auto& ineq: *formula.inequations) {
            hash = hash_combine(hash, hash_array(*ineq.coefs));
        }
        return hash;
    }
};

struct Conjunction_State {
    vector<s64> constants;

    Conjunction_State(const vector<s64>& value): constants(value) {}

    Conjunction_State(const Conjunction_State& other): constants(other.constants) {}

    void post(unordered_set<Conjunction_State>& known_states, vector<Conjunction_State>& dest);
    bool operator==(const Conjunction_State& other) const;
    bool operator<(const Conjunction_State& other) const;
};

template <>
struct std::hash<Conjunction_State> {
    std::size_t operator() (const Conjunction_State& state) const {
        std::size_t hash = 0;

        for (u64 i = 0u; i < state.constants.size(); i++) {
            std::size_t atom_val_hash = std::hash<s64>{}(state.constants[i]);
            hash = hash + 0x9e3779b9 + (atom_val_hash << 6) + (atom_val_hash >> 2);
        }

        return hash;
    }
};

std::ostream& operator<<(std::ostream& output, const Presburger_Atom& atom);
std::ostream& operator<<(std::ostream& output, const Quantified_Atom_Conjunction& formula);
std::ostream& operator<<(std::ostream& output, const Conjunction_State& atom);


struct Alphabet_Iterator {
    u64 free_bits_val;
    u64 free_bits_inc_count;

    u64 quantified_bits_val;
    u64 quantified_bits_inc_count;

    bool finished;
    bool has_more_quantif_symbols;

    u64 quantified_bits_mask;
    u64 quantified_bits_inc_limit;
    u64 free_bits_inc_limit;

    Alphabet_Iterator(u64 var_count, const vector<u64>& quantified_vars):
        free_bits_val(0u),
        free_bits_inc_count(0u),
        quantified_bits_val(0u),
        quantified_bits_inc_count(0u),
        finished(false),
        has_more_quantif_symbols(false)
    {
        quantified_bits_mask = 0u;
        for (auto quantified_var: quantified_vars) {
            quantified_bits_mask |= (1u << quantified_var);
        }

        quantified_bits_inc_limit = 1u << quantified_vars.size();
        free_bits_inc_limit = 1u << (var_count - quantified_vars.size());
    };

    u64 next_symbol() {
        quantified_bits_val = ((quantified_bits_val | ~quantified_bits_mask) + 1u) & quantified_bits_mask;
        ++quantified_bits_inc_count;

        if (quantified_bits_inc_count >= quantified_bits_inc_limit) {
            has_more_quantif_symbols = false;

            free_bits_val = ((free_bits_val | quantified_bits_mask) + 1u) & ~quantified_bits_mask;
            ++free_bits_inc_count;
            if (free_bits_inc_count >= free_bits_inc_limit) finished = true;
        }

        return free_bits_val | quantified_bits_val;
    }

    u64 init_quantif() {
        quantified_bits_val = 0u;
        quantified_bits_inc_count = 0u;
        has_more_quantif_symbols = true;

        return (free_bits_val | quantified_bits_val);
    }

    void reset() {
        free_bits_val = 0u;
        free_bits_inc_count = 0u;
        quantified_bits_val = 0u;
        quantified_bits_inc_count = 0u;
        finished = false;
    }
};


struct Formula_Pool { // Formula memory management
    unordered_set<Quantified_Atom_Conjunction> formulae;

    Formula_Pool() {
        auto top = Quantified_Atom_Conjunction(true);
        auto bottom = Quantified_Atom_Conjunction(false);
        formulae.insert(top);
        formulae.insert(bottom);
    }

    const Quantified_Atom_Conjunction* store_formula(Quantified_Atom_Conjunction& formula);
};

struct Structured_Macrostate {
    bool is_accepting = false;
    u64 handle;
    map<const Quantified_Atom_Conjunction*, list<Conjunction_State>> formulae;

    bool operator==(const Structured_Macrostate& other) const;
};

template <>
struct std::hash<Structured_Macrostate> {
    std::size_t operator() (const Structured_Macrostate& macrostate) const {
        std::size_t hash = 0u;

        // Hash the states constituting macrostate
        for (auto& [formula, states] : macrostate.formulae) {
            std::size_t formula_hash = std::hash<const Quantified_Atom_Conjunction*>{}(formula);

            // @Simplicity: Is hashing really useful in this case?
            for (auto& state: states) {
                std::size_t state_hash = 0u;

                for (auto state_constant: state.constants) {
                    std::size_t state_constant_hash = ((state_constant >> 16) ^ state_constant) * 0x45d9f3b;
                    state_constant = ((state_constant >> 16) ^ state_constant) * 0x45d9f3b;
                    state_constant = (state_constant >> 16) ^ state_constant;

                    state_hash += 0x9e3779b9 + (state_constant << 6) + (state_constant >> 2);
                }

                // @Note: hash combination must be order independent here, because we do not impose any ordering on the states
                formula_hash ^= state_hash;
            }

            hash += 0x9e3779b9 + (formula_hash << 6) + (formula_hash >> 2);
        }

        hash += macrostate.is_accepting * 33;

        return hash;
    }
};

char convert_cube_bit_to_char(u8 cube_bit);
void show_transitions_from_state(std::stringstream& output, const NFA& nfa, State origin, sylvan::MTBDD mtbdd);
NFA build_nfa_with_formula_entailement(const Quantified_Atom_Conjunction* formula, Conjunction_State& init_state, sylvan::BDDSET bdd_vars, Formula_Pool& formula_pool);
void init_mtbdd_libs();

struct Lazy_Construction_State {
    Formula_Pool& formula_pool;
    vector<Structured_Macrostate>& output_queue;
    unordered_map<Structured_Macrostate, u64> known_macrostates;
    unordered_set<u64> accepting_macrostates;

    bool is_trap_state_needed;
    State trap_state_handle;
    bool is_top_state_needed;
    State topstate_handle;
};


struct Stateful_Formula {
    Conjunction_State state;
    Quantified_Atom_Conjunction formula;
};

enum class Bound_Type : unsigned {
    NONE  = 0x00,
    UPPER = 0x01,
    LOWER = 0x02,
};

Dep_Graph build_dep_graph(const Quantified_Atom_Conjunction& conj);
void write_dep_graph_dot(std::ostream& output, Dep_Graph& graph);
void identify_potential_variables(Dep_Graph& graph);
Conjunction_State simplify_graph_using_value(Dep_Graph& graph, Conjunction_State& state, u64 var, s64 val);
void simplify_graph_with_unbound_var(Dep_Graph& graph, u64 var);
bool simplify_graph(Dep_Graph& graph, Conjunction_State& state);
Stateful_Formula convert_graph_into_formula(Dep_Graph& graph, Conjunction_State& state);
std::pair<const Quantified_Atom_Conjunction*, Conjunction_State> simplify_stateful_formula(const Quantified_Atom_Conjunction* formula, Conjunction_State& state, Formula_Pool& pool);

template <typename T>
void vector_remove(vector<T>& vec, T& elem) {
    auto pos = std::find(vec.begin(), vec.end(), elem);
    if (pos != vec.end()) {
        vec.erase(pos);
    }
}

template <typename T>
bool vector_contains(vector<T>& vec, T& elem) {
    auto pos = std::find(vec.begin(), vec.end(), elem);
    return (pos != vec.end());
}

s64 compute_multiplicative_inverse(s64 modulus, s64 a);

#endif
