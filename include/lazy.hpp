#ifndef AMAYA_LAZY_H
#define AMAYA_LAZY_H

#include "base.hpp"
#include "custom_leaf.hpp"
#include "operations.hpp"

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

struct Variable_Bound {
    bool is_present;
    u64 atom_idx;
};

struct Variable_Bounds {
    Variable_Bound lower;
    Variable_Bound upper;
    Preferred_Var_Value_Type preferred_value;
};

struct Variable_Bound_Analysis_Result {
    /* unordered_set<u64> vars_with_both_bounds; */
    bool has_var_with_both_bounds;
    vector<Variable_Bounds> bounds;
    vector<vector<u64>> congruences_per_var;
};

struct Conjuction_State;

struct Quantified_Atom_Conjunction {
    vector<Presburger_Atom> atoms;
    vector<u64>             bound_vars;
    u64 var_count;
    Variable_Bound_Analysis_Result* bounds_analysis_result;
    bool is_bottom = false;

    bool operator==(const Quantified_Atom_Conjunction& other) const;

    std::string fmt_with_state(Conjuction_State& state) const;
};

template <>
struct std::hash<Quantified_Atom_Conjunction> {
    std::size_t operator() (const Quantified_Atom_Conjunction& formula) const {
        std::size_t hash = formula.is_bottom ? 0 : 33;
        for (auto& atom: formula.atoms) {
            std::size_t atom_hash = std::hash<Presburger_Atom>{}(atom);
            hash = hash + 0x9e3779b9 + (atom_hash << 6) + (atom_hash >> 2);
        }
        return hash;
    }
};

struct Conjuction_State {
    const Quantified_Atom_Conjunction* formula;
    vector<s64> constants;

    Conjuction_State(const Quantified_Atom_Conjunction* formula, vector<s64> value): formula(formula), constants(value) {}

    Conjuction_State(const Conjuction_State& other): formula(other.formula), constants(other.constants) {}

    void post(unordered_set<Conjuction_State>& known_states, vector<Conjuction_State>& dest);
    optional<Conjuction_State> successor_along_symbol(u64 symbol);
    bool accepts_last_symbol(u64 symbol);
    bool operator==(const Conjuction_State& other) const;
};

template <>
struct std::hash<Conjuction_State> {
    std::size_t operator() (const Conjuction_State& state) const {
        std::size_t hash = std::hash<const Quantified_Atom_Conjunction*>{}(state.formula);

        for (u64 i = 0u; i < state.formula->atoms.size(); i++) {
            std::size_t atom_val_hash = std::hash<s64>{}(state.constants[i]);
            hash = hash + 0x9e3779b9 + (atom_val_hash << 6) + (atom_val_hash >> 2);
        }

        return hash;
    }
};

std::ostream& operator<<(std::ostream& output, const Presburger_Atom& atom);
std::ostream& operator<<(std::ostream& output, const Quantified_Atom_Conjunction& formula);
std::ostream& operator<<(std::ostream& output, const Conjuction_State& atom);

struct Entaiment_Status {
    bool has_no_integer_solution;
    u64 removed_atom_count;
    optional<Conjuction_State> state;
};


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
    Quantified_Atom_Conjunction top;    // Formula with solution space spanning entire space
    Quantified_Atom_Conjunction bottom; // Formula with no solution space

    unordered_set<Quantified_Atom_Conjunction> formulae;

    Formula_Pool() {
        top    = Quantified_Atom_Conjunction{.is_bottom = false};
        bottom = Quantified_Atom_Conjunction{.is_bottom = true};
    }

    const Quantified_Atom_Conjunction* store_formula(Quantified_Atom_Conjunction& formula);
};

struct Structured_Macrostate {
    bool is_accepting = false;
    u64 handle;
    map<const Quantified_Atom_Conjunction*, list<Conjuction_State>> formulae;

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

                formula_hash += 0x9e3779b9 + (state_hash << 6) + (state_hash >> 2);
            }

            hash += 0x9e3779b9 + (formula_hash << 6) + (formula_hash >> 2);
        }

        hash += macrostate.is_accepting * 33;

        return hash;
    }
};

char convert_cube_bit_to_char(u8 cube_bit);
void show_transitions_from_state(std::stringstream& output, const NFA& nfa, State origin, sylvan::MTBDD mtbdd);
NFA build_nfa_with_formula_entailement(Formula_Pool& formula_pool, Conjuction_State& init_state, sylvan::BDDSET bdd_vars);
void init_mtbdd_libs();

struct Lazy_Construction_State {
    Formula_Pool& formula_pool;
    vector<Structured_Macrostate>& output_queue;
    unordered_map<Structured_Macrostate, u64> known_macrostates;
    unordered_set<u64> accepting_macrostates;

    bool is_trap_state_needed;
    State trap_state_handle;
};

#endif
