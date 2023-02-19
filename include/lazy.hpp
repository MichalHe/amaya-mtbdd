#ifndef AMAYA_LAZY_H
#define AMAYA_LAZY_H
#include "base.hpp"

#include <list>
#include <cassert>
#include <optional>
#include <cstring>
#include <iostream>
#include <string>
#include <map>
#include <unordered_set>

typedef uint64_t u64;
typedef uint8_t  u8;
typedef int64_t  s64;

using std::optional;
using std::vector;
using std::unordered_set;

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

    s64 compute_post_along_sym(s64 constant, u64 symbol_bits) const;
    std::string fmt_with_rhs(s64 rhs) const;
    bool operator==(const Presburger_Atom& other) const;
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

struct Conjuction_State {
    const Quantified_Atom_Conjunction* formula;
    vector<s64> constants;

    Conjuction_State(const Quantified_Atom_Conjunction* formula, vector<s64> value): formula(formula), constants(value) {}

    Conjuction_State(const Conjuction_State& other): formula(other.formula), constants(other.constants) {}

    void post(unordered_set<Conjuction_State>& known_states, vector<Conjuction_State>& dest);
    Conjuction_State successor_along_symbol(u64 symbol);
    bool operator==(const Conjuction_State& other) const;
};


std::ostream& operator<<(std::ostream& output, const Presburger_Atom& atom);
std::ostream& operator<<(std::ostream& output, const Quantified_Atom_Conjunction& formula);
std::ostream& operator<<(std::ostream& output, const Conjuction_State& atom);

struct Entaiment_Status {
    bool has_no_integer_solution;
    u64 removed_atom_count;
    optional<Conjuction_State> state;
};

namespace std {
    template <>
    struct hash<Conjuction_State> {
        std::size_t operator() (const Conjuction_State& state) const {
            std::size_t hash = std::hash<const Quantified_Atom_Conjunction*>{}(state.formula);

            for (u64 i = 0u; i < state.formula->atoms.size(); i++) {
                std::size_t atom_val_hash = std::hash<s64>{}(state.constants[i]);
                hash = hash + 0x9e3779b9 + (atom_val_hash << 6) + (atom_val_hash >> 2);
            }

            return hash;
        }
    };

    template <>
    struct hash<Presburger_Atom> {
        std::size_t operator() (const Presburger_Atom& state) const {
            std::size_t hash = 0;
            for (auto coef: state.coefs) {
                std::size_t coef_hash = std::hash<s64>{}(coef);
                hash = hash + 0x9e3779b9 + (coef_hash << 6) + (coef_hash >> 2);
            }
            return hash;
        }
    };

    template <>
    struct hash<Quantified_Atom_Conjunction> {
        std::size_t operator() (const Quantified_Atom_Conjunction& formula) const {
            std::size_t hash = formula.is_bottom ? 0 : 33;
            for (auto& atom: formula.atoms) {
                std::size_t atom_hash = std::hash<Presburger_Atom>{}(atom);
                hash = hash + 0x9e3779b9 + (atom_hash << 6) + (atom_hash >> 2);
            }
            return hash;
        }
    };

    template <>
    struct hash< map<const Quantified_Atom_Conjunction*, list<Conjuction_State>> > {
        std::size_t operator() (const map<const Quantified_Atom_Conjunction*, list<Conjuction_State>>& post) const {
            std::size_t hash = 0u;

            for (auto& [formula, states] : post) {
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
            return hash;
        }
    };
}


struct Alphabet_Iterator {
    u64 free_bits_val;
    u64 free_bits_inc_count;

    u64 quantified_bits_val;
    u64 quantified_bits_inc_count;

    bool finished;

    u64 quantified_bits_mask;
    u64 quantified_bits_inc_limit;
    u64 free_bits_inc_limit;

    Alphabet_Iterator(u64 var_count, const vector<u64>& quantified_vars):
        free_bits_val(0u),
        free_bits_inc_count(0u),
        quantified_bits_val(0u),
        quantified_bits_inc_count(0u),
        finished(false)
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
            quantified_bits_inc_count = 0u;
            quantified_bits_val = 0u;

            free_bits_val = ((free_bits_val | quantified_bits_mask) + 1u) & ~quantified_bits_mask;
            ++free_bits_inc_count;

            if (free_bits_inc_count >= free_bits_inc_limit) finished = true;
        }

        return free_bits_val | quantified_bits_val;
    }

    u64 init() {
        free_bits_val = 0u;
        free_bits_inc_count = 0u;
        quantified_bits_val = 0u;
        quantified_bits_inc_count = 0u;
        finished = false;

        return 0u; // First symbol is always 0
    }
};


struct FormulaPool { // Formula memory management
    Quantified_Atom_Conjunction top;
    Quantified_Atom_Conjunction bottom;

    unordered_set<Quantified_Atom_Conjunction> formulae;

    FormulaPool() {
        top    = Quantified_Atom_Conjunction{.is_bottom = false};
        bottom = Quantified_Atom_Conjunction{.is_bottom = true};
    }

    const Quantified_Atom_Conjunction* store_formula(Quantified_Atom_Conjunction& formula);
};
#endif
