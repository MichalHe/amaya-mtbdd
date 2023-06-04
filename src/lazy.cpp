#include "../include/lazy.hpp"

#include <stdlib.h>

#include <bitset>
#include <chrono>
#include <iostream>
#include <list>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

using std::list;
using std::map;
using std::vector;
using std::unordered_set;
using std::unordered_map;

using std::pair;
using std::optional;

#define DEBUG_RUNTIME 0

typedef Quantified_Atom_Conjunction Formula;

optional<s64> Presburger_Atom::compute_post_along_sym(s64 state, u64 symbol_bits) const {
    s64 dot = 0;
    for (u64 var_i = 0; var_i < coefs.size(); var_i++) {
        s64 is_bit_set = (symbol_bits & (1u << var_i)) > 0;
        dot += is_bit_set * coefs[var_i];
    }

    if (this->type == PR_ATOM_CONGRUENCE) {
        s64 post = state - dot;
        post += modulus * ((post % 2) != 0);
        post /= 2;
        post = post % modulus;
        post += modulus * (post < 0);
        return post;
    }

    s64 post = (state - dot);

    s64 post_div_2 = post / 2;
    s64 post_mod_2 = post % 2;

    if (this->type == PR_ATOM_EQ && post_mod_2) {
        return std::nullopt;
    }

    // Floor the division
    post_div_2 -= (post_mod_2 != 0) * (post < 0);
    return post_div_2;
}

bool Presburger_Atom::accepts_last_symbol(s64 state, u64 symbol) const {
    s64 dot = 0;
    for (u64 var_i = 0; var_i < coefs.size(); var_i++) {
        s64 is_bit_set = (symbol & (1u << var_i)) > 0;
        dot += is_bit_set * coefs[var_i];
    }

    // @Optimize: This can be converted to branchless code by having a local copy of modulus and setting it to 1 if
    //            the one in atom is 0 (thus the atom is congruence). Then the result is an or between the remaining
    //            branches

    if (type == PR_ATOM_CONGRUENCE) {
        return ((state + dot) % modulus) == 0;
    } else if (type == PR_ATOM_EQ) {
        return (state + dot) == 0;
    } else {
        return (state + dot) >= 0;
    }
}


bool Conjunction_State::operator==(const Conjunction_State& other) const {
    return constants == other.constants;
}

bool Presburger_Atom::operator==(const Presburger_Atom& other) const {
    return coefs == other.coefs;
}

bool Quantified_Atom_Conjunction::operator==(const Quantified_Atom_Conjunction& other) const {
    return other.bottom == bottom && other.var_count == var_count && other.atoms == atoms && other.bound_vars == other.bound_vars;
}

optional<Conjunction_State> compute_successor(const Formula* formula, Conjunction_State& state, u64 symbol) {
    vector<s64> successor_values(state.constants.size());
    for (u64 atom_i = 0; atom_i < formula->atoms.size(); atom_i++) {
        s64 current_state = state.constants[atom_i];
        auto successor_value = formula->atoms[atom_i].compute_post_along_sym(current_state, symbol);
        if (!successor_value.has_value())
            return std::nullopt;
        successor_values[atom_i] = successor_value.value();
    }

    return Conjunction_State {.constants = successor_values};
}

bool accepts_last_symbol(const Formula* formula, Conjunction_State& state, u64 symbol) {
    bool accepts = true;
    for (u64 atom_i = 0; atom_i < formula->atoms.size(); atom_i++) {
        s64 current_state = state.constants[atom_i];
        accepts &= formula->atoms[atom_i].accepts_last_symbol(current_state, symbol);
    }
    return accepts;
}


vector<Presburger_Atom> broadcast_atoms_to_same_dimension(const vector<Sparse_Presburger_Atom>& atoms) {
    unordered_set<u64> referenced_vars;

    u64 max_var_id = 0;
    for (auto& atom: atoms) {
        for (u64 var_i = 0u; var_i < atom.var_count; var_i++) {
            max_var_id = max_var_id > atom.variables[var_i] ? max_var_id : atom.variables[var_i];
            referenced_vars.insert(atom.variables[var_i]);
        }
    }

    if (max_var_id == 0) return {};

    u64 var_count = max_var_id;
    vector<s64> dense_atom_coefs(var_count);

    vector<Presburger_Atom> dense_atoms;
    dense_atoms.reserve(atoms.size());

    for (auto& atom: atoms) {
        for (auto sparse_var_i = 0; sparse_var_i < atom.var_count; sparse_var_i++) {
            auto var = atom.variables[sparse_var_i];
            dense_atom_coefs[var-1] = atom.coefs[sparse_var_i];
        }

        Presburger_Atom dense_atom = Presburger_Atom(atom.type, dense_atom_coefs);
        dense_atoms.push_back(dense_atom);
    }

    return dense_atoms;
}

constexpr
Preferred_Var_Value_Type operator|(Preferred_Var_Value_Type left, Preferred_Var_Value_Type right)
{
    auto left_raw = static_cast<std::underlying_type_t<Preferred_Var_Value_Type>>(left);
    auto right_raw = static_cast<std::underlying_type_t<Preferred_Var_Value_Type>>(right);
    return static_cast<Preferred_Var_Value_Type>(left_raw | right_raw);
}

constexpr
Preferred_Var_Value_Type operator&(Preferred_Var_Value_Type left, Preferred_Var_Value_Type right)
{
    auto left_raw = static_cast<std::underlying_type_t<Preferred_Var_Value_Type>>(left);
    auto right_raw = static_cast<std::underlying_type_t<Preferred_Var_Value_Type>>(right);
    return static_cast<Preferred_Var_Value_Type>(left_raw & right_raw);
}

std::ostream& operator<<(std::ostream& output, const Presburger_Atom& atom) {
    bool is_first_nonzero_coef = true;
    for (u64 var_i = 0; var_i < atom.coefs.size(); var_i++) {
        if (atom.coefs[var_i]) {
            if (!is_first_nonzero_coef) output << " ";
            else is_first_nonzero_coef = false;

            output << atom.coefs[var_i] << "*x" << var_i;
        }
    }

    switch (atom.type) {
        case PR_ATOM_INEQ:
            output << " <= ?";
            break;
        case PR_ATOM_EQ:
            output << " = ?";
            break;
        case PR_ATOM_CONGRUENCE:
            output << " ~ ?";
            break;
        default:
            output << "  ???? invalid atom";
    }
    return output;
}

std::ostream& operator<<(std::ostream& output, const Quantified_Atom_Conjunction& formula) {
    output << "(";
    if (!formula.bound_vars.empty()) {
        output << "exists (";
        auto quantified_vars_iter = formula.bound_vars.begin();
        output << "x" << *quantified_vars_iter;
        ++quantified_vars_iter;
        for (; quantified_vars_iter != formula.bound_vars.end(); ++quantified_vars_iter)
            output << ", " << "x" << *quantified_vars_iter;
        output << ")";
    }

    for (auto& atom: formula.atoms) {
        output << "("<< atom <<") ";
    }

    output << ")";
    return output;
}

std::ostream& operator<<(std::ostream& output, const Conjunction_State& atom) {
    output << "(";
    if (!atom.constants.empty()) {
        auto constants_it = atom.constants.begin();
        output << *constants_it;
        ++constants_it;

        for (; constants_it != atom.constants.end(); ++constants_it) {
            output << ", " << *constants_it;
        }
    }

    output << ")"; // "[" << atom.formula << "]";
    return output;
}

std::ostream& operator<<(std::ostream& output, const Structured_Macrostate& macrostate) {
    output << "{";
    u64 written_formulae = 0;
    for (auto& [formula, states]: macrostate.formulae) {
        if (written_formulae > 0) output << ", ";
        output << "`" << *formula << "`: [";
        auto state_it = states.begin();
        output << *state_it;
        ++state_it;

        for (; state_it != states.end(); ++state_it) {
            output << ", " << *state_it;
        }
        output << "]";
        written_formulae += 1;
    }
    if (macrostate.is_accepting) {
        output << " (+ ACCEPTING)";
    }
    output << "}";
    return output;
}


bool Structured_Macrostate::operator==(const Structured_Macrostate& other) const {
    if (is_accepting != other.is_accepting) return false;
    bool is_equal = formulae == other.formulae;
    return is_equal;
};

s64 div_bound_by_coef(s64 bound, s64 coef) {
    s64 d = bound / coef;
    s64 m = bound % coef; // Needed to tell -3x <= -3 apart from -3x <= -2 so we can make corrections
    // Negative modulo means:
    // 1) lower_bound)  (coef < 0 && bound < 0): e.g. -2x <= -3  ->  x >= 3/2   == x >= 1, therefore we must make a correction or
    // 2) upper_bound)  (coef > 0 && bound < 0): e.g.  2x <= -3  ->  x <= -3/2  == x <= -1, therefore we must make a correction
    d += (m < 0);
    return d;
}

s64 normalize_reminder(s64 reminder, s64 modulus) {
    reminder %= modulus;
    reminder += (reminder < 0) * modulus;
    return reminder;
}

Dep_Graph build_dep_graph(const Quantified_Atom_Conjunction& conj) {
    Dep_Graph graph;
    graph.var_nodes = vector<Var_Node>(conj.var_count);
    graph.atom_nodes.reserve(conj.atoms.size());

    for (u64 atom_i = 0; atom_i < conj.atoms.size(); atom_i++) {
        auto& atom = conj.atoms[atom_i];
        graph.atom_nodes.push_back({.atom = atom, .atom_i = atom_i, .vars = {}});
        Atom_Node& atom_node = graph.atom_nodes[graph.atom_nodes.size() - 1];  // Get a reference to the inserted node

        for (u64 var = 0; var < conj.var_count; var++) {
            if (atom.coefs[var] != 0) {
                auto& var_node = graph.var_nodes[var];
                if (atom.type == PR_ATOM_INEQ) {
                    if (atom.coefs[var] < 0) var_node.lower_bounds.push_back(atom_i);
                    else var_node.upper_bounds.push_back(atom_i);
                } else if (atom.type == PR_ATOM_EQ) {
                    var_node.upper_bounds.push_back(atom_i);
                    var_node.lower_bounds.push_back(atom_i);
                    var_node.hard_lower_bound = {.is_present = true, .atom_i = atom_i};
                    var_node.hard_upper_bound = {.is_present = true, .atom_i = atom_i};
                } else {
                    var_node.congruences.push_back(atom_i);
                }
                atom_node.vars.push_back(var);
            }
        }
    }

    // Detect hard bounds
    for (auto& atom_node: graph.atom_nodes) {
        if (atom_node.vars.size() == 1) {
            auto var = atom_node.vars[0];
            auto coef = atom_node.atom.coefs[var];
            if (atom_node.atom.type == PR_ATOM_INEQ) {
                Hard_Bound bound_info = {.is_present = true, .atom_i = atom_node.atom_i};
                if (coef > 0) graph.var_nodes[var].hard_upper_bound = bound_info;
                else graph.var_nodes[var].hard_lower_bound = bound_info;
            }
        }
    }

    graph.quantified_vars = conj.bound_vars;
    return graph;
}


unordered_set<u64> compute_free_vars_affected_via_congruences(Dep_Graph& graph, u64 var) {
    auto& var_node = graph.var_nodes[var];
    unordered_set<u64> affected_vars;
    unordered_set<u64> known_nodes (var_node.congruences.begin(), var_node.congruences.end());
    vector<u64> work_list(var_node.congruences);

    while (!work_list.empty()) {
        auto current_atom_i = work_list.back();
        auto& current_atom = graph.atom_nodes[current_atom_i];
        work_list.pop_back();

        for (auto atom_var: current_atom.vars) {
            if  (atom_var == var) continue;

            if (std::find(graph.quantified_vars.begin(), graph.quantified_vars.end(), atom_var) == graph.quantified_vars.end()) {
                affected_vars.insert(atom_var);
            }

            auto& node = graph.var_nodes[atom_var];
            for (auto atom: node.congruences) {
                if (!known_nodes.contains(atom)) {
                    known_nodes.insert(atom);
                    work_list.push_back(atom);
                }
            }
            for (auto atom: node.lower_bounds) {
                if (!known_nodes.contains(atom)) {
                    known_nodes.insert(atom);
                    work_list.push_back(atom);
                }
            }
            for (auto atom: node.upper_bounds) {
                if (!known_nodes.contains(atom)) {
                    known_nodes.insert(atom);
                    work_list.push_back(atom);
                }
            }
        }
    }

    return affected_vars;
}

unordered_set<u64> compute_free_vars_affected_directly(Dep_Graph& graph, u64 var) {
    unordered_set<u64> affected_vars;
    auto& var_node = graph.var_nodes[var];
    for (auto atom_i: var_node.lower_bounds) {
        auto& atom = graph.atom_nodes[atom_i];
        for (auto var: atom.vars) {
            if (!vector_contains(graph.quantified_vars, var)) {
                affected_vars.insert(var);
            }
        }
    }
    for (auto atom_i: var_node.upper_bounds) {
        auto& atom = graph.atom_nodes[atom_i];
        for (auto var: atom.vars) {
            if (!vector_contains(graph.quantified_vars, var)) {
                affected_vars.insert(var);
            }
        }
    }
    return affected_vars;
}

bool is_safe_to_instantiate_var_wrt_deps(Dep_Graph& graph, u64 var) {
    // Check if picking a value and instantiating will influence multiple free vars as that can cause problems
    unordered_set<u64> affected_via_congruences = compute_free_vars_affected_via_congruences(graph, var);
    unordered_set<u64> affected_directly = compute_free_vars_affected_directly(graph, var);
    return (affected_directly.empty() || affected_via_congruences.empty());
}


bool is_safe_to_instantiate(Dep_Graph& graph, u64 var) {
    auto& node = graph.var_nodes[var];

    bool has_both_hard_bounds = node.hard_lower_bound.is_present && node.hard_upper_bound.is_present;
    if (has_both_hard_bounds) return true; // Can be simplified in the future when the bounds imply a constant value

    // The variable is missing at least one hard bound, check if its ussages all desire the same kind of value
    bool has_soft_upper_bounds = node.upper_bounds.size() > (node.hard_upper_bound.is_present);
    bool has_soft_lower_bounds = node.lower_bounds.size() > (node.hard_lower_bound.is_present);

    if (has_soft_lower_bounds && has_soft_upper_bounds) return false; // Both low and high values are desired

    // Var is unbound from at least one direction, all contexts agree on the same value kind
    bool var_can_be_neg_inf = node.lower_bounds.empty() && !node.hard_lower_bound.is_present;
    bool var_can_be_pos_inf = node.upper_bounds.empty() && !node.hard_upper_bound.is_present;
    if (var_can_be_neg_inf || var_can_be_pos_inf) {
        // Problem is when both are non-empty simulatneously
        return is_safe_to_instantiate_var_wrt_deps(graph, var);
    }

    // Variable has only a hard bound - its context requires opposite value, and thus, we can instantiate it
    // - pick value near bound
    if (node.hard_lower_bound.is_present && !has_soft_lower_bounds) return true;
    if (node.hard_upper_bound.is_present && !has_soft_upper_bounds) return true;

    return false;
}

void identify_potential_variables(Dep_Graph& graph) {
    graph.potential_vars.clear();
    for (auto var: graph.quantified_vars) {
        auto& node = graph.var_nodes[var];
        if (is_safe_to_instantiate(graph, var)) graph.potential_vars.push_back(var);
    }
}

void write_dep_graph_dot(std::ostream& output, Dep_Graph& graph) {
    output << "graph deps {\n";

    for (u64 var = 0; var < graph.var_nodes.size(); var++) {
        bool is_quantified = (std::find(graph.quantified_vars.begin(), graph.quantified_vars.end(), var) != graph.quantified_vars.end());
        std::string quantif_attrs = is_quantified  ? ", color=green" : "";
        output << "  v" << var << " [label=\"x" << var << "\",shape=box" <<  quantif_attrs << "]\n";
    }

    for (u64 atom_i = 0; atom_i < graph.atom_nodes.size(); atom_i++) {
        auto& node = graph.atom_nodes[atom_i];

        u64 is_hard_bound = 0;
        for (auto var: node.vars) {
            auto& var_node = graph.var_nodes[var];
            if (var_node.is_hard_lower_bound(atom_i)) is_hard_bound = 1;
            if (var_node.is_hard_upper_bound(atom_i)) is_hard_bound = 2;
        }

        std::string hard_bound_attrs;
        switch (is_hard_bound) {
            case 1: // Lower
                hard_bound_attrs = ", style=filled, fillcolor=\"#6bb7ea\"";
                break;
            case 2: // Upper
                hard_bound_attrs = ", style=filled, fillcolor=\"#ea6b78\"";
                break;
            default:
                hard_bound_attrs = "";
                break;
        }

        output << "  a" << atom_i << " [label=\"" << node.atom << "\"" << hard_bound_attrs << "]\n";
    }

    for (u64 var = 0; var < graph.var_nodes.size(); var++) {
        auto& var_node = graph.var_nodes[var];
        for (auto& atom_i: var_node.congruences) {
            output << "  v" << var << " -- " << "a" << atom_i << "\n";
        }

        for (auto atom_i: var_node.lower_bounds) {
            output << "  v" << var << " -- " << "a" << atom_i << " [color=\"blue\"]" << "\n";
        }
        for (auto atom_i: var_node.upper_bounds) {
            output << "  v" << var << " -- " << "a" << atom_i << " [color=\"red\"]" << "\n";
        }
    }

    for (auto quantif_var: graph.quantified_vars) {
        auto affected_vars = compute_free_vars_affected_via_congruences(graph, quantif_var);
        for (auto affected_var: affected_vars) {
            output << "  v" << quantif_var << " -- v" << affected_var << " [style=dotted]\n";
        }
    }

    output << "}\n";
}

bool get_constant_value_implied_by_bounds(Dep_Graph& graph, u64 var, Conjunction_State& state, s64* val) {
    auto& var_node = graph.var_nodes[var];
    if (!var_node.hard_lower_bound.is_present || !var_node.hard_upper_bound.is_present) {
        return false;
    }

    auto& upper_bound_i = var_node.hard_upper_bound.atom_i;
    auto& upper_bound_node = graph.atom_nodes[upper_bound_i];
    s64 upper_bound_val = div_bound_by_coef(state.constants[upper_bound_i],
                                            upper_bound_node.atom.coefs[var]);

    auto& lower_bound_i = var_node.hard_lower_bound.atom_i;
    auto& lower_bound_node = graph.atom_nodes[lower_bound_i];
    s64 lower_bound_val = div_bound_by_coef(state.constants[lower_bound_i],
                                            upper_bound_node.atom.coefs[var]);

    if (upper_bound_val == lower_bound_val) {
        *val = upper_bound_val;
        return true;
    }
    return false;
}

Conjunction_State simplify_graph_using_value(Dep_Graph& graph, Conjunction_State& state, u64 var, s64 val) {
    auto& var_node = graph.var_nodes[var];

    for (auto atom_i: var_node.lower_bounds) {
        auto& atom_node = graph.atom_nodes[atom_i];
        state.constants[atom_i] -= atom_node.atom.coefs[var] * val;
        atom_node.atom.coefs[var] = 0;
        vector_remove(atom_node.vars, var);
        if (atom_node.vars.empty()) atom_node.is_satisfied = true;
    }
    for (auto atom_i: var_node.upper_bounds) {
        auto& atom_node = graph.atom_nodes[atom_i];
        state.constants[atom_i] -= atom_node.atom.coefs[var] * val;
        atom_node.atom.coefs[var] = 0;
        vector_remove(atom_node.vars, var);
        if (atom_node.vars.empty()) atom_node.is_satisfied = true;
    }
    for (auto atom_i: var_node.congruences) {
        auto& atom_node = graph.atom_nodes[atom_i];
        state.constants[atom_i] -= atom_node.atom.coefs[var] * val;
        state.constants[atom_i] = normalize_reminder(state.constants[atom_i], atom_node.atom.modulus);
        atom_node.atom.coefs[var] = 0;
        vector_remove(atom_node.vars, var);
        if (atom_node.vars.empty()) atom_node.is_satisfied = true;
    }

    return state;
}

void mark_lin_atom_as_satisfied(Dep_Graph& graph, u64 atom_i, u64 unbound_var) {
    auto& atom_node = graph.atom_nodes[atom_i];
    atom_node.is_satisfied = true;
    for (auto atom_var: atom_node.vars) {
        if (atom_var == unbound_var) continue;  // The node for var causing atom to be satisfied will be "removed", do not modify it
        auto& dep_var_node = graph.var_nodes[atom_var];
        vector_remove(dep_var_node.lower_bounds, atom_i);
        vector_remove(dep_var_node.upper_bounds, atom_i);
    }
    atom_node.vars.clear();
}

void simplify_graph_with_unbound_var(Dep_Graph& graph, u64 var) {
    auto& var_node = graph.var_nodes[var];

    for (auto atom_i: var_node.lower_bounds) {
        mark_lin_atom_as_satisfied(graph, atom_i, var);
    }

    for (auto atom_i: var_node.upper_bounds) {
        mark_lin_atom_as_satisfied(graph, atom_i, var);
    }

    for (auto atom_i: var_node.congruences) {
        auto& atom_node = graph.atom_nodes[atom_i];
        atom_node.is_satisfied = true;
        for (auto atom_var: atom_node.vars) {
            if (atom_var == var) continue;
            auto& dep_var_node = graph.var_nodes[atom_var];
            vector_remove(dep_var_node.congruences, atom_i);
        }
        atom_node.vars.clear();
    }

    var_node.congruences.clear();
    var_node.lower_bounds.clear();
    var_node.upper_bounds.clear();
}

s64 compute_multiplicative_inverse(s64 modulus, s64 a) {
    s64 prev_reminder = modulus;
    s64 reminder = a;
    s64 prev_s = 1;
    s64 s = 0;
    s64 prev_t = 0;
    s64 t = 1;

    s64 aux;
    while (reminder != 0) {
        s64 q = prev_reminder / reminder;

        aux = reminder;
        reminder = prev_reminder - q * reminder;
        prev_reminder = aux;

        aux = s;
        s = prev_s - q * s;
        prev_s = aux;

        aux = t;
        t = prev_t - q * t;
        prev_t = aux;
    }

    // Normalize prev_t to be in [0..modulus)
    prev_t += (prev_t < 0) * modulus;
    return prev_t;
}


bool get_value_close_to_bounds(Dep_Graph& graph, Conjunction_State& state, u64 var, s64* value) {
    if (!is_safe_to_instantiate(graph, var)) {
        return false;
    }

    auto& var_node = graph.var_nodes[var];
    u64 bound_node_i = 0;
    Bound_Type bound_type = Bound_Type::NONE;
    bool has_only_hard_upper_bound = var_node.upper_bounds.size() <= (var_node.hard_upper_bound.is_present);
    bool has_only_hard_lower_bound = var_node.lower_bounds.size() <= (var_node.hard_lower_bound.is_present);

    if (var_node.hard_lower_bound.is_present && has_only_hard_lower_bound) {
        // All reminding bounds are upper, preferring the smallest value possible
        bound_node_i = var_node.hard_lower_bound.atom_i;
        bound_type = Bound_Type::LOWER;
    } else if (var_node.hard_upper_bound.is_present && has_only_hard_upper_bound) {
        // All reminding bounds are lower, preferring the largest value possible
        bound_node_i = var_node.hard_upper_bound.atom_i;
        bound_type = Bound_Type::UPPER;
    }

    if (bound_type == Bound_Type::NONE) return false;

    auto& bound_node = graph.atom_nodes[bound_node_i];

    // Var has a clear lower bound and its all of its ussages would like it to be some low value
    s64 bound_value = div_bound_by_coef(state.constants[bound_node_i], bound_node.atom.coefs[var]);
    if (var_node.congruences.empty()) {
        *value = bound_value;
        return true;
    }

    if (var_node.congruences.size() > 1) return false;

    auto congruence_i = var_node.congruences[0];
    auto& congruence_node = graph.atom_nodes[congruence_i];

    if (congruence_node.vars.size() > 1) {
        // @Research: In case there are multiple variables in a congruence, but we affect only free vars,
        // we can  maybe still instantiate the value if the congruence variables are not restristed too much.
        return false;   // We don't know how to infer value when there are multiple vars in congruence
    }

    auto& congruence = graph.atom_nodes[congruence_i].atom;

    u64 nonzero_coefs = 0;
    for (auto coef: congruence.coefs) nonzero_coefs += (coef != 0);
    if (nonzero_coefs > 1) return false;

    // Find a solution in the unshifted congruence range, e.g., -y ~ 100 (mod 303)
    s64 modulus = congruence.modulus;

    s64 rhs = normalize_reminder(state.constants[congruence_i], modulus);
    s64 coef = normalize_reminder(congruence.coefs[var], modulus);

    s64 mult_inv = compute_multiplicative_inverse(modulus, coef);
    rhs = (rhs * mult_inv) % modulus; // e.g., y ~ 200 (mod 303), rhs is now a solution

    s64 shift_coef = bound_value / congruence.modulus;
    shift_coef -= (bound_value < 0); // A signed floor division
    s64 congruence_shift = congruence.modulus * shift_coef;

    s64 instantiated_value = rhs + congruence_shift;

    // Both bound and the proposed solution lie in the same block of reminders [kM ... (k+1)M]
    // but the proposed solution might be slightly above/bellow the bound, correct it.
    if (bound_type == Bound_Type::LOWER) {
        // (x >= 5) and the suggested solution is x = 3
        instantiated_value += (instantiated_value < bound_value) * modulus;
    } else {
        // (x <= 5) and the suggested solution is x = 8
        instantiated_value -= (bound_value < instantiated_value) * modulus;
    }

    *value = instantiated_value;
    return true;
}


bool simplify_graph(Dep_Graph& graph, Conjunction_State& state) {
    bool all_vars_probed = false;
    bool was_graph_simplified = false;
    while (!all_vars_probed) {
        bool was_graph_invalidated = false;
        for (auto potential_var: graph.potential_vars) {
            s64 inst_val = 0;
            bool is_val_implied = get_constant_value_implied_by_bounds(graph, potential_var, state, &inst_val);
            if (is_val_implied) {
                PRINT_DEBUG("Simplifying graph - bounds imply value " << inst_val << " for x" << potential_var);
                state = simplify_graph_using_value(graph, state, potential_var, inst_val);
                vector_remove(graph.quantified_vars, potential_var);
                was_graph_invalidated = true;
                was_graph_simplified = true;
                break;
            }

            // Try simplifying unbound vars
            auto& var_node = graph.var_nodes[potential_var];
            bool can_be_neg_inf = var_node.lower_bounds.empty();
            bool can_be_pos_inf = var_node.upper_bounds.empty();
            if (can_be_neg_inf || can_be_pos_inf) {
                PRINT_DEBUG("Simplifying graph on unbound var: x" << potential_var);
                simplify_graph_with_unbound_var(graph, potential_var);
                vector_remove(graph.quantified_vars, potential_var);
                was_graph_simplified = true;
                was_graph_invalidated = true;
                break;
            }

            // Try simplifying using a var with a hard bound such that the var context directly
            // tells what such a value should be
            bool can_be_instantiated = get_value_close_to_bounds(graph, state, potential_var, &inst_val);
            if (can_be_instantiated) {
                PRINT_DEBUG("Simplifying graph - instantiating var: x" << potential_var << " with value " << inst_val);
                state = simplify_graph_using_value(graph, state, potential_var, inst_val);
                vector_remove(graph.quantified_vars, potential_var);
                was_graph_simplified = true;
                was_graph_invalidated = true;
                break;
            }
        }
        if (was_graph_invalidated) identify_potential_variables(graph);
        else all_vars_probed = true;
    }
    return was_graph_simplified;
}

Stateful_Formula convert_graph_into_formula(Dep_Graph& graph, Conjunction_State& state) {
    // @Note: In theory, we could return a conjunction of the same number of atoms with satisfied atoms replaced
    //        by TRUE and reuse the graph troughout the entire execution, however, this would require successor
    //        computation to be aware about this and ignore such atoms.
    vector<Presburger_Atom> atoms;
    vector<s64> formula_state;
    for (u64 atom_i = 0; atom_i < graph.atom_nodes.size(); atom_i++) {
        auto& atom_node = graph.atom_nodes[atom_i];
        if (!atom_node.is_satisfied) {
            formula_state.push_back(state.constants[atom_i]);
            atoms.push_back(atom_node.atom);
        }
    }

    Quantified_Atom_Conjunction conjunction(atoms, graph.quantified_vars, graph.var_nodes.size()); // new dep graph is automatically created
    Conjunction_State new_state(formula_state);
    return {.state = new_state, .formula = conjunction};
}

std::pair<const Formula*, Conjunction_State> simplify_stateful_formula(const Formula* formula, Conjunction_State& state, Formula_Pool& pool) {
    Dep_Graph graph_copy(formula->dep_graph);
    bool was_graph_simplified = simplify_graph(graph_copy, state);  // If the graph is simplified, then state will be modified

    if (was_graph_simplified) {
        auto new_stateful_formula = convert_graph_into_formula(graph_copy, state);
        auto formula_ptr = pool.store_formula(new_stateful_formula.formula);
        return {formula_ptr, new_stateful_formula.state};
    }
    return {formula, state};
}


void add_var_coef_term(std::stringstream& dest, u64 var, s64 coef) {
    switch (coef) {
        case 1:
            dest << "x" << var;
            break;
        case -1:
            dest << "(- x" << var << ")";
            break;
        default:
            dest << "(* " << coef << " x" << var << ")";
    }
}

std::string Presburger_Atom::fmt_with_rhs(s64 rhs) const {
    std::stringstream str_builder;

    str_builder << "(";
    if (type == Presburger_Atom_Type::PR_ATOM_INEQ) str_builder << "<= ";
    else if (type == Presburger_Atom_Type::PR_ATOM_CONGRUENCE || Presburger_Atom_Type::PR_ATOM_EQ) str_builder << "= ";
    else assert(0 && "The value of atom->type should be set!");

    if (type == Presburger_Atom_Type::PR_ATOM_CONGRUENCE)
        str_builder << "(mod ";

    u64 nonzero_coefficients = 0;
    bool is_first_coef_write_queued = false;
    u64 first_coef_i = 0;
    for (u64 var_i = 0u; var_i < coefs.size(); var_i++) {
        if (!coefs[var_i]) continue;

        // Postpone the writing of the first coefficient so that we know whether or not to print "(+"
        if (!nonzero_coefficients) {
            nonzero_coefficients++;
            is_first_coef_write_queued = true;
            first_coef_i = var_i;
            continue;
        }

        nonzero_coefficients++;

        if (is_first_coef_write_queued) {
            str_builder << "(+ ";
            str_builder << "x" << first_coef_i;
            is_first_coef_write_queued = false;
        }

        str_builder << " ";
        add_var_coef_term(str_builder, var_i, coefs[var_i]);
    }

    if (is_first_coef_write_queued)
        add_var_coef_term(str_builder, first_coef_i, coefs[first_coef_i]);

    if (nonzero_coefficients > 1)
        str_builder << ")";

    if (type == Presburger_Atom_Type::PR_ATOM_CONGRUENCE)
        str_builder << " " << modulus << ")";

    str_builder << " " << rhs << ")";

    return str_builder.str();
}

std::string Quantified_Atom_Conjunction::fmt_with_state(Conjunction_State& state) const {
    if (this->is_bottom())
        return "false";
    else if (this->is_top()) {
        return "true";
    }

    std::stringstream str_builder;

    str_builder << "(exists (";
    if (!bound_vars.empty()) {
        auto quantified_vars_iter = bound_vars.begin();
        str_builder << "(x" << *quantified_vars_iter << " Int)";
        ++quantified_vars_iter;
        for (; quantified_vars_iter != bound_vars.end(); ++quantified_vars_iter) {
            str_builder << " (x" << *quantified_vars_iter << " Int)";
        }
    }
    str_builder << ") ";

    if (atoms.size() == 1) {
        str_builder << atoms[0].fmt_with_rhs(state.constants[0]);
    } else {
        str_builder << "(land ";
        str_builder << atoms[0].fmt_with_rhs(state.constants[0]);
        for (u64 atom_i = 1u; atom_i < atoms.size(); atom_i++) {
            str_builder << " " << atoms[atom_i].fmt_with_rhs(state.constants[atom_i]);
        }
        str_builder << ")";
    }

    str_builder << ")";  // Closing the (exists
    return str_builder.str();
}


void insert_successor_into_post_if_valuable(Structured_Macrostate& post, const Formula* formula, Conjunction_State& successor) {
    auto& bucket = post.formulae[formula];

    bool should_be_inserted = true;
    bool insert_position_found = false;

    /*
     * @Optimize: The bucket is kept in a sorted order, so once we arrive at a state that is lexicographically
     *            larger, we don't have to iterate further.
     */

    auto bucket_iter = bucket.begin();
    list<Conjunction_State>::iterator insert_position;
    while (bucket_iter != bucket.end()) {
        auto& other_successor = *bucket_iter;

        // Compute pareto optimality
        bool is_smaller = false; // successor <= other_successor at some fragment of the state (constant)
        bool is_larger  = false; // successor >= other_successor

        bool is_other_smaller_then_inserted = true;

        for (u64 state_fragment_i = 0u; state_fragment_i < successor.constants.size(); state_fragment_i++) {
            auto& atom = formula->atoms[state_fragment_i];

            bool is_simulation_plain = (atom.type == PR_ATOM_EQ || atom.type == PR_ATOM_CONGRUENCE);

            if (is_simulation_plain) {
                if (successor.constants[state_fragment_i] != other_successor.constants[state_fragment_i]) {
                    // Successor and other successor are incomparable
                    is_smaller = true;
                    is_larger = true;
                    break;
                }
            }
            else {
                // @Optimize: Those IFs can be avoided by bitpacking the use booleans into one machine word and using bit OPs
                if (successor.constants[state_fragment_i] < other_successor.constants[state_fragment_i]) {
                    is_smaller = true; // The successor is simulated at the current field
                }
                else if (successor.constants[state_fragment_i] > other_successor.constants[state_fragment_i]) {
                    is_larger = true;
                }
            }

            if (other_successor.constants[state_fragment_i] > successor.constants[state_fragment_i]) {
                is_other_smaller_then_inserted = false;
            }
        }

        // The current successor should be inserted after all <= states were exhaused (insertion sort)
        if (!is_other_smaller_then_inserted && !insert_position_found) {
            // @Note: The usage of the insert_position iterator should not get invalidated future list erasures, because
            //        the iterator points to a state that is lexicographically larger than the inserted state. As the state
            //        is larger, it cannot be that it is simulated by the inserted - either it is larger in some plain constant,
            //        or larger in some inequation -> no simulation.
            insert_position = bucket_iter;
            insert_position_found = true;
        }

        if (is_larger && is_smaller) { // Incomparable
            ++bucket_iter;
            continue;
        }
        else {
            if (is_larger) {
                bucket.erase(bucket_iter++);
            } else {
                // Either the inserted atom is smaller (is_smaller = true), or they are the same (is_smaller = is_larger = false)
                should_be_inserted = false;
                break;
            }
        }

    }

    if (should_be_inserted) {
        if (!insert_position_found) insert_position = bucket.end();
        bucket.insert(insert_position, successor);
    }
}

void make_macrostate_canoical(Structured_Macrostate& macrostate) {
    auto comparator = [](const Conjunction_State& left, const Conjunction_State& right){
        // Is left < right ?
        for (u64 i = 0; i < left.constants.size(); i++) {
            if (left.constants[i] > right.constants[i]) return false;
            else if (left.constants[i] < right.constants[i]) return true;
        }
        return false;
    };

    for (auto& [formula, atoms]: macrostate.formulae) {
        atoms.sort(comparator);
    }
}

Structured_Macrostate make_trivial_macrostate(Formula_Pool& pool, bool polarity) {
    Formula formula(polarity);
    auto formula_ptr = pool.store_formula(formula);
    Conjunction_State state ({});
    list<Conjunction_State> state_list = {state};
    map<const Formula*, list<Conjunction_State>> contents = {{formula_ptr, state_list}};

    Structured_Macrostate macrostate = {.is_accepting = polarity, .formulae = contents};
    return macrostate;
}


void explore_macrostate(NFA& constructed_nfa,
                        Structured_Macrostate& macrostate,
                        Alphabet_Iterator& alphabet_iter,
                        Lazy_Construction_State& constr_state)
{
    alphabet_iter.reset();

    PRINT_DEBUG("Exploring " << macrostate);

    while (!alphabet_iter.finished) {
        Structured_Macrostate post;
        bool is_accepting = false;

        u64 transition_symbol = alphabet_iter.init_quantif();  // The quantified bits will be masked away, so it is sufficient to take the first one

        bool is_post_top = false;
        for (u64 symbol = transition_symbol; alphabet_iter.has_more_quantif_symbols; symbol = alphabet_iter.next_symbol()) {
#if DEBUG_RUNTIME
            auto start = std::chrono::system_clock::now();
#endif
            for (auto& [formula, states]: macrostate.formulae) {
                for (auto& state: states) {
                    auto maybe_successor = compute_successor(formula, state, symbol);

                    if (!maybe_successor.has_value()) continue;

                    auto successor = maybe_successor.value();

                    auto new_formula_state_pair = simplify_stateful_formula(formula, successor, constr_state.formula_pool);
                    successor = new_formula_state_pair.second;
                    auto simplified_formula = new_formula_state_pair.first;

                    if (simplified_formula->is_top()) {
                        is_post_top = true;
                        break;
                    } else if (!simplified_formula->is_bottom()) {
                        insert_successor_into_post_if_valuable(post, simplified_formula, successor);
                    }

                    post.is_accepting |= accepts_last_symbol(simplified_formula, successor, symbol);
                }

                if (is_post_top) break;
            }

            make_macrostate_canoical(post);

#if DEBUG_RUNTIME
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            std::cout << "Number of formulae: " << macrostate.formulae.size()
                      << "; insertion took: " << elapsed_seconds.count() << " seconds." << std::endl;
#endif
        }

        if (post.formulae.empty()) {
            constr_state.is_trap_state_needed = true;
            constructed_nfa.add_transition(macrostate.handle, constr_state.trap_state_handle, transition_symbol, alphabet_iter.quantified_bits_mask);
            continue;
        };

        if (is_post_top) {
            constr_state.is_trap_state_needed = true;
            constructed_nfa.add_transition(macrostate.handle, constr_state.trap_state_handle, transition_symbol, alphabet_iter.quantified_bits_mask);
            continue;
        }

        // Assign a unique integer to every state
        post.handle = constr_state.known_macrostates.size();
        auto [element_iter, was_inserted] = constr_state.known_macrostates.emplace(post, post.handle);
        if (was_inserted) {
            // @Simplicity: Maybe the known_macrostates should be a set instead of a map since we are storing the handle
            //              inside the macrostate either way.
            constr_state.output_queue.push_back(post);

            if (post.is_accepting) { // Do the hash-query only if we see the macrostate for the first time
                constr_state.accepting_macrostates.emplace(post.handle);
            }
        } else {
            post.handle = element_iter->second;  // Use the already existing handle
        }

        constructed_nfa.add_transition(macrostate.handle, post.handle, transition_symbol, alphabet_iter.quantified_bits_mask);
    }
}

// @Cleanup: Figure out where should we put this
void init_mtbdd_libs() {
    int n_workers = 1;
    size_t dequeue_size = 10000000;

    lace_init(n_workers, dequeue_size);
    const size_t stack_size = 1LL << 20;
    lace_startup(0, NULL, NULL);
    sylvan::sylvan_set_limits(500LL*1024*1024, 3, 5);
    sylvan::sylvan_init_package();
    sylvan::sylvan_init_mtbdd();

    mtbdd_leaf_type_set = sylvan::sylvan_mt_create_type();
    sylvan::sylvan_mt_set_hash(mtbdd_leaf_type_set, set_leaf_hash);
    sylvan::sylvan_mt_set_equals(mtbdd_leaf_type_set, set_leaf_equals);
    sylvan::sylvan_mt_set_create(mtbdd_leaf_type_set, mk_set_leaf);
    sylvan::sylvan_mt_set_destroy(mtbdd_leaf_type_set, destroy_set_leaf);
    sylvan::sylvan_mt_set_to_str(mtbdd_leaf_type_set, set_leaf_to_str);
}


// @Cleanup: Move this into base.cpp

char convert_cube_bit_to_char(u8 cube_bit) {
    switch (cube_bit) {
        case 0:
            return '0';
        case 1:
            return '1';
        default:
            return '*';
    }
}


void show_transitions_from_state(std::stringstream& output, const NFA& nfa, State origin, sylvan::MTBDD mtbdd) {
    u8 symbol[nfa.var_count];

    sylvan::MTBDD leaf = sylvan::mtbdd_enum_first(mtbdd, nfa.vars, symbol, NULL);
    while (leaf != sylvan::mtbdd_false) {
        output << origin << " (" << convert_cube_bit_to_char(symbol[0]);
        for (u64 symbol_i = 1; symbol_i < nfa.var_count; symbol_i++) {
            output << ", " << convert_cube_bit_to_char(symbol[symbol_i]);
        }
        output << ") ";

        auto leaf_contents = reinterpret_cast<Transition_Destination_Set*>(sylvan::mtbdd_getvalue(leaf));

        if (leaf_contents->destination_set.empty()) {
            output << "{}";
            continue;
        }

        auto contents_iter = leaf_contents->destination_set.begin();
        output << "{" << *contents_iter;
        ++contents_iter;
        for (; contents_iter != leaf_contents->destination_set.end(); ++contents_iter) {
            output << ", " << *contents_iter;
        }
        output << "}";

 	    leaf = sylvan::mtbdd_enum_next(mtbdd, nfa.vars, symbol, NULL);
        if (leaf != sylvan::mtbdd_false) output << "\n";
    }
}

std::string NFA::show_transitions() const {
    std::stringstream str_builder;

    if (transitions.empty()) {
        return "{}\n";
    }

    std::cout << "Number of states: " << transitions.size() << std::endl;
    for (auto& [origin_state, mtbdd]: this->transitions) {
        show_transitions_from_state(str_builder, *this, origin_state, mtbdd);
        str_builder << std::endl;
    }

    return str_builder.str();
}

NFA build_nfa_with_formula_entailement(const Formula* formula, Conjunction_State& init_state, sylvan::BDDSET bdd_vars, Formula_Pool& formula_pool) {
    vector<Structured_Macrostate> work_queue;
    Lazy_Construction_State constr_state = {.formula_pool = formula_pool, .output_queue = work_queue};

    // Prepare a trapstate in case it is needed
    {
        auto trapstate = make_trivial_macrostate(formula_pool, false);
        auto [container_pos, was_insterted] = constr_state.known_macrostates.emplace(trapstate, constr_state.known_macrostates.size());
        constr_state.trap_state_handle = container_pos->second;
    }

    // Prepare TRUE state
    {
        auto topstate = make_trivial_macrostate(formula_pool, true);
        auto [container_pos, was_insterted] = constr_state.known_macrostates.emplace(topstate, constr_state.known_macrostates.size());
        constr_state.topstate_handle = container_pos->second;
    }

    PRINT_DEBUG("Performing initial simplification of formula");
    auto initial_simplification = simplify_stateful_formula(formula, init_state, formula_pool);
    init_state = initial_simplification.second;
    formula = initial_simplification.first;

    if (formula->is_top()) {
        State state = 0;
        NFA nfa(bdd_vars, formula->var_count, {state}, {state}, {state});
        nfa.add_transition(state, state, 0u, static_cast<u64>(-1));
        return nfa;
    }

    // Populate the queue with the initial states
    const u64 init_state_handle = constr_state.known_macrostates.size();
    {
        list<Conjunction_State> init_list = {init_state};
        map<const Formula*, list<Conjunction_State>> init_macrostate_formulae = { {formula, init_list} };
        Structured_Macrostate init_macrostate = { .is_accepting = false, .handle = init_state_handle, .formulae = init_macrostate_formulae};
        work_queue.push_back(init_macrostate);

        auto [container_pos, was_inserted] = constr_state.known_macrostates.emplace(init_macrostate, constr_state.known_macrostates.size());
    }


    NFA nfa(bdd_vars, formula->var_count, {}, {}, {static_cast<State>(init_state_handle)});

    Alphabet_Iterator alphabet_iter = Alphabet_Iterator(formula->var_count, formula->bound_vars);
    PRINT_DEBUG("The following vars are quantified: " << formula->bound_vars);
    PRINT_DEBUG("Using the following mask: " << std::bitset<8>(alphabet_iter.quantified_bits_mask));

    u64 processed = 0;
    while (!work_queue.empty()) {
        auto macrostate = work_queue.back();
        work_queue.pop_back();
        // PRINT_DEBUG("Proccessing" << macrostate << " (Queue size " << work_queue.size() << ')');

        nfa.states.insert(macrostate.handle);
        if (macrostate.is_accepting) {
            nfa.final_states.insert(macrostate.handle);
        }

        explore_macrostate(nfa, macrostate, alphabet_iter, constr_state);
        processed += 1;
    }

    const u64 all_bits_dont_care_mask = static_cast<u64>(-1);
    if (constr_state.is_trap_state_needed) {
        nfa.states.insert(constr_state.trap_state_handle);
        nfa.add_transition(constr_state.trap_state_handle, constr_state.trap_state_handle, 0u, all_bits_dont_care_mask);
    }

    if (constr_state.is_top_state_needed) {
        nfa.states.insert(constr_state.topstate_handle);
        nfa.add_transition(constr_state.topstate_handle, constr_state.topstate_handle, 0u, all_bits_dont_care_mask);
    }


    PRINTF_DEBUG("The constructed NFA has %lu states", nfa.states.size());

    for (auto& [macrostate, handle]: constr_state.known_macrostates) {
        std::cout << handle << " :: " << macrostate << std::endl;
    }

    nfa.perform_pad_closure();
    return determinize_nfa(nfa);
}


const Quantified_Atom_Conjunction* Formula_Pool::store_formula(Quantified_Atom_Conjunction& formula) {
    auto [it, did_insertion_happen] = formulae.emplace(formula);
    const Quantified_Atom_Conjunction* stored_formula_ptr = &(*it);
    return &(*it);
}

Quantified_Atom_Conjunction::Quantified_Atom_Conjunction(const vector<Presburger_Atom>& atoms, const vector<u64>& bound_vars, u64 var_count) {
    this->atoms = atoms;
    this->bound_vars = bound_vars;
    this->var_count = var_count;
    this->dep_graph = build_dep_graph(*this);
    identify_potential_variables(this->dep_graph);
}
