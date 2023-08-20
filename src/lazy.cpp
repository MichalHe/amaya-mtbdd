#include "../include/lazy.hpp"

#include <cmath>
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
#include <utility>

using std::list;
using std::map;
using std::vector;
using std::unordered_set;

using std::pair;
using std::optional;

#define DEBUG_RUNTIME 0

typedef Quantified_Atom_Conjunction Formula;

std::size_t hash_array(Sized_Array<s64>& arr) {
    std::size_t hash = 0;
    for (u64 i = 0; i < arr.size; i++) {
        std::size_t coef_hash = std::hash<s64>{}(arr.items[i]);
        hash = hash + 0x9e3779b9 + (coef_hash << 6) + (coef_hash >> 2);
    }
    return hash;
}

s64 vector_dot(const Sized_Array<s64>& coefs, const u64 symbol) {
    s64 dot = 0;
    for (u64 var_i = 0; var_i < coefs.size; var_i++) {
        s64 is_bit_set = (symbol & (1u << var_i)) > 0;
        dot += is_bit_set * coefs.items[var_i];
    }
    return dot;
};

s64 combine_moduli(s64 mod_2pow, s64 mod_odd) {
    s64 modulus = (1 << mod_2pow) * mod_odd;
    return modulus;
}

typedef s64 Modulus;

optional<s64> congruence_compute_post(Congruence& congruence, s64 state, u64 symbol) {
    s64 dot  = vector_dot(congruence.coefs, symbol);
    s64 post = state - dot;

    if (congruence.modulus_2pow > 0) {
        if (post % 2 == 0) {
            post /= 2;
            s64 new_modulus = combine_moduli(congruence.modulus_2pow - 1, congruence.modulus_odd);
            post %= new_modulus;
            post += (post < 0) * new_modulus;
            return post;
        }
        return std::nullopt;
    }

    post += congruence.modulus_odd * ((post % 2) != 0);
    post /= 2;
    post = post % congruence.modulus_odd;
    post += congruence.modulus_odd * (post < 0);
    return post;
};

bool congruence_accepts_symbol(Congruence& congruence, s64 state, u64 symbol) {
    s64 dot = vector_dot(congruence.coefs, symbol);
    s64 modulus = combine_moduli(congruence.modulus_2pow, congruence.modulus_odd);
    return ((state + dot) % modulus) == 0;
}

std::optional<s64> equation_compute_post(Equation& eq, s64 state, u64 symbol) {
    s64 dot = vector_dot(eq.coefs, symbol);
    s64 post = state - dot;
    s64 post_div_2 = post / 2;
    s64 post_mod_2 = post % 2;

    if (post_mod_2) return std::nullopt;

    post_div_2 -= (post_mod_2 != 0) * (post < 0); // Floor the division
    return post_div_2;
}

bool equation_accepts_symbol(Equation& eq, s64 state, u64 symbol) {
    s64 dot = vector_dot(eq.coefs, symbol);
    return (state + dot) == 0;
}

s64 inequation_compute_post(Inequation& ineq, s64 state, u64 symbol) {
    s64 dot = vector_dot(ineq.coefs, symbol);
    s64 post = state - dot;
    s64 post_div_2 = post / 2;
    s64 post_mod_2 = post % 2;
    post_div_2 -= (post_mod_2 != 0) * (post < 0); // Floor the division
    return post_div_2;
}

bool inequation_accepts_symbol(Inequation& ineq, s64 state, u64 symbol) {
    s64 dot = vector_dot(ineq.coefs, symbol);
    return (state + dot) >= 0;
}

bool Conjunction_State::operator==(const Conjunction_State& other) const {
    return constants == other.constants;
}

bool Quantified_Atom_Conjunction::operator==(const Quantified_Atom_Conjunction& other) const {
    return other.bottom == bottom &&
           other.var_count == var_count &&
           other.congruences == congruences &&
           other.inequations == inequations &&
           other.equations == equations &&
           other.bound_vars == other.bound_vars;
}

struct Successor {
    const Formula* formula = nullptr;
    Conjunction_State state;
};

optional<Successor> compute_successor(Formula_Pool& formula_pool, const Formula* formula, s64* state_data, u64 symbol) {
    vector<s64> successor_values(formula->atom_count());

    u64 atom_idx = 0;
    for (auto& cong: formula->congruences) {
        auto maybe_successor = congruence_compute_post(cong, state_data[atom_idx], symbol);

        if (!maybe_successor.has_value()) {
            return std::nullopt;
        }

        s64 new_state = maybe_successor.value();

        successor_values[atom_idx] = new_state;
        atom_idx += 1;
    }

    for (auto& eq: formula->equations) {
        auto successor_value = equation_compute_post(eq, state_data[atom_idx], symbol);

        if (!successor_value.has_value()) return std::nullopt;

        successor_values[atom_idx] = successor_value.value();
        atom_idx += 1;
    }

    for (auto& ineq: formula->inequations) {
        s64 post = inequation_compute_post(ineq, state_data[atom_idx], symbol);
        successor_values[atom_idx] = post;
        atom_idx += 1;
    };

    Conjunction_State post(successor_values);

    // Determine if we actually will need to allocate memory
    if (!formula->post_formula) {
        // All components, including all congruences had post (none returned empty optional). As
        // FA are synchronous, it does not matter what symbol did we compute post with since with
        // every other symbol with successful post computation we will modify the congruence moduli
        // in the same way although the target state might be different. The was_post_computed
        // field is not hashed, so it is OK to drop const and set it.
        bool should_allocate = false;
        for (auto& cong: formula->congruences) {
            if (cong.modulus_2pow > 0) {
                should_allocate = true;
                break;
            }
        }

        if (should_allocate) {
            Sized_Array<Congruence> new_congruences = formula_pool.allocator.alloc_congruences(formula->congruences.size);

            for (u64 congruence_i = 0; congruence_i < formula->congruences.size; congruence_i++) {
                Congruence& old_congruence = formula->congruences.items[congruence_i];
                Congruence new_congruence = {
                    .coefs = old_congruence.coefs,
                    .modulus_odd = old_congruence.modulus_odd,
                    .modulus_2pow = old_congruence.modulus_2pow - 1
                };
                new_congruences.items[congruence_i] = new_congruence;
            }

            Formula new_formula = Formula(new_congruences, formula->equations, formula->inequations,
                                          formula->bound_vars, formula->var_count);
            const Formula* new_formula_canon_ptr = formula_pool.store_formula(new_formula);

            { // Make sure we do not need to recompute the thing
                const_cast<Formula*>(formula)->post_formula = new_formula_canon_ptr;
            }

            Successor successor = {.formula = new_formula_canon_ptr, .state = post};
            return successor;
        } else {
            // All congruences have odd moduli, make sure we don't try allocating a new formula due
            // to modulus change ever again (self loop in dependency graph)
            const_cast<Formula*>(formula)->post_formula = formula;
        }
    }

    Successor successor = {.formula = formula->post_formula, .state = post};
    return successor;
}

bool accepts_last_symbol(const Formula* formula, s64* state_data, u64 symbol) {
    bool accepts = true;

    u64 atom_i = 0;
    for (auto& eq: formula->equations) {
        accepts &= equation_accepts_symbol(eq, state_data[atom_i], symbol);
        atom_i += 1;
    }

    for (auto& congr: formula->congruences) {
        accepts &= congruence_accepts_symbol(congr, state_data[atom_i], symbol);
        atom_i += 1;
    }

    for (auto& ineq: formula->inequations) {
        accepts &= inequation_accepts_symbol(ineq, state_data[atom_i], symbol);
        atom_i += 1;
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
            output << " ~ ? (mod " << atom.modulus << ")";
            break;
        default:
            output << "  ???? invalid atom";
    }
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


void write_memview(std::ostream& output, s64* data, s64 count) {
    for (s64 i = 0; i < count; i++) {
        if (i > 0) output << ", ";
        output << data[i];
    }
}

std::ostream& operator<<(std::ostream& output, const Finalized_Macrostate& macrostate) {
    output << "{";
    s64 constants_written = 0;
    for (s64 formula_idx = 0; formula_idx < macrostate.header.size; formula_idx++) {
        const Macrostate_Header_Elem& header = macrostate.header.items[formula_idx];

        if (formula_idx > 0) output << ", ";
        output << "`" << *header.formula << "`: [";

        for (s64 state_idx = 0; state_idx < header.state_cnt; state_idx++) {
            if (state_idx > 0) output << ", ";

            output << "(";
            write_memview(output, macrostate.state_data.items + constants_written, header.formula->atom_count());
            output << ")";
            constants_written += header.formula->atom_count();
        }

        output << "]";
    }

    if (macrostate.is_accepting) {
        output << " (+ ACCEPTING)";
    }
    output << "}";
    return output;
}

bool Finalized_Macrostate::operator==(const Finalized_Macrostate& other) const {
    if (is_accepting != other.is_accepting) return false;
    if (header.size  != other.header.size) return false;
    if (state_data.size != other.state_data.size) return false;

    // @Optimize: Ensure that these loops get properly vectorized
    for (s64 i = 0; i < header.size; i++) {
        bool is_equal = (header.items[i].formula == other.header.items[i].formula) && (header.items[i].state_cnt == other.header.items[i].state_cnt);
        if (!is_equal) return false;
    }

    for (s64 i = 0; i < state_data.size; i++) {
        if (state_data.items[i] != other.state_data.items[i]) return false;
    }

    return true;
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
    graph.linear_nodes.reserve(conj.equations.size + conj.inequations.size);
    graph.congruence_nodes.reserve(conj.congruences.size);

    u64 lin_atom_i = 0;
    for (auto& eq: conj.equations) {
        vector<s64> coefs(eq.coefs.begin(), eq.coefs.end());
        Linear_Node atom_node = {.type = Linear_Node_Type::EQ, .coefs = coefs, .vars = {}, .is_satisfied = false};

        for (auto var = 0; var < conj.var_count; var++) {
            Var_Node& var_node = graph.var_nodes[var];
            s64 coef           = eq.coefs.items[var];
            if (coef) {
                var_node.upper_bounds.push_back(lin_atom_i);
                var_node.lower_bounds.push_back(lin_atom_i);
                atom_node.vars.push_back(var);
            }
        }

        graph.linear_nodes.push_back(atom_node);
        lin_atom_i += 1;
    }

    for (auto& ineq: conj.inequations) {
        vector<s64> coefs(ineq.coefs.begin(), ineq.coefs.end());
        Linear_Node atom_node = {.type = Linear_Node_Type::INEQ, .coefs = coefs, .vars = {}, .is_satisfied = false};

        for (auto var = 0; var < conj.var_count; var++) {
            Var_Node& var_node = graph.var_nodes[var];
            s64 coef           = ineq.coefs.items[var];
            if (coef) {
                if (coef < 0) var_node.lower_bounds.push_back(lin_atom_i);
                else          var_node.upper_bounds.push_back(lin_atom_i);
                atom_node.vars.push_back(var);
            }
        }

        graph.linear_nodes.push_back(atom_node);
        lin_atom_i += 1;
    }

    u64 congruence_i = 0;
    for (auto& congr: conj.congruences) {
        vector<s64> coefs(congr.coefs.begin(), congr.coefs.end());
        Congruence_Node congruence_node = {
                .coefs = coefs,
                .vars = {},
                .modulus_2pow = congr.modulus_2pow,
                .modulus_odd = congr.modulus_odd,
                .is_satisfied = false,
        };

        for (auto var = 0; var < conj.var_count; var++) {
            Var_Node& var_node = graph.var_nodes[var];
            s64 coef           = congr.coefs.items[var];
            if (coef) {
                var_node.congruences.push_back(congruence_i);
                congruence_node.vars.push_back(var);
            }
        }

        graph.congruence_nodes.push_back(congruence_node);
        congruence_i += 1;
    }

    // Detect hard bounds
    for (u64 lin_atom_i = 0; lin_atom_i < graph.linear_nodes.size(); lin_atom_i++) {
        Linear_Node node = graph.linear_nodes[lin_atom_i];
        if (node.vars.size() == 1) {
            u64 var            = node.vars[0];
            s64 coef           = node.coefs[var];
            Var_Node& var_node = graph.var_nodes[var];

            if (node.type == Linear_Node_Type::INEQ) {
                Hard_Bound bound_info = {.is_present = true, .atom_i = lin_atom_i};
                if (coef > 0) var_node.hard_upper_bound = bound_info;
                else          var_node.hard_lower_bound = bound_info;
            } else { // Single variable equation
                var_node.hard_upper_bound.is_present = true;
                var_node.hard_upper_bound.atom_i = lin_atom_i;

                var_node.hard_lower_bound.is_present = true;
                var_node.hard_lower_bound.atom_i = lin_atom_i;
            }
        }
    }

    graph.quantified_vars = conj.bound_vars;
    return graph;
}

void add_new_affected_vars_to_worklist(vector<u64>& potential_vars, vector<u64>& worklist, unordered_set<u64>& known_affected_vars) {
    for (u64 var: potential_vars) {
        if (!known_affected_vars.contains(var)) {
            known_affected_vars.insert(var);
            worklist.push_back(var);
        }
    }
}

unordered_set<u64> compute_overall_affected_vars(Dep_Graph& graph, u64 var) {
    typedef vector<u64> Atom_Vars;

    Var_Node& var_node = graph.var_nodes[var];
    unordered_set<u64> affected_free_vars;
    unordered_set<u64> known_affected_vars;

    // We care only about affected vars and therefore there is no need to distinguish between
    // different atom node types and we can work with atom vars directly
    vector<u64> work_list = {var};

    while (!work_list.empty()) {
        u64 affected_var = work_list.back();
        auto& var_node   = graph.var_nodes[affected_var];
        work_list.pop_back();

        if (!vector_contains(graph.quantified_vars, affected_var)) {
            affected_free_vars.insert(affected_var);
        }

        for (u64 lin_node_i: var_node.lower_bounds) {
            add_new_affected_vars_to_worklist(graph.linear_nodes[lin_node_i].vars, work_list, known_affected_vars);
        }

        for (u64 lin_node_i: var_node.upper_bounds) {
            add_new_affected_vars_to_worklist(graph.linear_nodes[lin_node_i].vars, work_list, known_affected_vars);
        }

        for (u64 congruence_i: var_node.congruences) {
            add_new_affected_vars_to_worklist(graph.congruence_nodes[congruence_i].vars, work_list, known_affected_vars);
        }
    }

    return affected_free_vars;
}

unordered_set<u64> compute_free_vars_affected_linearly(Dep_Graph& graph, u64 var) {
    unordered_set<u64> affected_vars;
    auto& var_node = graph.var_nodes[var];
    for (auto atom_i: var_node.lower_bounds) {
        Linear_Node& atom = graph.linear_nodes[atom_i];
        for (auto var: atom.vars) {
            if (!vector_contains(graph.quantified_vars, var)) {
                affected_vars.insert(var);
            }
        }
    }
    for (auto atom_i: var_node.upper_bounds) {
        Linear_Node& atom = graph.linear_nodes[atom_i];
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
    unordered_set<u64> all_affected_vars = compute_overall_affected_vars(graph, var);
    unordered_set<u64> vars_affected_linearly = compute_free_vars_affected_linearly(graph, var);

    u64 nonlinearly_affected_vars = all_affected_vars.size() - vars_affected_linearly.size();
    // We either affect unlimited number of variables directly (we can tell what value is preferred), or at most
    // one variable in nonlinearly, but then we cannot have any linearly dependent variables because a clear
    // preference of value cannot be determined.
    return (!nonlinearly_affected_vars || (nonlinearly_affected_vars == 1 && vars_affected_linearly.empty()));
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

std::string get_node_attrs_for_hard_bound(Bound_Type type) {
    switch (type) {
        case Bound_Type::LOWER:
            return ", style=filled, fillcolor=\"#6bb7ea\"";
        case Bound_Type::UPPER: // Upper
            return ", style=filled, fillcolor=\"#ea6b78\"";
        default:
            return "";
    }
}


void write_dep_graph_dot(std::ostream& output, Dep_Graph& graph) {
    output << "graph deps {\n";

    for (u64 var = 0; var < graph.var_nodes.size(); var++) {
        bool is_quantified = (std::find(graph.quantified_vars.begin(), graph.quantified_vars.end(), var) != graph.quantified_vars.end());
        std::string quantif_attrs = is_quantified  ? ", color=green" : "";
        output << "  v" << var << " [label=\"x" << var << "\",shape=box" <<  quantif_attrs << "]\n";
    }

    u64 atom_i = 0;
    for (auto& linear_node: graph.linear_nodes) {
        Bound_Type bound_type = Bound_Type::NONE;
        for (auto var: linear_node.vars) {
            auto& var_node = graph.var_nodes[var];
            if (var_node.is_hard_lower_bound(atom_i)) {
                bound_type = Bound_Type::LOWER;
                break;
            }
            if (var_node.is_hard_upper_bound(atom_i)) {
                bound_type = Bound_Type::UPPER;
                break;
            }
        }

        std::string hard_bound_attrs = get_node_attrs_for_hard_bound(bound_type);

        output << "  a" << atom_i << " [label=\"" << linear_node << "\"" << hard_bound_attrs << "]\n";
        atom_i += 1;
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
        unordered_set<u64> affected_vars          = compute_overall_affected_vars(graph, quantif_var);
        unordered_set<u64> linearly_affected_vars = compute_free_vars_affected_linearly(graph, quantif_var);
        for (auto affected_var: affected_vars) {
            if (!linearly_affected_vars.contains(affected_var)) {
                output << "  v" << quantif_var << " -- v" << affected_var << " [style=dotted]\n";
            }
        }
    }

    output << "}\n";
}

bool get_constant_value_implied_by_bounds(Dep_Graph& graph, u64 var, Conjunction_State& state, s64* val) {
    auto& var_node = graph.var_nodes[var];
    if (!var_node.hard_lower_bound.is_present || !var_node.hard_upper_bound.is_present) {
        return false;
    }

    u64 linear_states_offset = graph.congruence_nodes.size();

    u64 upper_bound_i             = var_node.hard_upper_bound.atom_i;
    Linear_Node& upper_bound_node = graph.linear_nodes[upper_bound_i];
    s64 upper_bound_val           = div_bound_by_coef(state.constants[upper_bound_i + linear_states_offset],
                                                      upper_bound_node.coefs[var]);

    u64 lower_bound_i             = var_node.hard_lower_bound.atom_i;
    Linear_Node& lower_bound_node = graph.linear_nodes[lower_bound_i];
    s64 lower_bound_val           = div_bound_by_coef(state.constants[lower_bound_i + linear_states_offset],
                                                      upper_bound_node.coefs[var]);

    if (upper_bound_val == lower_bound_val) {
        *val = upper_bound_val;
        return true;
    }
    return false;
}

void set_var_value_in_lin_atom(Dep_Graph& graph, Conjunction_State& state, u64 lin_atom_i, u64 var, s64 val) {
    auto& atom_node = graph.linear_nodes[lin_atom_i];

    u64 atom_i = lin_atom_i + graph.congruence_nodes.size();
    state.constants[atom_i] -= atom_node.coefs[var] * val;

    atom_node.coefs[var] = 0;
    vector_remove(atom_node.vars, var);
    if (atom_node.vars.empty()) atom_node.is_satisfied = true;
}

Conjunction_State simplify_graph_using_value(Dep_Graph& graph, Conjunction_State& state, u64 var, s64 val) {
    auto& var_node = graph.var_nodes[var];

    for (auto atom_i: var_node.lower_bounds)
        set_var_value_in_lin_atom(graph, state, atom_i, var, val);

    for (auto atom_i: var_node.upper_bounds)
        set_var_value_in_lin_atom(graph, state, atom_i, var, val);

    for (auto atom_i: var_node.congruences) {
        auto& congruence_node = graph.congruence_nodes[atom_i];

        // Congruences are at the front of the state, no need to introduce offset
        s64 modulus = combine_moduli(congruence_node.modulus_2pow, congruence_node.modulus_odd);
        state.constants[atom_i] -= congruence_node.coefs[var] * val;
        state.constants[atom_i] = normalize_reminder(state.constants[atom_i], modulus);

        congruence_node.coefs[var] = 0;
        vector_remove(congruence_node.vars, var);
        if (congruence_node.vars.empty()) congruence_node.is_satisfied = true;
    }

    return state;
}

void mark_lin_atom_as_satisfied(Dep_Graph& graph, u64 atom_i, u64 unbound_var) {
    auto& atom_node = graph.linear_nodes[atom_i];
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
        Congruence_Node& congruence_node = graph.congruence_nodes[atom_i];
        congruence_node.is_satisfied = true;
        for (auto atom_var: congruence_node.vars) {
            if (atom_var == var) continue;
            auto& dep_var_node = graph.var_nodes[atom_var];
            vector_remove(dep_var_node.congruences, atom_i);
        }
        congruence_node.vars.clear();
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

    auto& bound_node = graph.linear_nodes[bound_node_i];

    // Var has a clear lower bound and its all of its ussages would like it to be some low value
    s64 bound_value = div_bound_by_coef(state.constants[bound_node_i], bound_node.coefs[var]);
    if (var_node.congruences.empty()) {
        *value = bound_value;
        return true;
    }

    if (var_node.congruences.size() > 1) return false;

    u64 congruence_i                 = var_node.congruences[0];
    Congruence_Node& congruence_node = graph.congruence_nodes[congruence_i];

    if (congruence_node.vars.size() > 1) {
        // @Research: In case there are multiple variables in a congruence, but we affect only free vars,
        // we can  maybe still instantiate the value if the congruence variables are not restristed too much.
        return false;   // We don't know how to infer value when there are multiple vars in congruence
    }

    u64 nonzero_coefs = 0;
    for (auto coef: congruence_node.coefs) nonzero_coefs += (coef != 0);
    if (nonzero_coefs > 1) return false;

    // Find a solution in the unshifted congruence range, e.g., -y ~ 100 (mod 303)
    s64 modulus = combine_moduli(congruence_node.modulus_2pow, congruence_node.modulus_odd);

    s64 rhs  = normalize_reminder(state.constants[congruence_i], modulus);
    s64 coef = normalize_reminder(congruence_node.coefs[var], modulus);

    s64 mult_inv = compute_multiplicative_inverse(modulus, coef);
    rhs = (rhs * mult_inv) % modulus; // e.g., y ~ 200 (mod 303), rhs is now a solution

    s64 shift_coef = bound_value / modulus;
    shift_coef    -= (bound_value < 0); // A signed floor division
    s64 congruence_shift = modulus * shift_coef;

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

    // std::cout << "Instantiating x" << var << " with value " << instantiated_value << std::endl;
    // std::cout << congruence << std::endl;
    // std::cout << state.constants[congruence_i] << std::endl;
    // std::cout << state << std::endl;
    // std::cout << "Modulus " << modulus << std::endl;

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

Stateful_Formula convert_graph_into_formula(Dep_Graph& graph, Formula_Allocator& allocator, Conjunction_State& state) {
    vector<s64> formula_state;

    for (u64 congr_i = 0; congr_i < graph.congruence_nodes.size(); congr_i++) {
        Congruence_Node& congruence_node = graph.congruence_nodes[congr_i];
        if (congruence_node.is_satisfied) continue;

        formula_state.push_back(state.constants[congr_i]);

        Congruence* congruence = allocator.alloc_temporary_congruence();

        memcpy(congruence->coefs.items, congruence_node.coefs.data(), sizeof(s64) * congruence_node.coefs.size());

        congruence->modulus_2pow = congruence_node.modulus_2pow;
        congruence->modulus_odd = congruence_node.modulus_odd;
    }

    u64 linear_atoms_offset = graph.congruence_nodes.size();
    for (u64 lin_node_i = 0; lin_node_i < graph.linear_nodes.size(); lin_node_i++) {
        Linear_Node& atom_node = graph.linear_nodes[lin_node_i];
        if (atom_node.is_satisfied) continue;

        u64 state_coef_i = lin_node_i + linear_atoms_offset;
        formula_state.push_back(state.constants[state_coef_i]);

        if (atom_node.type == Linear_Node_Type::EQ) {
            Equation* equation = allocator.alloc_temporary_equation();
            memcpy(equation->coefs.items, atom_node.coefs.data(), sizeof(s64) * atom_node.coefs.size());
            equation->coefs.size = atom_node.coefs.size();
        } else {
            Inequation* inequation = allocator.alloc_temporary_inequation();
            memcpy(inequation->coefs.items, atom_node.coefs.data(), sizeof(s64) * atom_node.coefs.size());
            inequation->coefs.size = atom_node.coefs.size();
        }
    }

    Quantified_Atom_Conjunction conjunction = allocator.get_tmp_formula();
    conjunction.bound_vars = graph.quantified_vars;
    conjunction.var_count  = graph.var_nodes.size();

    Conjunction_State new_state(formula_state);
    return {.state = new_state, .formula = conjunction};
}

std::pair<const Formula*, Conjunction_State> simplify_stateful_formula(const Formula* formula, Conjunction_State& state, Formula_Pool& pool) {
    if (formula->dep_graph.potential_vars.empty()) return {formula, state};

    Dep_Graph graph_copy(formula->dep_graph);
    bool was_graph_simplified = simplify_graph(graph_copy, state);  // If the graph is simplified, then state will be modified

    if (was_graph_simplified) {
        auto new_stateful_formula = convert_graph_into_formula(graph_copy, pool.allocator, state);
        const auto& [formula_ptr, was_formula_new] = pool.store_formula_with_info(new_stateful_formula.formula);
        if (was_formula_new) {
            // The formula that was inserted has pointers referring to the temporary storage
            // that will be used in the future so that the formula coefficients might get
            // overwritten with something else. We have to replace the formula pointers
            // with the ones pointing to a commited memory.
            Formula commited_formula = pool.allocator.commit_tmp_space();
            // We are only changing pointers to point to a different memory with the exact same contents.
            // Since we are hashing the pointer contents and not pointers, temporary dropping const should be ok.
            Formula* modif_formula_ptr = const_cast<Formula*>(formula_ptr);
            modif_formula_ptr->congruences.items = commited_formula.congruences.items;
            modif_formula_ptr->equations.items   = commited_formula.equations.items;
            modif_formula_ptr->inequations.items = commited_formula.inequations.items;
        } else {
           pool.allocator.drop_tmp();
        }
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

void insert_pareto_optimal(std::list<Conjunction_State>& queue, const u64 ineq_offset, Conjunction_State& state) {
    bool should_be_inserted    = true;
    bool insert_position_found = false;

    /*
     * @Optimize: The bucket is kept in a sorted order, so once we arrive at a state that is lexicographically
     *            larger, we don't have to iterate further.
     */

    auto bucket_iter = queue.begin();
    list<Conjunction_State>::iterator insert_position;
    while (bucket_iter != queue.end()) {
        auto& other_state = *bucket_iter;

        // Compute pareto optimality
        bool is_smaller = false; // successor <= other_successor at some fragment of the state
        bool is_larger  = false; // successor >= other_successor at some fragment of the state
        bool is_other_smaller_than_inserted = true;

        for (u64 state_fragment_i = 0u; state_fragment_i < state.constants.size(); state_fragment_i++) {
            if (state.constants[state_fragment_i] < other_state.constants[state_fragment_i]) {
                is_smaller = true; // The successor is subsumed at the current field
                is_other_smaller_than_inserted = false;
            }
            else if (state.constants[state_fragment_i] > other_state.constants[state_fragment_i]) {
                is_larger = true;
            }
        }

        // The current successor should be inserted after all <= states were exhaused (insertion sort)
        if (!is_other_smaller_than_inserted && !insert_position_found) {
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
                queue.erase(bucket_iter++);
            } else {
                // Either the inserted atom is smaller (is_smaller = true), or they are the same (is_smaller = is_larger = false)
                should_be_inserted = false;
                break;
            }
        }

    }

    if (should_be_inserted) {
        if (!insert_position_found) insert_position = queue.end();
        queue.insert(insert_position, state);
    }
}

void insert_into_post_if_valuable2(Intermediate_Macrostate& post, const Formula* formula, Conjunction_State& successor) {
    Prefix_Table& formula_bucket = post.formulae[formula];

    u64 prefix_size = formula->congruences.size + formula->equations.size;
    vector<s64> prefix(prefix_size);
    for (u64 i = 0; i < prefix_size; i++) {
        prefix[i] = successor.constants[i];
    }

    // @Optimize: Currently we are storing the entire state (including congruence and equation data)
    //            which is redundant since we already have this information in post.

    list<Conjunction_State>& prefix_bucket = formula_bucket[prefix];
    insert_pareto_optimal(prefix_bucket, prefix_size, successor);
}

Finalized_Macrostate finalize_macrostate(Formula_Allocator& alloc, Intermediate_Macrostate& intermediate_macrostate) {
    // Finalize the intermediate macrostate into a final form that is canonical (atoms are sorted)
    auto comparator = [](const Conjunction_State& left, const Conjunction_State& right) {
        // Is left < right ?
        for (u64 i = 0; i < left.constants.size(); i++) {
            if (left.constants[i] > right.constants[i]) return false;
            else if (left.constants[i] < right.constants[i]) return true;
        }
        return false;
    };

    // formulae is a map, and therefore, the formula should  be kept in sorted order,
    // no need to sort them right now.
    Sized_Array<Macrostate_Header_Elem> headers = alloc.alloc_macrostate_headers(intermediate_macrostate.formulae.size());
    s64 header_idx = 0;
    s64 max_states_per_formula = 0;
    s64 required_constant_count = 0;
    for (const auto& [formula, prefix_table]: intermediate_macrostate.formulae) {
        headers.items[header_idx].formula   = formula;
        headers.items[header_idx].state_cnt = prefix_table.size();

        max_states_per_formula = prefix_table.size() > max_states_per_formula ? prefix_table.size() : max_states_per_formula;

        header_idx += 1;
    }

    vector<Conjunction_State> sort_space; // This can be reused by an allocator that caches the space pointer and does a malloc only if the current size < request
    sort_space.reserve(max_states_per_formula);

    // Count how many s64's will have to be in the allocated block
    for (auto& [formula, prefix_table]: intermediate_macrostate.formulae) {
        u64 states_with_this_formula = 0; // How many states are there in the prefix table
        for (auto& [prefix, prefix_entries]: prefix_table) {
            for (auto& state: prefix_entries) {
                states_with_this_formula += prefix_entries.size();
            }
        }
        required_constant_count += states_with_this_formula * formula->atom_count();
    }

    Sized_Array<s64> macrostate_data = alloc.alloc_macrostate_data(required_constant_count);
    s64 macrostate_data_next_free_slot = 0;

    for (auto& [formula, prefix_table]: intermediate_macrostate.formulae) {
        sort_space.clear();

        for (auto& [prefix, prefix_entries]: prefix_table) {
            for (auto& state: prefix_entries) {
                sort_space.push_back(state);
            }
        }
        std::sort(sort_space.begin(), sort_space.end(), comparator);

        for (const auto& state: sort_space) {
            // Memcpy
            for (const auto& constant: state.constants) {
                macrostate_data.items[macrostate_data_next_free_slot] = constant;
                macrostate_data_next_free_slot += 1;
            }
        }
    }

    return {.header = headers, .state_data = macrostate_data, .is_accepting = intermediate_macrostate.is_accepting};
}

Finalized_Macrostate make_trivial_macrostate(Formula_Pool& pool, bool polarity) {
    Formula formula(polarity);
    auto formula_ptr = pool.store_formula(formula);

    Finalized_Macrostate macrostate = {
        .header = {.items = nullptr, .size = 0},
        .state_data = {.items = nullptr, .size = 0},
        .is_accepting = polarity,
    };
    return macrostate;
}


void explore_macrostate(NFA& constructed_nfa,
                        Finalized_Macrostate& macrostate,
                        Alphabet_Iterator& alphabet_iter,
                        Lazy_Construction_State& constr_state)
{
    alphabet_iter.reset();

    // PRINT_DEBUG("Exploring " << macrostate);

    while (!alphabet_iter.finished) {
        Intermediate_Macrostate post;
        bool is_accepting = false;

        u64 transition_symbol = alphabet_iter.init_quantif();  // The quantified bits will be masked away, so it is sufficient to take the first one

        bool is_post_top = false;
        for (u64 symbol = transition_symbol; alphabet_iter.has_more_quantif_symbols; symbol = alphabet_iter.next_symbol()) {
#if DEBUG_RUNTIME
            auto start = std::chrono::system_clock::now();
#endif
            s64 state_data_read_idx = 0;
            for (auto& macrostate_header: macrostate.header) {
                const Formula* formula = macrostate_header.formula;
                for (s64 formula_state_idx = 0; formula_state_idx < macrostate_header.state_cnt; formula_state_idx++) {
                    s64* current_origin_state_data = macrostate.state_data.items + state_data_read_idx;
                    state_data_read_idx += formula->atom_count();

                    auto maybe_successor = compute_successor(constr_state.formula_pool, formula, current_origin_state_data, symbol);

                    if (!maybe_successor.has_value()) continue;

                    auto successor = maybe_successor.value();

                    // The old formula needs to be used to determine whether the post is accepting
                    post.is_accepting |= accepts_last_symbol(formula, current_origin_state_data, symbol);
                    formula = successor.formula;

                    auto new_formula_state_pair = simplify_stateful_formula(formula, successor.state, constr_state.formula_pool);
                    successor.state = new_formula_state_pair.second;
                    auto simplified_formula = new_formula_state_pair.first;

                    if (simplified_formula->is_top()) {
                        is_post_top = true;
                        break;
                    } else if (!simplified_formula->is_bottom()) {
                        insert_into_post_if_valuable2(post, simplified_formula, successor.state);
                    }
                }

                if (is_post_top) break;
            }


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

        Finalized_Macrostate finalized_post = finalize_macrostate(constr_state.formula_pool.allocator, post);

        // Assign a unique integer to every state
        finalized_post.handle = constr_state.known_macrostates.size();
        auto [element_iter, was_inserted] = constr_state.known_macrostates.emplace(finalized_post, finalized_post.handle);
        if (was_inserted) {
            // @Simplicity: Maybe the known_macrostates should be a set instead of a map since we are storing the handle
            //              inside the macrostate either way.
            constr_state.output_queue.push_back(finalized_post);

            if (post.is_accepting) { // Do the hash-query only if we see the macrostate for the first time
                constr_state.accepting_macrostates.emplace(finalized_post.handle);
            }
        } else {
            finalized_post.handle = element_iter->second;  // Use the already existing handle
        }

        constructed_nfa.add_transition(macrostate.handle, finalized_post.handle, transition_symbol, alphabet_iter.quantified_bits_mask);
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
    vector<Finalized_Macrostate> work_queue;
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
        State init_state = 0;
        State acc_state  = 1;
        NFA nfa(bdd_vars, formula->var_count, {init_state, acc_state}, {acc_state}, {init_state});
        nfa.add_transition(init_state, acc_state, 0u, static_cast<u64>(-1));
        nfa.add_transition(acc_state, acc_state, 0u, static_cast<u64>(-1));
        return nfa;
    }

    // Populate the queue with the initial states
    const u64 init_state_handle = constr_state.known_macrostates.size();
    {
        Sized_Array<Macrostate_Header_Elem> header = formula_pool.allocator.alloc_macrostate_headers(1);
        header.items[0].formula = formula;
        header.items[0].state_cnt = 1;

        Sized_Array<s64> macrostate_data = formula_pool.allocator.alloc_macrostate_data(formula->atom_count());
        for (s64 state_constant_idx = 0; state_constant_idx < init_state.constants.size(); state_constant_idx++) {
            macrostate_data.items[state_constant_idx]  = init_state.constants[state_constant_idx];
        }

        Finalized_Macrostate init_macrostate = {
            .header = header,
            .state_data = macrostate_data,
            .is_accepting = false,
            .handle = init_state_handle,
        };
        work_queue.push_back(init_macrostate);

        auto [container_pos, was_inserted] = constr_state.known_macrostates.emplace(init_macrostate, constr_state.known_macrostates.size());
    }


    NFA nfa(bdd_vars, formula->var_count, {}, {}, {static_cast<State>(init_state_handle)});

    Alphabet_Iterator alphabet_iter = Alphabet_Iterator(formula->var_count, formula->bound_vars);
    PRINT_DEBUG("The following vars are quantified: " << formula->bound_vars);
    PRINT_DEBUG("Using the following mask: " << std::bitset<8>(alphabet_iter.quantified_bits_mask));

    u64 processed = 0;
    u64 processed_batch = 0;
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
        processed_batch += 1;
        if (processed_batch >= 2000) {
            PRINTF_DEBUG("Processed %lu states\n", processed);
            processed_batch = 0;
        }
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

    PRINTF_DEBUG("The constructed NFA has %lu states\n", nfa.states.size());
    PRINT_DEBUG(nfa.states);
    for (auto& [macrostate, handle]: constr_state.known_macrostates) {
        PRINT_DEBUG(handle << " :: " << macrostate);
    }

    if (!formula->bound_vars.empty()) {
        nfa.perform_pad_closure();
        return determinize_nfa(nfa);
    }
    return nfa;
}


const Quantified_Atom_Conjunction* Formula_Pool::store_formula(Quantified_Atom_Conjunction& formula) {
    auto [it, did_insertion_happen] = formulae.emplace(formula);
    const Quantified_Atom_Conjunction* stored_formula_ptr = &(*it);
    return &(*it);
}

pair<const Formula*, bool> Formula_Pool::store_formula_with_info(Quantified_Atom_Conjunction& formula) {
    auto [it, did_insertion_happen] = formulae.emplace(formula);
    const Quantified_Atom_Conjunction* stored_formula_ptr = &(*it);
    return {&(*it), did_insertion_happen};
}

template<typename InputIt>
std::string fmt_coef_var_sum(InputIt coef_first, InputIt coef_end) {
    std::stringstream builder;

    u64 i = 0;
    for (; coef_first != coef_end; ++coef_first) {
        s64 coef  = *coef_first;
        char sign = coef >= 0 ? '+' : '-';
        if (i > 0) {
            builder << ' ' << sign << ' ';
            coef = abs(coef);
        }
        builder << coef << "*x" << i;
        i += 1;
    }

    return builder.str();
}

std::ostream& operator<<(std::ostream& output, const Linear_Node& lin_node) {
    output << fmt_coef_var_sum(lin_node.coefs.begin(), lin_node.coefs.end());
    const char* log_symbol = lin_node.type == Linear_Node_Type::EQ ? " = ?" : " <= ?";
    output << log_symbol;
    return output;
}

std::ostream& operator<<(std::ostream& output, const Congruence_Node& congruence_node) {
    s64 modulus = combine_moduli(congruence_node.modulus_2pow, congruence_node.modulus_odd);
    output << fmt_coef_var_sum(congruence_node.coefs.begin(), congruence_node.coefs.end());
    output << " ~ ? (mod " << modulus << ")";
    return output;
}

std::ostream& operator<<(std::ostream& output, const map<const Quantified_Atom_Conjunction*, Prefix_Table>& macrostate) {
    if (macrostate.empty()) {
        output << "{}";
        return output;
    }

    u64 row_cnt = 0;
    output << "{\n";
    for (auto& [formula_ptr, prefix_table]: macrostate) {
        output << "  " << formula_ptr << " :: ";

        for (auto& [prefix, states]: prefix_table) {
            for (auto& state: states) output << state << ", ";
        }

        if (row_cnt < macrostate.size()) output << "\n";
        row_cnt += 1;
    }
    std::cout << "}";
    return output;

}
std::ostream& operator<<(std::ostream& output, const Congruence& congruence) {
    output << fmt_coef_var_sum(congruence.coefs.begin(), congruence.coefs.end())
           << " ~= (mod "
           << combine_moduli(congruence.modulus_2pow, congruence.modulus_odd) << ")";
    return output;
}

std::ostream& operator<<(std::ostream& output, const Quantified_Atom_Conjunction& formula) {
    if (!formula.has_atoms()) {
        const char* formula_str = formula.is_top() ? "(TRUE)" : "(FALSE)";
        output << formula_str;
        return output;
    }

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

    for (auto& atom: formula.congruences) {
        output << "(" << atom << ") ";
    }

    output << ")";
    return output;
}

