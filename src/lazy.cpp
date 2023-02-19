#include "../include/lazy.hpp"

#include <stdlib.h>

#include <bitset>
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

typedef Quantified_Atom_Conjunction Formula;

s64 Presburger_Atom::compute_post_along_sym(s64 state, u64 symbol_bits) const {
    s64 dot = 0;
    for (u64 var_i = 0; var_i < coefs.size(); var_i++) {
        s64 is_bit_set = (symbol_bits & (1u << var_i)) > 0;
        dot += is_bit_set * coefs[var_i];
    }

    if (this->type == PR_ATOM_CONGRUENCE) {
        s64 post = state - dot;
        post += modulus * ((state - dot % 2) != 0);
        post /= 2;
        post = post % modulus;
        return post;
    }

    s64 post = (state - dot);

    s64 post_div_2 = post / 2;
    s64 post_mod_2 = post % 2;

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
        return (state + dot) <= 0;
    }
}


bool Conjuction_State::operator==(const Conjuction_State& other) const {
    return formula == other.formula && constants == other.constants;
}

bool Presburger_Atom::operator==(const Presburger_Atom& other) const {
    return coefs == other.coefs;
}

bool Quantified_Atom_Conjunction::operator==(const Quantified_Atom_Conjunction& other) const {
    return other.is_bottom == is_bottom && other.var_count == var_count && other.atoms == atoms && other.bound_vars == other.bound_vars;
}

Conjuction_State Conjuction_State::successor_along_symbol(u64 symbol) {
    vector<s64> successor_values(constants.size());
    for (u64 atom_i = 0; atom_i < formula->atoms.size(); atom_i++) {
        s64 current_state = constants[atom_i];
        successor_values[atom_i] = formula->atoms[atom_i].compute_post_along_sym(current_state, symbol);
    }
    return Conjuction_State {.formula = formula, .constants = successor_values};
}

bool Conjuction_State::accepts_last_symbol(u64 symbol) {
    bool accepts = true;
    for (u64 atom_i = 0; atom_i < formula->atoms.size(); atom_i++) {
        s64 current_state = constants[atom_i];
        accepts &= formula->atoms[atom_i].compute_post_along_sym(current_state, symbol);
    }
    return accepts;
}

void Conjuction_State::post(unordered_set<Conjuction_State>& known_states, vector<Conjuction_State>& new_states) {
    const s64 max_symbol_val = (1 << formula->var_count);

    vector<s64> post(constants.size());

    for (s64 symbol_bits = 0; symbol_bits < max_symbol_val; symbol_bits++) {

        for (u64 atom_i = 0; atom_i < formula->atoms.size(); atom_i++) {
            s64 current_state = constants[atom_i];
            post[atom_i] = formula->atoms[atom_i].compute_post_along_sym(current_state, symbol_bits);
        }

        Conjuction_State dest = {.formula = formula, .constants = post};

        auto insertion_result = known_states.emplace(dest);
        auto did_insertion_happen = insertion_result.second;
        if (did_insertion_happen) {
            new_states.push_back(dest);
        }
    }
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

std::ostream& operator<<(std::ostream& output, const Conjuction_State& atom) {
    output << "(";
    if (!atom.constants.empty()) {
        auto constants_it = atom.constants.begin();
        output << *constants_it;
        ++constants_it;

        for (; constants_it != atom.constants.end(); ++constants_it) {
            output << ", " << *constants_it;
        }
    }

    output << ") [" << atom.formula << "]";
    return output;
}

bool Structured_Macrostate::operator==(const Structured_Macrostate& other) const {
    if (is_accepting != other.is_accepting) return false;
    return formulae == other.formulae;
};


Variable_Bound_Analysis_Result* compute_bounds_analysis(Quantified_Atom_Conjunction& conj) {
    vector<Variable_Bounds> var_bounds {conj.var_count};
    bool var_with_both_bounds_exists = false;

    vector<vector<u64>> congruences_per_var (conj.var_count);

    for (u64 atom_i = 0u; atom_i < conj.atoms.size(); atom_i++) {
        auto& atom = conj.atoms[atom_i];

        u64 bounded_var_count = 0u;
        u64 bounded_var_i;
        for (u64 var_i = 0u; var_i < conj.var_count; var_i++) {
            if (atom.coefs[var_i]) {
                if (bounded_var_count && atom.type == PR_ATOM_INEQ) {
                    // Make sure we compute preferred values only for atoms with more than one variable
                    auto preferred_value = atom.coefs[var_i] > 0 ? Preferred_Var_Value_Type::LOW : Preferred_Var_Value_Type::HIGH;
                    var_bounds[var_i].preferred_value = var_bounds[var_i].preferred_value | preferred_value;

                    if (bounded_var_count == 1) { // Update also the info for the first encountered variable
                        auto preferred_value = atom.coefs[bounded_var_i] > 0 ? Preferred_Var_Value_Type::LOW : Preferred_Var_Value_Type::HIGH;
                        var_bounds[bounded_var_i].preferred_value = var_bounds[bounded_var_i].preferred_value | preferred_value;
                    }
                }
                if (atom.type == PR_ATOM_CONGRUENCE) {
                    congruences_per_var[var_i].push_back(atom_i);
                }

                bounded_var_count++;
                bounded_var_i = var_i;
            }
        }

        if (bounded_var_count == 1) {
            if (atom.type == Presburger_Atom_Type::PR_ATOM_EQ) {
                auto var_coef = atom.coefs[bounded_var_i];

                var_with_both_bounds_exists = true;
                /* vars_with_both_bounds.insert(bounded_var_i); */

                var_bounds[bounded_var_i].lower.is_present = true;
                var_bounds[bounded_var_i].lower.atom_idx = atom_i;
                var_bounds[bounded_var_i].upper.is_present = true;
                var_bounds[bounded_var_i].upper.atom_idx = atom_i;

            } else if (atom.type == Presburger_Atom_Type::PR_ATOM_INEQ) {
                auto var_coef = atom.coefs[bounded_var_i];

                Variable_Bound *bound_to_set, *other_bound;
                if (var_coef > 0) {
                    bound_to_set = &var_bounds[bounded_var_i].upper;
                    other_bound  = &var_bounds[bounded_var_i].lower;
                } else { // -x <= 10  ===  x >= -10
                    bound_to_set = &var_bounds[bounded_var_i].lower;
                    other_bound  = &var_bounds[bounded_var_i].upper;
                }

                // @Robustness: Deal with multiple bounds of the same type (lower/upper) on the same value
                bound_to_set->is_present = true;
                bound_to_set->atom_idx = atom_i;

                if (other_bound->is_present) {
                    var_with_both_bounds_exists = true;
                    /* vars_with_both_bounds.insert(bounded_var_i); */
                }
            }
        }
    }

    return new Variable_Bound_Analysis_Result{
        .has_var_with_both_bounds = var_with_both_bounds_exists,
        .bounds = var_bounds,
        .congruences_per_var = congruences_per_var
    };
}


std::optional<std::pair<s64, u64>> compute_instantiating_value_for_var(const Quantified_Atom_Conjunction& formula,
                                                                       Conjuction_State& state,
                                                                       u64 var_to_instantiate,
                                                                       Variable_Bounds& bounds,
                                                                       vector<u64>& congruences_referencing_var)
{
    bool do_bounds_allow_instantiation = false;
    u64 used_atom_i = 0;

    if (bounds.preferred_value == Preferred_Var_Value_Type::LOW && bounds.lower.is_present) {
        do_bounds_allow_instantiation = true;
        used_atom_i = bounds.lower.atom_idx;
    }
    else if (bounds.preferred_value == Preferred_Var_Value_Type::HIGH && bounds.upper.is_present) {
        do_bounds_allow_instantiation = true;
        used_atom_i = bounds.upper.atom_idx;
    }

    if (!do_bounds_allow_instantiation || congruences_referencing_var.size() > 1)
        return std::nullopt;

    if (congruences_referencing_var.empty())
        return std::make_pair(state.constants[used_atom_i], used_atom_i);

    auto& congruence = formula.atoms[congruences_referencing_var[0]];
    s64 bound_point = state.constants[used_atom_i];

    bool found_congruence_var = false;
    u64 var_in_congruence;
    for (u64 var_i = 0u; var_i < congruence.coefs.size(); ++var_i) {
        if (congruence.coefs[var_i]) {
            if (found_congruence_var) {
                // @Research: If there are multiple variables in the congruence no idea what value shuold the var have.
                return std::nullopt;
            }
            found_congruence_var = true;
            var_in_congruence = var_i;
        }
    }

    assert(var_in_congruence == var_to_instantiate);

    // Find a solution in the unshifted congruence range, e.g., -y == 303 (mod 299_993)
    s64 unshifted_solution = state.constants[congruences_referencing_var[0]] / congruence.coefs[var_in_congruence];
    assert (unshifted_solution * congruence.coefs[var_in_congruence] == state.constants[congruences_referencing_var[0]]);
    unshifted_solution += (unshifted_solution < 0)*congruence.modulus;

    s64 shift_coef = bound_point / congruence.modulus;
    shift_coef -= (bound_point < 0); // Essentailly a signed floor division
    s64 congruence_shift = congruence.modulus * shift_coef;

    s64 instantiated_value = unshifted_solution + congruence_shift;
    return std::make_pair(instantiated_value, used_atom_i);
}


Entaiment_Status compute_entailed_formula(Formula_Pool& formula_pool, Conjuction_State& state) {
    auto formula = state.formula;

    if (!formula->bounds_analysis_result->has_var_with_both_bounds)
        return { .has_no_integer_solution = false, .removed_atom_count = 0, .state = std::nullopt };

    u64 atoms_implying_constant_var = 0;
    vector<bool> does_atom_imply_constant_var (formula->atoms.size());
    vector<pair<bool, s64>> constant_vars (formula->var_count);

    bool has_no_integer_solution = false;

    // Check whether any of the atoms in variable bounds imply a constant value
    auto& var_bounds = formula->bounds_analysis_result->bounds;
    for (u64 var_i = 0u; var_i < formula->var_count; var_i++) {
        auto& lower_bound = var_bounds[var_i].lower, upper_bound = var_bounds[var_i].upper;
        auto& lower_bound_atom = formula->atoms[lower_bound.atom_idx];
        auto& upper_bound_atom = formula->atoms[upper_bound.atom_idx];

        if (lower_bound.is_present && upper_bound.is_present) {
            s64 lower_bound_val = state.constants[lower_bound.atom_idx] / lower_bound_atom.coefs[var_i];
            s64 upper_bound_val = state.constants[upper_bound.atom_idx] / upper_bound_atom.coefs[var_i];

            if (lower_bound_val > upper_bound_val)
                return { .has_no_integer_solution = true, .removed_atom_count = 0, .state = Conjuction_State(&formula_pool.bottom, {}) };

            if (lower_bound_val == upper_bound_val) {
                constant_vars[var_i] = {true, lower_bound_val};

                if (!does_atom_imply_constant_var[lower_bound.atom_idx]) {
                    does_atom_imply_constant_var[lower_bound.atom_idx] = true;
                    atoms_implying_constant_var++;
                };
                if (!does_atom_imply_constant_var[upper_bound.atom_idx]) {
                    does_atom_imply_constant_var[upper_bound.atom_idx] = true;
                    atoms_implying_constant_var++;
                };
            }
        }
    }

    // @Robustness: We should also handle formulae without any hard bounds on them, e.g., (exists ((x Int) (y Int))
    //              (<= (+ x y) 0)). Handling it should be probably done when we shift from simple analysis to inter-
    //              atom influence analysis.

    if (!atoms_implying_constant_var)  // The current formula does not imply a constant variable value
        return { .has_no_integer_solution = false, .removed_atom_count = 0, .state = state};

    // We will remove at least a single atom implying a constant - a new formula will be produced
    Quantified_Atom_Conjunction entailed_formula = Quantified_Atom_Conjunction(*formula);
    entailed_formula.bounds_analysis_result = nullptr;
    Conjuction_State entailed_state(state);

    for (u64 atom_i = 0u; atom_i < formula->atoms.size(); atom_i++) {
        if (does_atom_imply_constant_var[atom_i]) continue;

        // Check whether the atom references any constant vars; if yes, substitute them into the atom
        auto& atom = entailed_formula.atoms[atom_i];
        for (u64 var_i = 0u; var_i < formula->var_count; var_i++) {
            if (constant_vars[var_i].first && atom.coefs[var_i]) {
                entailed_state.constants[var_i] -= constant_vars[var_i].second;
                atom.coefs[var_i] = 0;
            }
        }
    }

    // Try instantiating the quantifier if there is a clear value for the quantified variables to pick. E.g, in
    // \exists x (x <= 1 \land x + y >= 10) x should be 1 because x = 0 implies the sol(y) = {10, 11, 12}, whereas
    // x = 1 implies sol(y) = {9, 10, 11, 12} and we want to have the free variables unrestricted as much as possible
    {
        vector<u64> quantified_vars_not_instantiated;
        for (auto bound_var: formula->bound_vars) {
            auto& bounds = var_bounds[bound_var];

            if (constant_vars[bound_var].first) {
                // The variable already has been instantiated, continuing will not add to the
                // uninstantiated variables in the entailed formula
                continue;
            }

            auto& congruences_for_var = formula->bounds_analysis_result->congruences_per_var[bound_var];
            auto maybe_instantiation_result = compute_instantiating_value_for_var(entailed_formula, state, bound_var, bounds, congruences_for_var);

            if (maybe_instantiation_result.has_value()) {
                auto [instantiation_value, used_atom_i] = maybe_instantiation_result.value();
                for (u64 atom_i = 0u; atom_i < entailed_formula.atoms.size(); atom_i++) {
                    auto& atom = entailed_formula.atoms[atom_i];
                    entailed_state.constants[atom_i] -= atom.coefs[bound_var] * instantiation_value;
                    atom.coefs[bound_var] = 0;
                }

                does_atom_imply_constant_var[used_atom_i] = true;

            } else {
                quantified_vars_not_instantiated.push_back(bound_var);
            }
        }
        entailed_formula.bound_vars = quantified_vars_not_instantiated;
    }

    for (u64 atom_i = 0u ; atom_i < entailed_formula.atoms.size(); atom_i++) {
        auto& atom = entailed_formula.atoms[atom_i];
        bool all_coefs_are_zero = true;
        for (auto& coef: atom.coefs) {
            if (coef) all_coefs_are_zero = false;
        }
        if (all_coefs_are_zero) {
            does_atom_imply_constant_var[atom_i] = true;
            atoms_implying_constant_var++;
        }
    }

    // Remove atoms that have implied a constant variable value
    {
        u64 entailed_conj_size = entailed_formula.atoms.size() - atoms_implying_constant_var;
        vector<Presburger_Atom> new_atoms(entailed_conj_size);
        vector<s64> relevant_atom_constants(entailed_conj_size);
        u64 new_atom_insertion_cursor = 0;
        for (u64 atom_i = 0u; atom_i < entailed_formula.atoms.size(); atom_i++) {
            if (does_atom_imply_constant_var[atom_i]) continue;

            new_atoms[new_atom_insertion_cursor] = entailed_formula.atoms[atom_i];
            relevant_atom_constants[new_atom_insertion_cursor] = entailed_state.constants[atom_i];
            new_atom_insertion_cursor++;
        }

        entailed_formula.atoms = new_atoms;
        entailed_state.constants = relevant_atom_constants;
    }

    entailed_state.formula = formula_pool.store_formula(entailed_formula);

    return Entaiment_Status {
        .has_no_integer_solution = false,
        .removed_atom_count = atoms_implying_constant_var,
        .state = entailed_state,
    };
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

std::string Quantified_Atom_Conjunction::fmt_with_state(Conjuction_State& state) const {
    if (atoms.empty()) {
        if (is_bottom) return "false";
        else return "true";
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


void insert_successor_into_post_if_valuable(Structured_Macrostate& post, Conjuction_State& successor) {
    auto formula = successor.formula;
    auto& bucket = post.formulae[successor.formula];

    bool should_be_inserted = true;
    bool insert_position_found = false;

    /*
     * @Optimize: The bucket is kept in a sorted order, so once we arrive at a state that is lexicographically
     *            larger, we don't have to iterate further.
     */

    auto bucket_iter = bucket.begin();
    list<Conjuction_State>::iterator insert_position;
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

std::ostream& operator<<(std::ostream& output, const Structured_Macrostate& macrostate) {
    output << "{";
    for (auto& [formula, states]: macrostate.formulae) {
        auto state_it = states.begin();
        output << *state_it;
        ++state_it;

        for (; state_it != states.end(); ++state_it) {
            output << ", " << *state_it;
        }

    }
    if (macrostate.is_accepting) {
        output << " (+ ACCEPTING)";
    }
    output << "}";
    return output;
}

void explore_macrostate(Structured_Macrostate& macrostate,
                        Alphabet_Iterator& alphabet_iter,
                        Formula_Pool& formula_pool,
                        unordered_map<Structured_Macrostate, u64>& known_macrostates,
                        unordered_set<u64>& accepting_macrostates,
                        vector<Structured_Macrostate>& output_queue)
{
    alphabet_iter.reset();

    std::cout << "Exploring: " << macrostate << std::endl;

    while (!alphabet_iter.finished) {
        Structured_Macrostate post;
        bool is_accepting = false;
        for (u64 symbol = alphabet_iter.init_quantif(); alphabet_iter.has_more_quantif_symbols; symbol = alphabet_iter.next_symbol()) {
            for (auto& [formula, states]: macrostate.formulae) {
                for (auto& state: states) {
                    auto successor = state.successor_along_symbol(symbol);
                    auto entailment_result = compute_entailed_formula(formula_pool, successor);

                    if (entailment_result.state.has_value()) {
                        successor = entailment_result.state.value();
                    }

                    if (successor.formula == &formula_pool.top) {
                        assert(false); // @TODO: Make the entire post lead to \\top here
                        break;
                    } else if (successor.formula != &formula_pool.bottom) {
                        insert_successor_into_post_if_valuable(post, successor);
                    }

                    post.is_accepting |= state.accepts_last_symbol(symbol);
                }
            }
        }

        // Assign a unique integer to every state
        post.handle = known_macrostates.size();
        auto [element_iter, was_inserted] = known_macrostates.emplace(post, post.handle);
        if (was_inserted) {
            // @Simplicity: Maybe the known_macrostates should be a set instead of a map since we are storing the handle
            //              inside the macrostate either way.
            output_queue.push_back(post);

            if (post.is_accepting) { // Do the hash-query only if we see the macrostate for the first time
                accepting_macrostates.emplace(post.handle);
            }
        } else {
            post.handle = element_iter->second;  // Use the already existing handle
        }
        
        //NFA nfa;
        //nfa.add_transition(macrostate.handle, post.handle, 0u, 0u);
    }
}


// @Cleanup: Move this into base.cpp
void NFA::add_transition(State from, State to, u64 symbol, u64 quantified_bits_mask) {
    u8 rich_symbol[var_count];

    for (u64 bit_i = 0u; bit_i < var_count; bit_i++) {
        u64 current_bit = (1u << bit_i);
        bool dont_care = (current_bit & quantified_bits_mask) > 0;
        u8 care_val = (symbol & current_bit > 0);

        rich_symbol[bit_i] = dont_care ? 2 : care_val;
    }
    
    set<State> leaf_states = {to};
    Transition_Destination_Set leaf_contents = {.destination_set = &leaf_states};
    sylvan::MTBDD leaf = make_set_leaf(&leaf_contents);

    sylvan::MTBDD mtbdd = sylvan::mtbdd_cube(this->vars, rich_symbol, leaf);
}


void build_nfa_with_formula_entailement(Formula_Pool& formula_pool, Conjuction_State& init_state) {
    u64 max_symbol = (1u << init_state.formula->var_count);
    vector<Conjuction_State> produced_states;

    // Prepare bit masks and constants to iterate over quantified symbols efficiently
    u64 quantified_vars_mask = 0u;
    for (auto quantified_var: init_state.formula->bound_vars) {
        quantified_vars_mask |= (1u << quantified_var);
    }
    const u64 quantified_symbols_cnt = 1u << init_state.formula->bound_vars.size();
    const u64 free_symbols_cnt = 1u << (init_state.formula->var_count - init_state.formula->bound_vars.size());

    // Use map instead of unordered_map so that we can serialize its contents in an canonical fashion
    Structured_Macrostate post;

    unordered_map<Structured_Macrostate, u64> known_macrostates;
    unordered_set<u64> accepting_macrostates;

    vector<Structured_Macrostate> work_queue;
    {
        list<Conjuction_State> init_list = {init_state};
        map<const Formula*, list<Conjuction_State>> init_macrostate_formulae = { {init_state.formula, init_list} };
        Structured_Macrostate init_macrostate = { .is_accepting = false, .formulae = init_macrostate_formulae };
        work_queue.push_back(init_macrostate);

        known_macrostates.emplace(init_macrostate, known_macrostates.size());
    }

    Alphabet_Iterator alphabet_iter = Alphabet_Iterator(init_state.formula->var_count, init_state.formula->bound_vars);
    while (!work_queue.empty()) {
        auto macrostate = work_queue.back();
        work_queue.pop_back();

        explore_macrostate(macrostate, alphabet_iter, formula_pool, known_macrostates, accepting_macrostates, work_queue);
    }

    std::cout << "Known macrostates:" << std::endl;
    for (auto& [state, handle]: known_macrostates) {
        std::cout << "@" << handle << " :: " << state << std::endl;
    }

    std::cout << "Accepting macrostates: " << std::endl;
    for (auto macrostate_handle: accepting_macrostates) {
        std::cout << "  @" << macrostate_handle << std::endl;
    }

}

const Quantified_Atom_Conjunction* Formula_Pool::store_formula(Quantified_Atom_Conjunction& formula) {
    auto [it, did_insertion_happen] = formulae.emplace(formula);

    const Quantified_Atom_Conjunction* stored_formula_ptr = &(*it);
    if (did_insertion_happen && it->bounds_analysis_result == nullptr) {
        // It is OK to mutate this field as the bounds_analysis_result field is not hashed
        auto mut_formula = const_cast<Quantified_Atom_Conjunction*>(stored_formula_ptr);
        mut_formula->bounds_analysis_result = compute_bounds_analysis(formula);
    }

    return &(*it);
}

void test_nfa_construction_without_quantifiers() {
    Quantified_Atom_Conjunction formula = {
        .atoms = {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1}),
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0}, 3),
        },
        .bound_vars = {},
        .var_count = 2
    };

    Formula_Pool pool = Formula_Pool();
    auto formula_id = pool.store_formula(formula);

    Conjuction_State init_state = Conjuction_State{.formula = formula_id, .constants = {-1, 0}};

    std::cout << "Input formula   : " << formula.fmt_with_state(init_state) << std::endl;
    build_nfa_with_formula_entailement(pool, init_state);
}


void test_nfa_construction_with_single_quantifier() {
    Quantified_Atom_Conjunction formula = {
        .atoms = {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 1}),
        },
        .bound_vars = {1},
        .var_count = 2
    };

    Formula_Pool pool = Formula_Pool();
    auto formula_id = pool.store_formula(formula);

    Conjuction_State init_state = Conjuction_State{.formula = formula_id, .constants = {0}};

    std::cout << "Input formula   : " << formula.fmt_with_state(init_state) << std::endl;
    build_nfa_with_formula_entailement(pool, init_state);
}


void test_constr_on_real_formula() {
    /* (exists ((y Int), (m Int))
     *   (land
     *     (<= (+ x (- y))  -1)
     *     (<= (+ m (- z))  -1)
     *     (<= y -1)
     *     (<= (- m) 0)
     *     (<= m 0)
     *     (= (mod (+ m (- y)) 299_993) 303)
     *   )
     * )
     *
     * Variables are renamed as: m -> x0, x -> x1, y -> x2, z -> x3
     */
    Quantified_Atom_Conjunction real_formula = {
        .atoms = {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, -1, 0}),        // (<= (+ x (- y))  -1)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, -1}),        // (<= (+ m (- z))  -1)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 1, 0}),         // (<= y -1)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),        // (<= (- m) 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, 0}),         // (<= m 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1, 0}, 299993),  // (= (mod (+ m (- y)) 299_993) 303)
        },
        .bound_vars = {0, 2},
        .var_count = 4
    };

    real_formula.bounds_analysis_result = compute_bounds_analysis(real_formula);
    Conjuction_State real_state = Conjuction_State{.formula = &real_formula, .constants = {-1, -1, -1, 0, 0, 303}};

    std::cout << "Input formula   : " << real_formula.fmt_with_state(real_state) << std::endl;

    Formula_Pool pool = Formula_Pool();
    auto entailment_status = compute_entailed_formula(pool, real_state);
    std::cout << "Atoms removed   : " << entailment_status.removed_atom_count << std::endl;
    std::cout << "Entailed formula: " << entailment_status.state.value().formula->fmt_with_state(entailment_status.state.value()) << std::endl;

    build_nfa_with_formula_entailement(pool, real_state);
}


int main(void) {
    test_nfa_construction_with_single_quantifier();

    return 0;
}
