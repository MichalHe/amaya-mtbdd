#include "../include/lazy.hpp"

#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::pair;
using std::optional;


template <typename T>
std::string fmt_array(T* array, size_t array_size) {
    if (array_size == 0) return "[]";

    std::stringstream ss;
    ss << "[" << array[0];
    for (size_t i = 1; i < array_size; i++) {
        ss << ", " << array[i];
    }
    ss << "]";
    return ss.str();
}


inline int64_t dot_product(int64_t* coef, int64_t* symbol, size_t dimension) {
    assert(dimension > 0);
    int64_t dot = symbol[0] * coef[0];
    for (size_t i = 1; i < dimension; i++) {
        dot += symbol[i] * coef[i];
    }
    return dot;
}

s64 Presburger_Atom::compute_post_along_sym(s64 state, u64 symbol_bits) const {
    s64 post = 0;
    for (u64 var_i = 0; var_i < coefs.size(); var_i++) {
        s64 is_bit_set = (symbol_bits & (1 << var_i)) > 0;
        post += is_bit_set * coefs[var_i];
    }

    return (state - post) / 2;
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
    output << "State{";
    if (atom.constants.empty()) {
        output << "}";
        return output;
    }

    auto constants_it = atom.constants.begin();
    output << *constants_it;
    ++constants_it;

    for (; constants_it != atom.constants.end(); ++constants_it) {
        output << ", " << *constants_it;
    }

    output << "}";
    return output;
}

Variable_Bound_Analysis_Result* compute_bounds_analysis(Quantified_Atom_Conjunction& conj) {
    vector<Variable_Bounds> var_bounds {conj.var_count};
    /* unordered_set<u64> vars_with_both_bounds; */
    bool var_with_both_bounds_exists = false;

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

    return new Variable_Bound_Analysis_Result{.has_var_with_both_bounds = var_with_both_bounds_exists, .bounds = var_bounds};
}


Entaiment_Status compute_entailed_formula(FormulaPool& formula_pool, Conjuction_State& state) {
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
                return { .has_no_integer_solution = true, .removed_atom_count = 0, .state = std::nullopt };

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
            bool can_instantiate_var = false;
            s64 instantiation_bound_value;

            if (constant_vars[bound_var].first) {
                // The variable already has been instantiated, continuing will not add to the
                // uninstantiated variables in the entailed formula
                continue;
            }

            if (bounds.preferred_value == Preferred_Var_Value_Type::LOW) {
                if (bounds.lower.is_present) {
                    can_instantiate_var = true;
                    does_atom_imply_constant_var[bounds.lower.atom_idx] = true;
                    instantiation_bound_value = entailed_state.constants[bounds.lower.atom_idx];
                }
            } else if (bounds.preferred_value == Preferred_Var_Value_Type::HIGH) {
                if (bounds.upper.is_present) {
                    can_instantiate_var = true;
                    does_atom_imply_constant_var[bounds.upper.atom_idx] = true;
                    instantiation_bound_value = entailed_state.constants[bounds.upper.atom_idx];
                }
            }

            if (can_instantiate_var) {
                for (u64 atom_i = 0u; atom_i < entailed_formula.atoms.size(); atom_i++) {
                    auto& atom = entailed_formula.atoms[atom_i];

                    if (atom.coefs[bound_var])
                        std::cout << "Instantiating x" << bound_var << "=" << instantiation_bound_value << " into " << entailed_formula.atoms[atom_i].fmt_with_rhs(entailed_state.constants[atom_i]) << std::endl;

                    entailed_state.constants[atom_i] -= atom.coefs[bound_var] * instantiation_bound_value;
                    atom.coefs[bound_var] = 0;
                }
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
        if (is_bottom) return "\\bottom";
        else return "\\top";
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


void build_nfa_with_formula_entailement(FormulaPool& formula_pool, Conjuction_State& init_state) {
    unordered_set<Conjuction_State> seen_states; // Either the state was processed, or is in the queue
    vector<Conjuction_State> work_queue {init_state};

    u64 max_symbol = (1u << init_state.formula->var_count);
    vector<Conjuction_State> produced_states;
    while (!work_queue.empty()) {
        Conjuction_State state = work_queue.back();
        work_queue.pop_back();

        for (u64 symbol_bits = 0u; symbol_bits < max_symbol; symbol_bits++) {
            auto successor = state.successor_along_symbol(symbol_bits);

            auto entailment_status = compute_entailed_formula(formula_pool, successor);

            if (entailment_status.removed_atom_count > 0) {
                std::cout << "Was able to reduce formula `" << successor.formula->fmt_with_state(successor)
                          << " into " << entailment_status.state->formula->fmt_with_state(entailment_status.state.value()) << std::endl;
            }

            if (entailment_status.has_no_integer_solution) {
                successor = Conjuction_State{.formula=&formula_pool.bottom, .constants = {}};
            } else if (entailment_status.removed_atom_count > 0) {
                successor = entailment_status.state.value();
            }

            auto emplacement_result = seen_states.emplace(successor);
            if (!successor.formula->atoms.empty()) {
                if (emplacement_result.second) {
                    work_queue.push_back(successor);
                }
            }
        }
    }
}

const Quantified_Atom_Conjunction* FormulaPool::store_formula(Quantified_Atom_Conjunction& formula) {
    auto [it, did_insertion_happen] = formulae.emplace(formula);

    const Quantified_Atom_Conjunction* stored_formula_ptr = &(*it);
    if (did_insertion_happen && it->bounds_analysis_result == nullptr) {
        // It is OK to mutate this field as the bounds_analysis_result field is not hashed
        auto mut_formula = const_cast<Quantified_Atom_Conjunction*>(stored_formula_ptr);
        mut_formula->bounds_analysis_result = compute_bounds_analysis(formula);
    }

    return &(*it);
}


int main(void) {

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

    for (u64 var_i = 0u; var_i < real_formula.bounds_analysis_result->bounds.size(); var_i++) {
        auto& bounds = real_formula.bounds_analysis_result->bounds[var_i];
        std::cout << "x" << var_i << " has_lower=" << bounds.lower.is_present << " has_upper=" << bounds.upper.is_present << std::endl;
    }


    std::cout << real_formula.fmt_with_state(real_state) << std::endl;

    FormulaPool pool = FormulaPool();
    auto entailment_status = compute_entailed_formula(pool, real_state);
    std::cout << "Atoms removed: " << entailment_status.removed_atom_count << std::endl;
    std::cout << "Entailed formula: " << entailment_status.state.value().formula->fmt_with_state(entailment_status.state.value()) << std::endl;

    return 0;

    vector<Presburger_Atom> atoms = {
        Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0}),
        Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0}),
        Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 1}),
    };

    Quantified_Atom_Conjunction formula = {.atoms = atoms, .bound_vars = {1}, .var_count = 2};
    formula.bounds_analysis_result = compute_bounds_analysis(formula);

    // -y <= -1
    // x + y <= 1
    Conjuction_State state = Conjuction_State(&formula, {0, 1, 7});

    //FormulaPool pool = FormulaPool();
    //build_nfa_with_formula_entailement(pool, state);

    return 0;
}
