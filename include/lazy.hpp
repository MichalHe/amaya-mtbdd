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

    bool operator==(const Sized_Array<T>& other) const {
        if (other.size != this->size) return false;

        for (u64 i = 0; i < size; i++) {
            if (items[i] != other.items[i]) return false;
        }
        return true;
    }
};

std::size_t hash_array(Sized_Array<s64>& arr);

inline std::size_t hash_combine(std::size_t hash1, std::size_t hash2) {
    return hash1 + 0x9e3779b9 + (hash2 << 6) + (hash2 >> 2);
}

struct Congruence {
    Sized_Array<s64> coefs;
    s64 modulus_odd;
    s64 modulus_2pow;

    bool operator==(const Congruence& other) const {
        return coefs == other.coefs && modulus_odd == other.modulus_odd && modulus_2pow == other.modulus_2pow;
    }
};

struct Inequation {
    Sized_Array<s64> coefs;
    bool operator==(const Inequation& other) const {
        return coefs == other.coefs;
    }
};

struct Equation {
    Sized_Array<s64> coefs;
    bool operator==(const Equation& other) const {
        return coefs == other.coefs;
    }
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
    vector<Linear_Node>     linear_nodes; // First nodes are equations then inequations
    vector<Congruence_Node> congruence_nodes;
};

struct Quantified_Atom_Conjunction;
struct Formula_Allocator;

Dep_Graph build_dep_graph(const Quantified_Atom_Conjunction& conj);
Dep_Graph build_dep_graph(const Quantified_Atom_Conjunction& conj);
void write_dep_graph_dot(std::ostream& output, Dep_Graph& graph);
void identify_potential_variables(Dep_Graph& graph);
Conjunction_State simplify_graph_using_value(Dep_Graph& graph, Conjunction_State& state, u64 var, s64 val);
void simplify_graph_with_unbound_var(Dep_Graph& graph, u64 var);
bool simplify_graph(Dep_Graph& graph, Conjunction_State& state);
std::ostream& operator<<(std::ostream& output, const Congruence_Node& congruence_node);
std::ostream& operator<<(std::ostream& output, const Linear_Node& lin_node);


struct Quantified_Atom_Conjunction {
    // This requires certain order: Conjunction_State, conjunction state hashing
    Sized_Array<Congruence> congruences = {.items = nullptr, .size = 0};
    Sized_Array<Equation>   equations   = {.items = nullptr, .size = 0};
    Sized_Array<Inequation> inequations = {.items = nullptr, .size = 0};
    vector<u64> bound_vars;
    u64         var_count;
    Dep_Graph   dep_graph;
    bool        bottom = false;
    const Quantified_Atom_Conjunction* post_formula = nullptr;  // In case that the underlying modulus changes with post, we don't have to query formula pool

    Quantified_Atom_Conjunction() {}
    Quantified_Atom_Conjunction(bool top) : bottom(!top),
                                            var_count(0),
                                            bound_vars({}) {}
    Quantified_Atom_Conjunction(Sized_Array<Congruence> congruences, Sized_Array<Equation> eqs,
                                Sized_Array<Inequation> ineqs, const vector<u64>& bound_vars, u64 var_count) :
        congruences(congruences), equations(eqs), inequations(ineqs),
        bound_vars(bound_vars), var_count(var_count)
    {
        this->dep_graph = build_dep_graph(*this);
        identify_potential_variables(this->dep_graph);
    };

    bool operator==(const Quantified_Atom_Conjunction& other) const;

    u64 atom_count() const {
        return congruences.size + equations.size + inequations.size;
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
        for (auto& eq: formula.equations) {
            hash = hash_combine(hash, hash_array(eq.coefs));
        }

        for (auto& cg: formula.congruences) {
            hash = hash_combine(hash, hash_array(cg.coefs));
            hash = hash_combine(hash, cg.modulus_odd);
            hash = hash_combine(hash, cg.modulus_2pow);
        }
        for (auto& ineq: formula.inequations) {
            hash = hash_combine(hash, hash_array(ineq.coefs));
        }
        return hash;
    }
};

struct Conjunction_State {
    vector<s64> constants; // Congruences < Equations < Inequations

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

template <typename T>
struct Slab_Fragment {
    u64 next_free;
    u64 capacity;
    T*  items;

    Slab_Fragment(u64 capacity, T* items_buffer) : next_free(0), capacity(capacity), items(items_buffer) {};

    Slab_Fragment(u64 capacity) : capacity(capacity), next_free(0) {
        items = (T*) malloc(sizeof(Congruence) * capacity);
        assert(items);
    }

    Slab_Fragment(const Slab_Fragment& other) = delete;

    ~Slab_Fragment() {
        free(items);
    }
};


template <typename T>
Slab_Fragment<T> make_equaly_sized_fragment(Slab_Fragment<T>& buffer) {
    T* new_items = (T*) malloc(sizeof(T) * buffer.capacity);
    assert(new_items);
    return Slab_Fragment(buffer.capacity, new_items);
}


struct Allocator_Limits {
    u64 max_congruences, max_equations, max_inequations;
    u64 var_count;

    Allocator_Limits(const Quantified_Atom_Conjunction* root_formula):
        var_count(root_formula->var_count),
        max_congruences(root_formula->congruences.size),
        max_equations(root_formula->equations.size),
        max_inequations(root_formula->inequations.size)
    {};
};

struct Temporary_Space {
    Slab_Fragment<Congruence> congruences;  // Temporary space for congruences
    Slab_Fragment<Equation>   equations;    // Temporary space for equations
    Slab_Fragment<Inequation> inequations;  // Temporary space for inequations
    Slab_Fragment<s64>        coefs;        // Coefs for all the above temporary atoms

    Temporary_Space(const Quantified_Atom_Conjunction* root_formula) :
        congruences(root_formula->congruences.size),
        equations(root_formula->equations.size),
        inequations(root_formula->inequations.size),
        coefs(root_formula->var_count * (root_formula->congruences.size + root_formula->equations.size + root_formula->inequations.size)) {}
};

struct Formula_Allocator {
    u64 largest_congruence_block_size;
    vector<list<Slab_Fragment<Congruence>>> congruence_slabs;
    vector<list<Slab_Fragment<Equation>>>   equation_slabs;
    vector<list<Slab_Fragment<Inequation>>> inequation_slabs;

    list<s64*> coef_blocks;

    const u64 items_per_slab_buffer = 16;
    const u64 var_count;

    Temporary_Space tmp_space;
    Allocator_Limits limits;

    Formula_Allocator(const Quantified_Atom_Conjunction* root_formula) :
        largest_congruence_block_size(root_formula->congruences.size),
        var_count(root_formula->var_count),
        tmp_space(root_formula),
        limits(root_formula)
    {
        congruence_slabs.resize(limits.max_congruences);
        equation_slabs.resize(limits.max_equations);
        inequation_slabs.resize(limits.max_inequations);
    }

    // Formula_Allocator(const Allocator_Limits& limits, u64 var_count) :
    //     largest_congruence_block_size(limits.max_congruences),
    //     var_count(var_count),
    //     tmp_congruences(limits.max_congruences),
    //     tmp_equations(limits.max_equations),
    //     tmp_inequations(limits.max_inequations)
    // {
    //     congruence_slabs.resize(limits.max_congruences);

    //     // @Todo: Preallocate all slabs fragments with at least one buffer so that we can copy safely
    // }
    template <typename Atom_Type>
    Sized_Array<Atom_Type> alloc(list<Slab_Fragment<Atom_Type>>& slab_fragments, u64 count) {
        if (slab_fragments.empty()) {
            u64 slab_item_size = items_per_slab_buffer * count;
            slab_fragments.emplace_back(slab_item_size);
        }

        auto& last_slab = slab_fragments.back();

        if (last_slab.next_free >= last_slab.capacity) {
            slab_fragments.emplace_back(last_slab.capacity);
            last_slab = slab_fragments.back();
        }

        Atom_Type* allocated_atoms = &last_slab.items[last_slab.next_free];
        last_slab.next_free += count;

        return {.items = allocated_atoms, .size = count};
    }

    Sized_Array<Congruence> alloc_congruences(u64 count) {
        PRINTF_DEBUG("Allocating %lu congruence(s)\n", count);
        if (count == 0) return {.items = nullptr, .size = 0};

        auto& slabs_matching_size = congruence_slabs[count - 1];
        return alloc(slabs_matching_size, count);
    };

    Sized_Array<Equation> alloc_equations(u64 count) {
        PRINTF_DEBUG("Allocating %lu equation(s)\n", count);
        if (count == 0) return {.items = nullptr, .size = 0};

        auto& slabs_matching_size = equation_slabs[count - 1];
        return alloc(slabs_matching_size, count);
    }

    Sized_Array<Inequation> alloc_inequations(u64 count) {
        PRINTF_DEBUG("Allocating %lu equation(s)\n", count);
        if (count == 0) return {.items = nullptr, .size = 0};

        auto& slabs_matching_size = inequation_slabs[count - 1];
        return alloc(slabs_matching_size, count);
    }

    s64* alloc_coefs(u64 count) {
        s64* new_coef_block = new s64[count];
        coef_blocks.push_back(new_coef_block);
        return new_coef_block;
    }

    template <typename Atom_Type>
    Atom_Type* alloc_tmp_atom(Slab_Fragment<Atom_Type>& tmp_atom_storage) {
        u64 atom_slot = tmp_atom_storage.next_free;
        assert(atom_slot < tmp_atom_storage.capacity);
        tmp_atom_storage.next_free += 1;

        Atom_Type* atom = &tmp_atom_storage.items[atom_slot];

        s64* coef_block_start = tmp_space.coefs.items + tmp_space.coefs.next_free;
        tmp_space.coefs.next_free += limits.var_count;

        atom->coefs = {.items = coef_block_start, .size = limits.var_count};

        return atom;
    }

    Congruence* alloc_temporary_congruence() {
        return alloc_tmp_atom(this->tmp_space.congruences);
    }

    Equation* alloc_temporary_equation() {
        return alloc_tmp_atom(this->tmp_space.equations);
    }

    Inequation* alloc_temporary_inequation() {
        return alloc_tmp_atom(this->tmp_space.inequations);
    }

    void ensure_tmp_clean() { // Debug
        assert(tmp_space.coefs.next_free == 0);
        assert(tmp_space.congruences.next_free == 0);
        assert(tmp_space.equations.next_free == 0);
        assert(tmp_space.inequations.next_free == 0);
    }

    Quantified_Atom_Conjunction get_tmp_formula() {
        Quantified_Atom_Conjunction formula;
        formula.congruences = Sized_Array<Congruence>{.items = tmp_space.congruences.items, .size = tmp_space.congruences.next_free};
        formula.equations   = Sized_Array<Equation>{.items = tmp_space.equations.items, .size = tmp_space.equations.next_free};
        formula.inequations = Sized_Array<Inequation>{.items = tmp_space.inequations.items, .size = tmp_space.inequations.next_free};
        return formula;
    }

    Quantified_Atom_Conjunction commit_tmp_space() {
        Sized_Array<Congruence> congruences = alloc_congruences(tmp_space.congruences.next_free);
        Sized_Array<Equation>   equations   = alloc_equations(tmp_space.equations.next_free);
        Sized_Array<Inequation> inequations = alloc_inequations(tmp_space.inequations.next_free);

        u64 total_atoms = tmp_space.congruences.next_free + tmp_space.equations.next_free + tmp_space.inequations.next_free;
        s64* coefs = alloc_coefs(total_atoms * limits.var_count);

        memcpy(congruences.items, tmp_space.congruences.items, sizeof(Congruence) * tmp_space.congruences.next_free);
        memcpy(equations.items, tmp_space.equations.items, sizeof(Equation) * tmp_space.equations.next_free);
        memcpy(inequations.items, tmp_space.inequations.items, sizeof(Inequation) * tmp_space.inequations.next_free);

        assert(tmp_space.coefs.next_free == total_atoms * limits.var_count);
        memcpy(coefs, tmp_space.coefs.items, sizeof(s64) * tmp_space.coefs.next_free);

        for (Congruence& congruence: congruences) {
            congruence.coefs.items = coefs;
            coefs += limits.var_count;
        }

        for (Equation& equation: equations) {
            equation.coefs.items = coefs;
            coefs += limits.var_count;
        }

        for (Inequation& inequation: inequations) {
            inequation.coefs.items = coefs;
            coefs += limits.var_count;
        }

        this->drop_tmp();

        Quantified_Atom_Conjunction new_formula;
        new_formula.congruences = congruences;
        new_formula.equations   = equations;
        new_formula.inequations = inequations;

        return new_formula;
    }

    void drop_tmp() {  // "Free" the temporary storage
        tmp_space.congruences.next_free = 0;
        tmp_space.equations.next_free = 0;
        tmp_space.inequations.next_free = 0;
        tmp_space.coefs.next_free = 0;
    }

    ~Formula_Allocator() {
        for (auto coef_block: coef_blocks) delete[] coef_block;
    }
};



struct Formula_Pool { // Formula memory management
    unordered_set<Quantified_Atom_Conjunction> formulae;
    Formula_Allocator allocator;

    Formula_Pool(const Quantified_Atom_Conjunction* root_formula) : allocator(root_formula) {
        auto top = Quantified_Atom_Conjunction(true);
        auto bottom = Quantified_Atom_Conjunction(false);
        formulae.insert(top);
        formulae.insert(bottom);
    }

    const Quantified_Atom_Conjunction* store_formula(Quantified_Atom_Conjunction& formula);
    std::pair<const Quantified_Atom_Conjunction*, bool> store_formula_with_info(Quantified_Atom_Conjunction& formula);
};

typedef unordered_map<vector<s64>, list<Conjunction_State>> Prefix_Table;
std::ostream& operator<<(std::ostream& output, const map<const Quantified_Atom_Conjunction*, Prefix_Table>& macrostate);

struct Intermediate_Macrostate { // Macrostate that is being created, optimized for Pareto optimal insertions
    bool is_accepting = false;
    map<const Quantified_Atom_Conjunction*, Prefix_Table> formulae;
};

struct Macrostate_Entry_Family {
    const Quantified_Atom_Conjunction* formula;
    vector<Conjunction_State> states;

    bool operator==(const Macrostate_Entry_Family& other) const;
};

struct Finalized_Macrostate {
    vector<Macrostate_Entry_Family> entries; // In sorted order by formula pointeres.
    bool is_accepting = false;
    u64 handle;

    bool operator==(const Finalized_Macrostate& other) const;
};

template <>
struct std::hash<Finalized_Macrostate> {
    std::size_t operator() (const Finalized_Macrostate& macrostate) const {
        std::size_t hash = 0u;

        // Hash the states constituting macrostate
        for (auto& macrostate_entry : macrostate.entries) {
            std::size_t formula_hash = std::hash<const Quantified_Atom_Conjunction*>{}(macrostate_entry.formula);

            // @Simplicity: Is hashing really useful in this case?
            for (auto& state: macrostate_entry.states) {
                std::size_t state_hash = 0u;

                for (auto state_constant: state.constants) {
                    std::size_t state_constant_hash = ((state_constant >> 16) ^ state_constant) * 0x45d9f3b;
                    state_constant = ((state_constant >> 16) ^ state_constant) * 0x45d9f3b;
                    state_constant = (state_constant >> 16) ^ state_constant;

                    hash = hash_combine(hash, formula_hash);
                }

                // @Note: hash combination must be order independent here, because we do not impose any ordering on the states
                formula_hash ^= state_hash;
            }

            hash = hash_combine(hash, formula_hash);
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
    vector<Finalized_Macrostate>& output_queue;
    unordered_map<Finalized_Macrostate, u64> known_macrostates;
    unordered_set<u64> accepting_macrostates;

    bool  is_trap_state_needed;
    State trap_state_handle;
    bool  is_top_state_needed;
    State topstate_handle;
};


struct Stateful_Formula {
    Conjunction_State state;
    Quantified_Atom_Conjunction formula;
};
Stateful_Formula convert_graph_into_formula(Dep_Graph& graph, Formula_Allocator& allocator, Conjunction_State& state);
std::pair<const Quantified_Atom_Conjunction*, Conjunction_State> simplify_stateful_formula(const Quantified_Atom_Conjunction* formula, Conjunction_State& state, Formula_Pool& pool);

enum class Bound_Type : unsigned {
    NONE  = 0x00,
    UPPER = 0x01,
    LOWER = 0x02,
};

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

void insert_into_post_if_valuable2(Intermediate_Macrostate& post, const Quantified_Atom_Conjunction* formula, Conjunction_State& successor);

template <typename T>
Sized_Array<T> make_sized_array(const vector<T>& items) {
    T* buf = new T[items.size()];
    u64 item_idx = 0;
    for (auto& item: items) {
        buf[item_idx] = item;
        item_idx += 1;
    }
    return Sized_Array<T> {.items = buf, .size = items.size()};
}

#endif

