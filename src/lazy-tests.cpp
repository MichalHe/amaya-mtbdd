#define DOCTEST_CONFIG_IMPLEMENT
#include "../external/doctest.h"
#include "../include/lazy.hpp"
#include "../include/base.hpp"

#include <algorithm>
#include <vector>
#include <unordered_map>

typedef Quantified_Atom_Conjunction Formula;

using std::pair;
using std::vector;
using std::unordered_map;


#define DEBUG_TESTS 0
#if DEBUG_TESTS
#define TEST_PRINT_DEBUG(it) do { std::cerr << it << std::endl; } while (0)
#define TEST_PRINTF_DEBUG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define TEST_PRINT_DEBUG(it)
#define TEST_PRINTF_DEBUG(...)
#endif


struct Atom_Allocator {
    vector<s64*> coef_blocks;
    vector<Congruence*> congruence_blocks;
    vector<Equation*>   equation_blocks;
    vector<Inequation*> inequation_blocks;

    ~Atom_Allocator() {
        for (s64* coef_block: coef_blocks) delete[] coef_block;
        for (Congruence* congruence_block: congruence_blocks) delete[] congruence_block;
        for (Equation* equation_block: equation_blocks) delete[] equation_block;
        for (Inequation* inequation_block: inequation_blocks) delete[] inequation_block;
    }
};

Sized_Array<Congruence> alloc_congruences(Atom_Allocator& alloc, const vector<pair<vector<s64>, s64>>& congruence_coefs) {
    if (congruence_coefs.empty()) return {.items = nullptr, .size = 0};
    u64 var_count = congruence_coefs[0].first.size();

    s64* coef_block = new s64[var_count * congruence_coefs.size()];
    alloc.coef_blocks.push_back(coef_block);
    u64 next_free_coef_slot = 0;

    Congruence* congruence_block = new Congruence[congruence_coefs.size()];
    alloc.congruence_blocks.push_back(congruence_block);
    u64 next_free_congruence_slot = 0;

    for (auto& [coefs, modulus]: congruence_coefs) {
        // Make a sized array out of coefficients
        s64* congruence_coef_block = coef_block + next_free_coef_slot * var_count;

        for (auto coef: coefs) {
            coef_block[next_free_coef_slot] = coef;
            next_free_coef_slot += 1;
        }

        Sized_Array<s64> coefs_block = {.items = congruence_coef_block, .size=var_count};

        // Make a congruence and push it back
        auto [modulus_odd, modulus_pow2] = decompose_modulus(modulus);
        Congruence congruence = {.coefs = coefs_block, .modulus_odd = modulus_odd, .modulus_2pow = modulus_pow2};
        congruence_block[next_free_congruence_slot] = congruence;
        next_free_congruence_slot += 1;
    }

    return {.items = congruence_block, .size = congruence_coefs.size()};
};

Sized_Array<Equation> alloc_equations(Atom_Allocator& alloc, const vector<vector<s64>>& equation_coefs) {
    if (equation_coefs.empty()) return {.items = nullptr, .size = 0};

    u64 var_count = equation_coefs[0].size();

    s64* coef_block = new s64[var_count * equation_coefs.size()];
    alloc.coef_blocks.push_back(coef_block);
    u64 next_free_coef_slot = 0;

    Equation* equation_block = new Equation[equation_coefs.size()];
    alloc.equation_blocks.push_back(equation_block);
    u64 next_free_equation_slot = 0;

    for (auto& coefs: equation_coefs) {
        s64* eq_coef_block = coef_block + next_free_equation_slot * var_count;
        for (auto coef: coefs) {
            coef_block[next_free_coef_slot] = coef;
            next_free_coef_slot += 1;
        }
        Sized_Array<s64> coefs_block = {.items = eq_coef_block, .size=var_count};

        Equation equation = {.coefs = coefs_block};
        equation_block[next_free_equation_slot] = equation;
        next_free_equation_slot += 1;
        .congruences
    }

    return {.items = equation_block, .size = equation_coefs.size()};
};

Sized_Array<Inequation> alloc_inequations(Atom_Allocator& alloc, const vector<vector<s64>>& ineq_coefs) {
    if (ineq_coefs.empty()) return {.items = nullptr, .size = 0};

    u64 var_count = ineq_coefs[0].size();

    s64* coef_block = new s64[var_count * ineq_coefs.size()];
    alloc.coef_blocks.push_back(coef_block);
    u64 next_free_coef_slot = 0;

    Inequation* inequation_block = new Inequation[ineq_coefs.size()];
    alloc.inequation_blocks.push_back(inequation_block);
    u64 next_free_inequation_slot = 0;

    for (auto& coefs: ineq_coefs) {
        s64* eq_coef_block = coef_block + var_count * next_free_inequation_slot;
        for (auto coef: coefs) {
            coef_block[next_free_coef_slot] = coef;
            next_free_coef_slot += 1;
        }
        Sized_Array<s64> coefs_block = {.items = eq_coef_block, .size=var_count};

        Inequation ineq = {.coefs = coefs_block};
        inequation_block[next_free_inequation_slot] = ineq;
        next_free_inequation_slot += 1;
    }

    return {.items = inequation_block, .size = ineq_coefs.size()};
};

Formula make_formula(Atom_Allocator& allocator,
                     const vector<pair<vector<s64>, s64>>& congruences_coefs,
                     const vector<vector<s64>>& eqs,
                     const vector<vector<s64>>& ineqs,
                     const vector<u64> bound_vars) {
    Sized_Array<Congruence> congruences = alloc_congruences(allocator, congruences_coefs);
    Sized_Array<Equation> equations     = alloc_equations(allocator, eqs);
    Sized_Array<Inequation> inequations = alloc_inequations(allocator, ineqs);

    u64 var_count = 0;
    if (!congruences_coefs.empty()) var_count = congruences_coefs[0].first.size();
    else if (!ineqs.empty())        var_count = ineqs[0].size();
    else if (!eqs.empty())          var_count = eqs[0].size();

    return Formula(congruences, equations, inequations, bound_vars, var_count);
}

void assert_dfas_are_isomorphic(const NFA& expected, const NFA& actual) {
    CHECK(expected.initial_states.size() == 1);
    CHECK(actual.initial_states.size() == 1);

    CHECK(expected.states.size() == actual.states.size());
    CHECK(expected.final_states.size() == actual.final_states.size());

    unordered_map<State, State> isomorphism {{*expected.initial_states.begin(), *actual.initial_states.begin()}};

    vector<State> work_queue {*expected.initial_states.begin()};
    while (!work_queue.empty()) {
        auto expected_state = work_queue.back();
        auto actual_state_it = isomorphism.find(expected_state);

        work_queue.pop_back();

        // Make sure that the (expected, actual) mapping has been inserted when the expected state has been discovered
        CHECK(actual_state_it != isomorphism.end());
        auto actual_state = actual_state_it->second;
        TEST_PRINTF_DEBUG("Checking transitions between states: %lu (expected) <> %lu (actual)\n", expected_state, actual_state);

        auto expected_transitions = expected.get_symbolic_transitions_for_state(expected_state);
        auto actual_transitions   = actual.get_symbolic_transitions_for_state(actual_state);
        if (expected_transitions.size() != actual_transitions.size()) {
            TEST_PRINT_DEBUG("Found mismatch between symbolic transition sizes.\n");
            TEST_PRINT_DEBUG("  (expected): ");
            for (Transition& et: expected_transitions) TEST_PRINT_DEBUG(et << ", ");
            TEST_PRINT_DEBUG("\n");
            TEST_PRINT_DEBUG("  (actual): ");
            for (Transition& at: actual_transitions) TEST_PRINT_DEBUG(at << ", ");
            TEST_PRINT_DEBUG("\n");
            CHECK(expected_transitions.size() == actual_transitions.size());
        }

        for (auto& expected_transition: expected_transitions) {
            // Find a transition matching origin and the transition symbol; make sure that there is exactly 1 as the automaton is a DFA
            Transition& matching_transition = actual_transitions[0];
            u64 match_count = 0;
            for (auto& actual_transition: actual_transitions) {
                if (actual_transition.symbol == expected_transition.symbol) {
                    matching_transition = actual_transition;
                    match_count++;
                }
            }

            CHECK(match_count == 1);

            // In case there the already was an (expected, actual) mapping in the isomorphism check that it
            // is consistent with the currently matching transition
            auto [pos, was_inserted] = isomorphism.emplace(expected_transition.to, matching_transition.to);
            if (!was_inserted) {
                CHECK(matching_transition.to == pos->second);
            } else {
                work_queue.push_back(expected_transition.to);
            }
        }
    }

    for (auto expected_final_state: expected.final_states) {
        auto mapping_it = isomorphism.find(expected_final_state);
        CHECK(mapping_it != isomorphism.end());

        // if (mapping_it == isomorphism.end()) {
        //     std::cout << "Cloud not map expected final state " << expected_final_state << " to an actual one!" << std::endl;
        //     std::cout << "Mappings: " << std::endl;
        //     for (auto [key, value]: isomorphism) {
        //         std::cout << key << " -> " << value << std::endl;
        //     }

        //     std::cout << "Actual final states: " << std::endl;
        //     for (auto actual: actual.final_states) {
        //         std::cout << actual << std::endl;
        //     }
        // }

        CHECK(actual.final_states.find(mapping_it->second) != actual.final_states.end());
    }
}

TEST_CASE("lazy_construct `\\exists x (x + y <= 0)`")
{
    Atom_Allocator allocator;
    Formula formula = make_formula(allocator, {}, {}, {{1, 1}}, {0});;
    Formula_Pool pool = Formula_Pool(&formula);
    auto formula_id = pool.store_formula(formula);
    Conjunction_State init_state({0});

    sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
    vars = sylvan::mtbdd_set_add(vars, 1);
    vars = sylvan::mtbdd_set_add(vars, 2);

    auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

    NFA expected_nfa(vars, 2, {0, 1}, {1}, {0});

    vector<Transition> symbolic_transitions {
        {.from = 0, .to = 1, .symbol = {2, 2}},
        {.from = 1, .to = 1, .symbol = {2, 2}},
    };

    for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

    assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
}

TEST_CASE("lazy_construct simple atoms")
{
    SUBCASE("atom `x - y <= 2`") {
        Atom_Allocator alloc;
        Formula formula = make_formula(alloc, {}, {}, {{1, 1}}, {});
        Formula_Pool pool = Formula_Pool(&formula);
        auto formula_id = pool.store_formula(formula);
        Conjunction_State init_state({2});

        sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
        vars = sylvan::mtbdd_set_add(vars, 1);
        vars = sylvan::mtbdd_set_add(vars, 2);

        auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

        NFA expected_nfa(
            vars, 2,
            {
                0,  // {2}
                1,  // {1, F}
                2,  // {0, F}
                3,  // {-1, F}
                4,  // {-1}
                5,  // {-2, F}
                6,  // {-2}
            },
            {1, 2, 3, 5},
            {0}
        );

        vector<Transition> symbolic_transitions {
            // {2}
            {.from = 0, .to = 1, .symbol = {0, 0}},
            {.from = 0, .to = 2, .symbol = {0, 1}},
            {.from = 0, .to = 2, .symbol = {1, 2}},
            // {1, F}
            {.from = 1, .to = 2, .symbol = {0, 2}},
            {.from = 1, .to = 2, .symbol = {1, 0}},
            {.from = 1, .to = 3, .symbol = {1, 1}},
            // {0, F}
            {.from = 2, .to = 2, .symbol = {0, 0}},
            {.from = 2, .to = 3, .symbol = {0, 1}},
            {.from = 2, .to = 3, .symbol = {1, 2}},
            // {-1, F}
            {.from = 3, .to = 4, .symbol = {0, 0}},
            {.from = 3, .to = 3, .symbol = {0, 1}},
            {.from = 3, .to = 3, .symbol = {1, 0}},
            {.from = 3, .to = 5, .symbol = {1, 1}},
            // {-1}
            {.from = 4, .to = 4, .symbol = {0, 0}},
            {.from = 4, .to = 3, .symbol = {0, 1}},
            {.from = 4, .to = 3, .symbol = {1, 0}},
            {.from = 4, .to = 5, .symbol = {1, 1}},
            // {-2, F}
            {.from = 5, .to = 4, .symbol = {0, 0}},
            {.from = 5, .to = 6, .symbol = {0, 1}},
            {.from = 5, .to = 6, .symbol = {1, 0}},
            {.from = 5, .to = 5, .symbol = {1, 1}},
            // {-2}
            {.from = 6, .to = 4, .symbol = {0, 0}},
            {.from = 6, .to = 6, .symbol = {0, 1}},
            {.from = 6, .to = 6, .symbol = {1, 0}},
            {.from = 6, .to = 5, .symbol = {1, 1}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }

    SUBCASE("atom `2x - y = 0`") {
        Atom_Allocator alloc;
        Formula formula = make_formula(alloc, {}, {{2, -1}}, {}, {});
        Formula_Pool pool = Formula_Pool(&formula);
        auto formula_id = pool.store_formula(formula);
        Conjunction_State init_state({0});

        sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
        vars = sylvan::mtbdd_set_add(vars, 1);
        vars = sylvan::mtbdd_set_add(vars, 2);

        auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

        NFA expected_nfa(
            vars, 2,
            {
                0,  // {0}
                1,  // {0, F}
                2,  // {-1}
                3,  // {-1, F}
                4,  // Trap
            },
            {1, 3},
            {0}
        );

        vector<Transition> symbolic_transitions {
            // {0}
            {.from = 0, .to = 1, .symbol = {0, 0}},
            {.from = 0, .to = 4, .symbol = {0, 1}},
            {.from = 0, .to = 2, .symbol = {1, 0}},
            {.from = 0, .to = 4, .symbol = {1, 1}},
            // {0, F}
            {.from = 1, .to = 1, .symbol = {0, 0}},
            {.from = 1, .to = 4, .symbol = {0, 1}},
            {.from = 1, .to = 2, .symbol = {1, 0}},
            {.from = 1, .to = 4, .symbol = {1, 1}},
            // {-1}
            {.from = 2, .to = 4, .symbol = {0, 0}},
            {.from = 2, .to = 0, .symbol = {0, 1}},
            {.from = 2, .to = 4, .symbol = {1, 0}},
            {.from = 2, .to = 3, .symbol = {1, 1}},
            // {-1, F}
            {.from = 3, .to = 4, .symbol = {0, 0}},
            {.from = 3, .to = 0, .symbol = {0, 1}},
            {.from = 3, .to = 4, .symbol = {1, 0}},
            {.from = 3, .to = 3, .symbol = {1, 1}},
            // Trap
            {.from = 4, .to = 4, .symbol = {2, 2}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }

    SUBCASE("atom `x + 3y = 1 (mod 3)`") {
        Atom_Allocator alloc;
        Formula formula = make_formula(alloc, {{{1, 3}, 3}}, {}, {}, {});
        Formula_Pool pool = Formula_Pool(&formula);
        auto formula_id = pool.store_formula(formula);
        Conjunction_State init_state({1});

        sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
        vars = sylvan::mtbdd_set_add(vars, 1);
        vars = sylvan::mtbdd_set_add(vars, 2);

        auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

        NFA expected_nfa(
            vars, 2,
            {
                0,  // {1}
                1,  // {2}
                2,  // {2, F}
                3,  // {0}
                4,  // {0, F}
            },
            {2, 4},
            {0}
        );

        vector<Transition> symbolic_transitions {
            // {1}
            {.from = 0, .to = 1, .symbol = {0, 2}},
            {.from = 0, .to = 3, .symbol = {1, 2}},
            // {2}
            {.from = 1, .to = 0, .symbol = {0, 2}},
            {.from = 1, .to = 2, .symbol = {1, 2}},
            // {2, F}
            {.from = 2, .to = 0, .symbol = {0, 2}},
            {.from = 2, .to = 2, .symbol = {1, 2}},
            // {0}
            {.from = 3, .to = 4, .symbol = {0, 2}},
            {.from = 3, .to = 0, .symbol = {1, 2}},
            // {0, F}
            {.from = 4, .to = 4, .symbol = {0, 2}},
            {.from = 4, .to = 0, .symbol = {1, 2}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }
    SUBCASE("atom `x + y = 0 (mod 4)`") {
        Atom_Allocator alloc;
        Formula formula = make_formula(alloc, {{{1, 1}, 4}}, {}, {}, {});
        Formula_Pool pool = Formula_Pool(&formula);
        auto formula_id = pool.store_formula(formula);
        Conjunction_State init_state({0});

        sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
        vars = sylvan::mtbdd_set_add(vars, 1);
        vars = sylvan::mtbdd_set_add(vars, 2);

        auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

        NFA expected_nfa(
            vars, 2,
            {
                0,  // {0 (mod 4)}
                1,  // {0 (mod 2), TOP}
                2,  // {1 (mod 2)}
                3,  // {0 (mod 1), TOP}
                4,  // {BOTTOM}
            },
            {1, 3},
            {0}
        );

        vector<Transition> symbolic_transitions {
            // {0 (mod 4)}
            {.from = 0, .to = 1, .symbol = {0, 0}},
            {.from = 0, .to = 4, .symbol = {0, 1}},
            {.from = 0, .to = 4, .symbol = {1, 0}},
            {.from = 0, .to = 2, .symbol = {1, 1}},
            // {0 (mod 2), TOP}
            {.from = 1, .to = 3, .symbol = {0, 0}},
            {.from = 1, .to = 4, .symbol = {0, 1}},
            {.from = 1, .to = 4, .symbol = {1, 0}},
            {.from = 1, .to = 3, .symbol = {1, 1}},
            // {1 (mod 2), TOP}
            {.from = 2, .to = 4, .symbol = {0, 0}},
            {.from = 2, .to = 3, .symbol = {0, 1}},
            {.from = 2, .to = 3, .symbol = {1, 0}},
            {.from = 2, .to = 4, .symbol = {1, 1}},
            // {0 (mod 1), TOP}
            {.from = 3, .to = 3, .symbol = {2, 2}},
            // {BOTTOM}
            {.from = 4, .to = 4, .symbol = {2, 2}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }
}

/*
TEST_CASE("lazy_construct(1) `\\exists y,m (x - y <= -1 && y <= -1 && -m <= 0 && m <= 1 && m - y ~ 0 (mod 3))`")
{
    Formula formula = {
        .atoms = {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, -1}),        // (<= (+ x (- y)) -1)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 1}),         // (<= y -1)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0}),        // (<= (- m) 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0}),         // (<= m 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1}, 3),  // (= (mod (+ m (- y)) 3) 0)
        },
        .bound_vars = {0, 2},
        .var_count = 3
    };

    Formula_Pool pool = Formula_Pool();
    auto formula_id = pool.store_formula(formula);
    Conjunction_State init_state({-1, -1, 0, 1, 0});

    sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
    vars = sylvan::mtbdd_set_add(vars, 1);
    vars = sylvan::mtbdd_set_add(vars, 2);
    vars = sylvan::mtbdd_set_add(vars, 3);

    auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

    NFA expected_nfa(
        vars, 3,
        {
            0, // {(-1, -1, 0, 0, 0)}, formula: input
            1, // {(-2)}, formula: x <= ?
            2, // {(-1)}, formula: x <= ?
            3  // {(-1), F} formula: x <= ?
        },
        {3},
        {0}
    );

    // Symbols are (m, x, y)
    vector<Transition> symbolic_transitions {
        {.from = 0, .to = 1, .symbol = {2, 2, 2}},
        {.from = 1, .to = 2, .symbol = {2, 0, 2}},
        {.from = 1, .to = 1, .symbol = {2, 1, 2}},
        {.from = 2, .to = 2, .symbol = {2, 0, 2}},
        {.from = 2, .to = 3, .symbol = {2, 1, 2}},
        {.from = 3, .to = 2, .symbol = {2, 0, 2}},
        {.from = 3, .to = 3, .symbol = {2, 1, 2}},
    };

    for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

    assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
}
*/

TEST_CASE("complex formula :: lazy_construct `\\exists y,m (x - y <= 5 && z - m <= -1 && y <= -1 && -m <= 0 && m <= 1 && m - y ~ 0 (mod 3))`")
{
    // (exists ((y Int), (m Int))
    //   (land
    //     (<= (+ x (- y)) -1)
    //     (<= (- z m) -1)
    //     (<= y -1)
    //     (<= (- m) 0)
    //     (<= m 0)
    //     (= (mod (+ m (- y)) 299_993) 303)
    //   )
    // )
    //
    // Variables are renamed as: m -> x0, x -> x1, y -> x2, z -> x3
    //
    Atom_Allocator allocator;
    Formula formula = make_formula(allocator,
                                   {{{1, 0, -1, 0}, 3}},
                                   {},
                                   {
                                       {0, 1, -1, 0},   // (<= (- x y) 5)
                                       {-1, 0, 0, 1},   // (<= (- z m) -1)
                                       {0, 0, 1, 0},    // (<= y -1)
                                       {-1, 0, 0, 0},   // (<= (- m) 0)
                                       {1, 0, 0, 0},    // (<= m 0)
                                   },
                                   {0, 2});

    Conjunction_State init_state({0 /*congruence*/, 5, -1, -1, 0, 0});

    Formula_Pool pool = Formula_Pool(&formula);
    const Formula* canon_formula_ptr = pool.store_formula(formula);

    u32 var_ids[] = {1, 2, 3, 4};
    sylvan::BDDSET vars = sylvan::mtbdd_set_from_array(var_ids, 4);

    NFA actual_dfa = build_nfa_with_formula_entailement(canon_formula_ptr, init_state, vars, pool);

    // The formula should be simplified into: x <= 2 && z <= -1
    NFA expected_dfa(vars,
                     4,
                     {
                        0,  // (2, -1)
                        1,  // (1, -1) + ACC
                        2,  // (1, -1)
                        3,  // (0, -1)
                        4,  // (0, -1) + ACC
                        5,  // (-1, -1)
                        6   // (-1, -1) + ACC
                     },
                     {1, 4, 6}, {0});

    vector<Transition> expected_transitions {
        // (2, -1)
        {.from = 0, .to = 2, .symbol = {2, 0, 2, 0}},
        {.from = 0, .to = 1, .symbol = {2, 0, 2, 1}},
        {.from = 0, .to = 3, .symbol = {2, 1, 2, 0}},
        {.from = 0, .to = 4, .symbol = {2, 1, 2, 1}},
        // (1, -1) + ACC
        {.from = 1, .to = 3, .symbol = {2, 2, 2, 0}},
        {.from = 1, .to = 4, .symbol = {2, 2, 2, 1}},
        // (1, -1)
        {.from = 2, .to = 3, .symbol = {2, 2, 2, 0}},
        {.from = 2, .to = 4, .symbol = {2, 2, 2, 1}},
        // (0, -1)
        {.from = 3, .to = 3, .symbol = {2, 0, 2, 0}},
        {.from = 3, .to = 4, .symbol = {2, 0, 2, 1}},
        {.from = 3, .to = 5, .symbol = {2, 1, 2, 0}},
        {.from = 3, .to = 6, .symbol = {2, 1, 2, 1}},
        // (0, -1) + ACC
        {.from = 4, .to = 3, .symbol = {2, 0, 2, 0}},
        {.from = 4, .to = 4, .symbol = {2, 0, 2, 1}},
        {.from = 4, .to = 5, .symbol = {2, 1, 2, 0}},
        {.from = 4, .to = 6, .symbol = {2, 1, 2, 1}},
        // (-1, -1)
        {.from = 5, .to = 5, .symbol = {2, 0, 2, 2}},
        {.from = 5, .to = 5, .symbol = {2, 1, 2, 0}},
        {.from = 5, .to = 6, .symbol = {2, 1, 2, 1}},
        // (-1, -1) + ACC
        {.from = 6, .to = 5, .symbol = {2, 0, 2, 2}},
        {.from = 6, .to = 5, .symbol = {2, 1, 2, 0}},
        {.from = 6, .to = 6, .symbol = {2, 1, 2, 1}},
    };

    for (Transition& transition: expected_transitions)
        expected_dfa.add_transition(transition.from, transition.to, transition.symbol.data());

    assert_dfas_are_isomorphic(expected_dfa, actual_dfa);
}

TEST_CASE("NFA::pad_closure (simple)") {
    sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
    vars = sylvan::mtbdd_set_add(vars, 1);

    NFA nfa(vars, 1, {1, 2, 3}, {3}, {1});

    u8 symbol[] = {1};
    nfa.add_transition(1, 2, symbol);
    nfa.add_transition(2, 3, symbol);

    std::vector<Transition> expected_transitions = nfa_unpack_transitions(nfa);
    expected_transitions.push_back({.from=1, .to=4, .symbol={1}});

    nfa.perform_pad_closure();
    auto transitions = nfa_unpack_transitions(nfa);
    CHECK(transitions.size() == 3);
    for (auto& expected_transition: expected_transitions) {
        CHECK(std::find(transitions.begin(), transitions.end(), expected_transition) != transitions.end());
    }
}

TEST_CASE("NFA::determinize (simple)") {
    sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
    vars = sylvan::mtbdd_set_add(vars, 1);
    NFA nfa(vars, 1, {1, 2, 3}, {3}, {1});

    nfa.add_transition(1, 1, vector<u8>{0});
    nfa.add_transition(1, 2, vector<u8>{0});
    nfa.add_transition(2, 3, vector<u8>{0});

    auto actual_dfa = determinize_nfa(nfa);

    NFA expected_dfa(
        vars, 1,
        {
            1, // {1}
            2, // {1, 2}
            3, // {1, 2, 3}
            4, // trap
        },
        {3},
        {1}
    );
    expected_dfa.add_transition(1, 2, vector<u8>{0});
    expected_dfa.add_transition(1, 4, vector<u8>{1});

    expected_dfa.add_transition(2, 3, vector<u8>{0});
    expected_dfa.add_transition(2, 4, vector<u8>{1});

    expected_dfa.add_transition(3, 3, vector<u8>{0});
    expected_dfa.add_transition(3, 4, vector<u8>{1});

    expected_dfa.add_transition(4, 4, vector<u8>{2});

    assert_dfas_are_isomorphic(expected_dfa, actual_dfa);
}


TEST_CASE("NFA::intersection (simple)") {
    u32 var_ids[] = {1, 2};
    sylvan::BDDSET vars = sylvan::mtbdd_set_from_array(var_ids, 2);

    NFA left_nfa(vars), right_nfa(vars);
    {
        left_nfa.var_count = 2;
        left_nfa.states = {1, 2, 3, 4, 5};
        left_nfa.initial_states = {1};
        left_nfa.final_states = {5};

        left_nfa.add_transition(1, 2, {1, 1});
        left_nfa.add_transition(1, 3, {0, 0});
        left_nfa.add_transition(3, 4, {1, 1});
        left_nfa.add_transition(3, 5, {0, 0});
    }
    {
        right_nfa.var_count = 2;
        right_nfa.states = {1, 2, 3, 4, 5};
        right_nfa.initial_states = {1};
        right_nfa.final_states = {5};

        right_nfa.add_transition(1, 2, {0, 1}); // < differs in 1st bit from left
        right_nfa.add_transition(1, 3, {0, 0});
        right_nfa.add_transition(3, 4, {0, 1}); // < differs in 1st bit from left
        right_nfa.add_transition(3, 5, {0, 0});
    }

    auto result = compute_nfa_intersection(left_nfa, right_nfa);

    NFA expected_nfa(vars, 2u, {1, 2, 3}, {3}, {1});

    expected_nfa.add_transition(1, 2, {0, 0});
    expected_nfa.add_transition(2, 3, {0, 0});
    assert_dfas_are_isomorphic(expected_nfa, result);
}

TEST_CASE("remove_nonfinishing_states :: simple") {
    u32 var_arr[] = {1, 2};
    sylvan::BDDSET vars = sylvan::mtbdd_set_from_array(var_arr, 2);

    NFA dfa(vars, 2, {1, 2, 3, 4, 5, 6}, {3}, {1});
    dfa.add_transition(1, 2, {0, 0});
    dfa.add_transition(2, 3, {0, 0});

    // Transitions to states not reaching any of the final states
    dfa.add_transition(1, 4, {0, 1});
    dfa.add_transition(1, 5, {1, 0});
    dfa.add_transition(2, 5, {1, 1});
    dfa.add_transition(4, 5, {1, 1});

    remove_nonfinishing_states(dfa);

    NFA expected_dfa(vars, 2, {1, 2, 3}, {3}, {1});

    expected_dfa.add_transition(1, 2, {0, 0});
    expected_dfa.add_transition(2, 3, {0, 0});

    assert_dfas_are_isomorphic(expected_dfa, dfa);
}

TEST_CASE("Minimization - already minimal DFA ") {
    u32 var_arr[] = {1};
    sylvan::BDDSET vars = sylvan::mtbdd_set_from_array(var_arr, 1);
    NFA already_minimal_dfa(vars, 1, {1, 2, 3}, {2}, {1});

    already_minimal_dfa.add_transition(1, 1, (std::vector<u8>){0});
    already_minimal_dfa.add_transition(1, 2, (std::vector<u8>){1});
    already_minimal_dfa.add_transition(2, 2, (std::vector<u8>){0});
    already_minimal_dfa.add_transition(2, 3, (std::vector<u8>){1});
    already_minimal_dfa.add_transition(3, 3, (std::vector<u8>){2});

    NFA result = minimize_hopcroft(already_minimal_dfa);

    assert_dfas_are_isomorphic(already_minimal_dfa, result);
}

TEST_CASE("Minimization - Wiki automaton") {
    u32 var_arr[] = {1};
    sylvan::BDDSET vars = sylvan::mtbdd_set_from_array(var_arr, 1);

    NFA input(
        vars, 1,
        {
            1, // a
            2, // b
            3, // c
            4, // d
            5, // e
            6, // f
        },
        {3, 4, 5},
        {1}
    );

    input.add_transition(1, 2, (std::vector<u8>){0});
    input.add_transition(1, 3, (std::vector<u8>){1});

    input.add_transition(2, 1, (std::vector<u8>){0});
    input.add_transition(2, 4, (std::vector<u8>){1});

    input.add_transition(3, 5, (std::vector<u8>){0});
    input.add_transition(3, 6, (std::vector<u8>){1});

    input.add_transition(4, 5, (std::vector<u8>){0});
    input.add_transition(4, 6, (std::vector<u8>){1});

    input.add_transition(5, 5, (std::vector<u8>){0});
    input.add_transition(5, 6, (std::vector<u8>){1});

    input.add_transition(6, 6, (std::vector<u8>){2});

    NFA expected_result(vars, 1, {1, 2, 3}, {2}, {1});

    expected_result.add_transition(1, 1, (std::vector<u8>){0});
    expected_result.add_transition(1, 2, (std::vector<u8>){1});
    expected_result.add_transition(2, 2, (std::vector<u8>){0});
    expected_result.add_transition(2, 3, (std::vector<u8>){1});
    expected_result.add_transition(3, 3, (std::vector<u8>){2});

    NFA result = minimize_hopcroft(input);

    assert_dfas_are_isomorphic(expected_result, result);
}

TEST_CASE("Dep. analysis :: potential var identifiction") {
    // m = vars[0] > x(1) > y(2) > z(3)
    Atom_Allocator allocator;
    Formula formula = make_formula(
        allocator,
        {{{1, 0, -1, 0}, 10}}, // m - y ~ 0
        {},
        {
            {0, 1, -1, 0},  // x <= y
            {1, 0, 0, -1},  // m <= z
            {-1, 0, 0, 0},  // m >= 0
            {1, 0, 0, 0}    // m <= 0
        },
        {0, 2}
    );

    auto graph = build_dep_graph(formula);
    identify_potential_variables(graph);

    vector<u64> expected_potent_vars = {0};
    CHECK(graph.potential_vars == expected_potent_vars);
}

TEST_CASE("Dep. analysis :: simplify (const var)") {
    // m(0) > x(1) > y(2) > z(3)
    Atom_Allocator alloc;
    Formula formula = make_formula(
        alloc,
        {{{1, 0, -1, 0}, 10}},  // m - y ~ 0
        {},
        {
            {0, 1, -1, 0},  // x <= y
            {1, 0, 0, -1},  // m <= z
            {-1, 0, 0, 0},  // m >= 0
            {1, 0, 0,  0},  // m >= 0
        },
        {0, 2}
    );

    auto graph = build_dep_graph(formula);
    identify_potential_variables(graph);

    Conjunction_State state({0, 0, 0, 1, 1});
    bool was_simplified = simplify_graph(graph, state);

    CHECK(was_simplified);

    // The graph should be simplified to
    // 0. x <= y, m ~ y, m <= z, m >= 1, m <= 1
    // 1. x <= y, 1 ~ y, 1 <= z
    // 2. 1 <= z
    CHECK( graph.linear_nodes[0].is_satisfied);
    CHECK(!graph.linear_nodes[1].is_satisfied);
    CHECK( graph.linear_nodes[2].is_satisfied);
    CHECK( graph.linear_nodes[3].is_satisfied);
    CHECK( graph.congruence_nodes[0].is_satisfied);

    auto& actual_atom = graph.linear_nodes[1];
    Presburger_Atom expected_atom(PR_ATOM_INEQ, {0, 0, 0, -1}); // 0 <= z
    vector<s64> expected_coefs = {0, 0, 0, -1};
    CHECK(actual_atom.coefs == expected_coefs);
    CHECK(state.constants[2] == -1);
}

TEST_CASE("Dep. analysis :: simplify (unbound vars)") {
    // vars = {[0]=x1, x2 < x3 < x4}
    Atom_Allocator alloc;
    Formula formula = make_formula(
        alloc,
        {{{1, 0, 0, -1}, 13}},
        {},
        {
            {-1, 0, 0, 0},  // 0  <= x1
            {0, -1, 0, 0},  // 23 <= x2
            {-1, 5, 0, 0},  // 5*x2 - x1 <= 0
            {0, 0, -1, 1},  // x4 <= x3
            {0, 0, 0, -1},  // x4 >= 0
            {0, 0, 0,  1}   // x4 <= 12
        },
        {0, 1, 3}// x1 ~ x4
    );

    auto graph = build_dep_graph(formula);
    identify_potential_variables(graph);

    Conjunction_State state({0, -23, 0, 0, 0, 12, 0});

    auto was_simplified = simplify_graph(graph, state);
    CHECK(was_simplified);

    // 0)  0 <= x1, 23 <= x2, 5x2 <= x1, x4 <= x3, x4 >= 0, x4 <= 12, x1 ~ x4
    // 1)  23 <= x2, x4 <= x3, x4 >= 0, x4 <= 12
    // 2)  x4 <= x3, x4 >= 0, x4 <= 12
    // 3)  0 <= x3
    CHECK( graph.linear_nodes[0].is_satisfied);
    CHECK( graph.linear_nodes[1].is_satisfied);
    CHECK( graph.linear_nodes[2].is_satisfied);
    CHECK(!graph.linear_nodes[3].is_satisfied);
    CHECK( graph.linear_nodes[4].is_satisfied);
    CHECK( graph.linear_nodes[5].is_satisfied);
    CHECK(graph.congruence_nodes[0].is_satisfied);

    CHECK(state.constants[3] == 0);
}

TEST_CASE("Dep. analysis :: simplify (presentation formula)") {
    Atom_Allocator allocator;
    // m(0) < x < y < z
    const u64 modulus = 299993;
    Formula formula = make_formula(
        allocator,
        {{{1, 0, -1, 0}, modulus}},// m - y ~ 303
        {},
        {
            {0, 1, -1, 0},   // x - y <= -1
            {1, 0, 0, -1},   // m - z <= -1
            {0, 0, 1,  0},    // y <= -1
            {-1, 0, 0, 0},   // -m <= 0
            {1, 0, 0,  0},   // m <= 12
        },
        {0, 2}
    );

    auto graph = build_dep_graph(formula);
    identify_potential_variables(graph);

    Conjunction_State state({303, -1, -1, -1, 0, 0});

    auto was_simplified = simplify_graph(graph, state);
    CHECK(was_simplified);
    // 0)  x - y <= -1, m - z <= -1, y <= -1, -m <= 0, m <= M, m - y ~ 303
    // 0)  x - y <= -1, -z <= -1, y <= -1, -y ~ 303   (instantiate y=-303)
    // 0)  x <= -304, -z <= -1
    CHECK(!graph.linear_nodes[0].is_satisfied);
    CHECK(!graph.linear_nodes[1].is_satisfied);
    CHECK( graph.linear_nodes[2].is_satisfied);
    CHECK( graph.linear_nodes[3].is_satisfied);
    CHECK( graph.linear_nodes[4].is_satisfied);
    CHECK( graph.congruence_nodes[0].is_satisfied);

    CHECK(state.constants[1] == -304);
    CHECK(state.constants[2] == -1);

    Formula_Allocator formula_allocator (&formula);

    auto stateful_formula = convert_graph_into_formula(graph, formula_allocator, state);
    auto& simplified_formula = stateful_formula.formula;

    CHECK(simplified_formula.inequations.size == 2);
    CHECK(simplified_formula.equations.size   == 0);
    CHECK(simplified_formula.congruences.size == 0);
    CHECK(stateful_formula.state.constants.size() == 2);

    s64 simplfied_ineq_coefs1[] = {0, 1, 0,  0};
    s64 simplfied_ineq_coefs2[] = {0, 0, 0, -1};
    Sized_Array<s64> ineq1_coefs = {.items = simplfied_ineq_coefs1, .size = 4};
    Sized_Array<s64> ineq2_coefs = {.items = simplfied_ineq_coefs2, .size = 4};
    CHECK(simplified_formula.inequations.items[0].coefs == ineq1_coefs);
    CHECK(simplified_formula.inequations.items[1].coefs == ineq2_coefs);

    CHECK(simplified_formula.bound_vars.empty());

    vector<s64> expected_state = {-304, -1};
    CHECK(stateful_formula.state.constants == expected_state);
}

TEST_CASE("Test Structured_Macrostate insertions (pareto optimality)") {
    Atom_Allocator alloc;
    Formula formula1 = make_formula(alloc, {{{1, 2}, 3}}, {{1, 2}}, {{1, 2}}, {});

    Intermediate_Macrostate macrostate;

    Conjunction_State post1({1, 2, 3});
    insert_into_post_if_valuable2(macrostate, &formula1, post1);
    CHECK(macrostate.formulae[&formula1].size() == 1);

    insert_into_post_if_valuable2(macrostate, &formula1, post1);
    CHECK(macrostate.formulae[&formula1].size() == 1);

    Conjunction_State post2({1, 2, 4});
    insert_into_post_if_valuable2(macrostate, &formula1, post2);
    CHECK(macrostate.formulae[&formula1].size() == 1);

    insert_into_post_if_valuable2(macrostate, &formula1, post1);
    CHECK(macrostate.formulae[&formula1].size() == 1);


    Conjunction_State post3({1, 3, 4});
    insert_into_post_if_valuable2(macrostate, &formula1, post3);
    CHECK(macrostate.formulae[&formula1].size() == 2);
    CHECK(macrostate.formulae[&formula1][{1, 3}].size() == 1);
    CHECK(macrostate.formulae[&formula1][{1, 2}].size() == 1);


    Formula formula2 = make_formula(alloc, {}, {}, {{1, 2}, {1, 3}}, {});
    Conjunction_State post2_1({1, 1});
    Conjunction_State post2_2({1, 2});
    Conjunction_State post2_3({2, 1});
    Conjunction_State post2_4({2, 2});
    insert_into_post_if_valuable2(macrostate, &formula2, post2_1);
    insert_into_post_if_valuable2(macrostate, &formula2, post2_2);
    insert_into_post_if_valuable2(macrostate, &formula2, post2_3);
    insert_into_post_if_valuable2(macrostate, &formula2, post2_4);

    CHECK(macrostate.formulae[&formula2].size() == 1);
}

TEST_CASE("Test extended Euclidean") {
    CHECK(compute_multiplicative_inverse(900, 37) == 73);
    CHECK(compute_multiplicative_inverse(13, 12) == 12);
}

TEST_CASE("Test evaluate_mod_congruence_at_point") {
    SUBCASE("atom `y - m = 0 (mod 4)`") {
        Atom_Allocator alloc;
        Formula formula1 = make_formula(alloc, {{{1, -1}, 4}}, {}, {}, {});

        Formula_Allocator formula_allocator (&formula1);

        CHECK(formula1.dep_graph.congruence_nodes.size() == 1);
        Congruence_Node& node = formula1.dep_graph.congruence_nodes[0];

        Captured_Modulus mod = {.leading_var = 0, .subordinate_var = 1};
        s64 value = eval_mod_congruence_at_point(node, 0, mod, 0);
        CHECK(value == 0);

        value = eval_mod_congruence_at_point(node, 0, mod, 4);
        CHECK(value == 0);

        value = eval_mod_congruence_at_point(node, 0, mod, -1);
        CHECK(value == 3);
    }

    SUBCASE("atom `3*y - m = 0 (mod 5)`") {
        Atom_Allocator alloc;
        Formula formula1 = make_formula(alloc, {{{3, -1}, 5}}, {}, {}, {});

        Formula_Allocator formula_allocator (&formula1);

        CHECK(formula1.dep_graph.congruence_nodes.size() == 1);
        Congruence_Node& node = formula1.dep_graph.congruence_nodes[0];

        Captured_Modulus mod = {.leading_var = 0, .subordinate_var = 1};
        s64 value = eval_mod_congruence_at_point(node, 0, mod, 0);
        CHECK(value == 0);

        value = eval_mod_congruence_at_point(node, 0, mod, 4);
        CHECK(value == 2);

        value = eval_mod_congruence_at_point(node, 0, mod, -1);
        CHECK(value == 2);
    }
}

TEST_CASE("Test get_point_for_mod_congruence_to_obtain_value") {
    SUBCASE("atom `3*y - m = 0 (mod 5)`") {
        Atom_Allocator alloc;
        Formula formula1 = make_formula(alloc, {{{3, -1}, 5}}, {}, {}, {});

        Formula_Allocator formula_allocator (&formula1);

        CHECK(formula1.dep_graph.congruence_nodes.size() == 1);
        Congruence_Node& node = formula1.dep_graph.congruence_nodes[0];

        Captured_Modulus mod = {.leading_var = 0, .subordinate_var = 1};
        s64 value = get_point_for_mod_congruence_to_obtain_value(node, 0, mod, 0);
        CHECK(value == 0);

        value = get_point_for_mod_congruence_to_obtain_value(node, 0, mod, 3);
        CHECK(value == 1);

        value = get_point_for_mod_congruence_to_obtain_value(node, 0, mod, 4);
        CHECK(value == 3);
    }
    SUBCASE("atom `y - m = 0 (mod 5)` (real formula)") {
        Atom_Allocator alloc;
        Formula formula1 = make_formula(alloc, {{{1, -1}, 5}}, {}, {}, {});
        Formula_Allocator formula_allocator (&formula1);

        CHECK(formula1.dep_graph.congruence_nodes.size() == 1);
        Congruence_Node& node = formula1.dep_graph.congruence_nodes[0];

        Captured_Modulus mod = {.leading_var = 0, .subordinate_var = 1};
        s64 value = get_point_for_mod_congruence_to_obtain_value(node, 0, mod, 1);
        CHECK(value == 1);

        value = get_point_for_mod_congruence_to_obtain_value(node, 0, mod, 4);
        CHECK(value == 4);
    }
}

TEST_CASE("Test linearize_congruence - `y - m ~ 0 && 1 <= m <= 4` (real formula)") {
    Atom_Allocator alloc;
    Formula formula1 = make_formula(alloc, {{{1, -1}, 5}}, {}, {{0, -1}, {0, 1}}, {});
    Conjunction_State state ({0, 1, 4});

    Captured_Modulus captured_mod = {.leading_var = 0, .subordinate_var = 1};
    Var_Preference preference = {.type = Var_Preference_Type::C_INCREASING, .c = -1};

    Linear_Function fn = linearize_congruence(formula1.dep_graph, state, 0, captured_mod, preference);
    CHECK(fn.a == 1);
    CHECK(fn.b == 5);
    CHECK(fn.valid);
}

TEST_CASE("Test linearize_formula - `\\exists y (x - y <= 0 && m - z <= 60_000 && y <= -1 && y - m ~ 0 (mod 299_993) && 1 <= m <= 299_992)` (real formula)") {
    Atom_Allocator alloc;
    // Vars: m -> x0, x -> x1, y -> x2, z -> x3
    Formula formula = make_formula(alloc,
                                   {{{-1, 0, 1, 0}, 299993}},  // y - m ~ 0 (mod 299_993)
                                   {},                   // No equations
                                   {
                                    {0, 1, -1, 0},       // x - y <= 0
                                    {1, 0, 0, -1},       // m - z <= 60_000
                                    {0, 0, 1, 0},        // y <= -1
                                    {-1, 0, 0, 0},       // m >= 1
                                    {1, 0, 0, 0},        // m <= 299_992
                                   },
                                   {});
    Formula_Allocator allocator(&formula);
    Conjunction_State state ({0 /*congruence*/, 0, 60000, -1, 1, 299992});

    auto result = linearize_formula(allocator, &formula, state);
    CHECK(result.has_value());
    auto linearized_formula = result.value().formula;
    CHECK(linearized_formula.congruences.size == 0);
    CHECK(linearized_formula.equations.size == 1);
    CHECK(linearized_formula.inequations.size == formula.inequations.size);

    auto& equation = linearized_formula.equations.items[0];

    s64 expected_coef_data[] = {-1, 0, 1, 0};
    Sized_Array<s64> expected_coefs = {.items = expected_coef_data, .size = 4};
    CHECK(equation.coefs == expected_coefs);

    auto& modified_state = result.value().state;
    std::cout << "modified state: " << modified_state << "\n";
    vector<s64> expected_state_data = {-299993 /*equation*/, 0, 60000, -1, 1, 299992};
    CHECK(modified_state.constants == expected_state_data);
}

TEST_CASE("Test shift_interval") {
    SUBCASE("[0, 3] <= -1 with (mod 4)") {
        Interval interval = {.low = 0, .high = 3};
        Var_Preference preference = {.type = Var_Preference_Type::C_INCREASING, .c = -1};
        shift_interval(interval, preference, 4);

        CHECK(interval.high == -1);
        CHECK(interval.low  == -4);
    }

    SUBCASE("[0, 3] <= 10 with (mod 4)") {
        Interval interval = {.low = 0, .high = 3};
        Var_Preference preference = {.type = Var_Preference_Type::C_INCREASING, .c = -1};
        shift_interval(interval, preference, 4);

        CHECK(interval.high == -1);
        CHECK(interval.low  == -4);
    }
}

int main(int argc, char* argv[]) {
    init_mtbdd_libs();

    doctest::Context context(argc, argv);

    int test_overall_rc = context.run();

    if(context.shouldExit()) {
        return test_overall_rc;
    }

    return 0;
}
