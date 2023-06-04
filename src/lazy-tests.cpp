#define DOCTEST_CONFIG_IMPLEMENT
#include "../external/doctest.h"
#include "../include/lazy.hpp"
#include "../include/base.hpp"

#include <algorithm>
#include <vector>
#include <unordered_map>

typedef Quantified_Atom_Conjunction Formula;

using std::vector;
using std::unordered_map;

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

        auto expected_transitions = expected.get_symbolic_transitions_for_state(expected_state);
        auto actual_transitions   = actual.get_symbolic_transitions_for_state(actual_state);

        CHECK(expected_transitions.size() == actual_transitions.size());

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

        /*
        if (mapping_it == isomorphism.end()) {
            std::cout << "Cloud not map expected final state " << expected_final_state << " to an actual one!" << std::endl;
            std::cout << "Mappings: " << std::endl;
            for (auto [key, value]: isomorphism) {
                std::cout << key << " -> " << value << std::endl;
            }

            std::cout << "Actual final states: " << std::endl;
            for (auto actual: actual.final_states) {
                std::cout << actual << std::endl;
            }
        }
        */

        CHECK(actual.final_states.find(mapping_it->second) != actual.final_states.end());
    }
}


TEST_CASE("lazy_construct `\\exists x (x + y <= 0)`")
{
    Formula formula = { .atoms = { Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 1})}, .bound_vars = {0}, .var_count = 2};
    Formula_Pool pool = Formula_Pool();
    auto formula_id = pool.store_formula(formula);
    Conjunction_State init_state({0});

    sylvan::BDDSET vars = sylvan::mtbdd_set_empty();
    vars = sylvan::mtbdd_set_add(vars, 1);
    vars = sylvan::mtbdd_set_add(vars, 2);

    auto actual_nfa = build_nfa_with_formula_entailement(formula_id, init_state, vars, pool);

    NFA expected_nfa(vars, 2, {0}, {0}, {0});

    vector<Transition> symbolic_transitions {
        {.from = 0, .to = 0, .symbol = {2, 2}},
    };

    for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

    assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
}

TEST_CASE("lazy_construct simple atoms")
{
    SUBCASE("atom `x - y <= 2`") {
        Formula formula = { .atoms = { Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 1})}, .bound_vars = {}, .var_count = 2};
        Formula_Pool pool = Formula_Pool();
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
            /* {2} */
            {.from = 0, .to = 1, .symbol = {0, 0}},
            {.from = 0, .to = 2, .symbol = {0, 1}},
            {.from = 0, .to = 2, .symbol = {1, 2}},
            /* {1, F} */
            {.from = 1, .to = 2, .symbol = {0, 2}},
            {.from = 1, .to = 2, .symbol = {1, 0}},
            {.from = 1, .to = 3, .symbol = {1, 1}},
            /* {0, F} */
            {.from = 2, .to = 2, .symbol = {0, 0}},
            {.from = 2, .to = 3, .symbol = {0, 1}},
            {.from = 2, .to = 3, .symbol = {1, 2}},
            /* {-1, F} */
            {.from = 3, .to = 4, .symbol = {0, 0}},
            {.from = 3, .to = 3, .symbol = {0, 1}},
            {.from = 3, .to = 3, .symbol = {1, 0}},
            {.from = 3, .to = 5, .symbol = {1, 1}},
            /* {-1} */
            {.from = 4, .to = 4, .symbol = {0, 0}},
            {.from = 4, .to = 3, .symbol = {0, 1}},
            {.from = 4, .to = 3, .symbol = {1, 0}},
            {.from = 4, .to = 5, .symbol = {1, 1}},
            /* {-2, F} */
            {.from = 5, .to = 4, .symbol = {0, 0}},
            {.from = 5, .to = 6, .symbol = {0, 1}},
            {.from = 5, .to = 6, .symbol = {1, 0}},
            {.from = 5, .to = 5, .symbol = {1, 1}},
            /* {-2} */
            {.from = 6, .to = 4, .symbol = {0, 0}},
            {.from = 6, .to = 6, .symbol = {0, 1}},
            {.from = 6, .to = 6, .symbol = {1, 0}},
            {.from = 6, .to = 5, .symbol = {1, 1}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }

    SUBCASE("atom `2x - y = 0`") {
        Formula formula = { .atoms = { Presburger_Atom(Presburger_Atom_Type::PR_ATOM_EQ, {2, -1})}, .bound_vars = {}, .var_count = 2};
        Formula_Pool pool = Formula_Pool();
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
            /* {0} */
            {.from = 0, .to = 1, .symbol = {0, 0}},
            {.from = 0, .to = 4, .symbol = {0, 1}},
            {.from = 0, .to = 2, .symbol = {1, 0}},
            {.from = 0, .to = 4, .symbol = {1, 1}},
            /* {0, F} */
            {.from = 1, .to = 1, .symbol = {0, 0}},
            {.from = 1, .to = 4, .symbol = {0, 1}},
            {.from = 1, .to = 2, .symbol = {1, 0}},
            {.from = 1, .to = 4, .symbol = {1, 1}},
            /* {-1} */
            {.from = 2, .to = 4, .symbol = {0, 0}},
            {.from = 2, .to = 0, .symbol = {0, 1}},
            {.from = 2, .to = 4, .symbol = {1, 0}},
            {.from = 2, .to = 3, .symbol = {1, 1}},
            /* {-1, F} */
            {.from = 3, .to = 4, .symbol = {0, 0}},
            {.from = 3, .to = 0, .symbol = {0, 1}},
            {.from = 3, .to = 4, .symbol = {1, 0}},
            {.from = 3, .to = 3, .symbol = {1, 1}},
            /* Trap */
            {.from = 4, .to = 4, .symbol = {2, 2}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }

    SUBCASE("atom `x + 3y = 1 (mod 3)`") {
        Formula formula = { .atoms = { Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 3}, 3)}, .bound_vars = {}, .var_count = 2};
        Formula_Pool pool = Formula_Pool();
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
            /* {1} */
            {.from = 0, .to = 1, .symbol = {0, 2}},
            {.from = 0, .to = 3, .symbol = {1, 2}},
            /* {2} */
            {.from = 1, .to = 0, .symbol = {0, 2}},
            {.from = 1, .to = 2, .symbol = {1, 2}},
            /* {2, F} */
            {.from = 2, .to = 0, .symbol = {0, 2}},
            {.from = 2, .to = 2, .symbol = {1, 2}},
            /* {0} */
            {.from = 3, .to = 4, .symbol = {0, 2}},
            {.from = 3, .to = 0, .symbol = {1, 2}},
            /* {0, F} */
            {.from = 4, .to = 4, .symbol = {0, 2}},
            {.from = 4, .to = 0, .symbol = {1, 2}},
        };

        for (auto it: symbolic_transitions) expected_nfa.add_transition(it.from, it.to, it.symbol.data());

        assert_dfas_are_isomorphic(expected_nfa, actual_nfa);
    }
}

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

TEST_CASE("lazy_construct `\\exists y,m (x - y <= -1 && && m - z <= -1 && y <= -1 && -m <= 0 && m <= 1 && m - y ~ 0 (mod 299993))`")
{
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
    Quantified_Atom_Conjunction real_formula(
        {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),                // (<= (- m) 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, 0}),                 // (<= m 299_992)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1, 0}, 299993),  // (= (mod (+ m (- y)) 299_993) 303)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, -1, 0, 0}),                // (<= (- x) 23)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, -1, 0}),                // (<= (- y) 0)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 5, -1, 0}),                // (<= (+ (- y) (* 5 x) 10)
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, -1}),                // (<= (+ m (- z))  7)
        },
        {0, 2},
        4
    );

    Conjunction_State real_state({-1, -1, -1, 0, 0, 303});

    Formula_Pool pool = Formula_Pool();
    pool.store_formula(real_formula);

    // auto entailment_status = simplify_stateful_formula(&real_formula, real_state, pool);
    //build_nfa_with_formula_entailement(pool, real_state);
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
    // m(3) > x(2) > y(1) > z(0)
    Quantified_Atom_Conjunction real_formula(
        {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, -1, 0}),             // x <= y
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1, 0}, 10),   // m - y ~ 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, -1}),             // m <= y
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),             // m >= 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, 0}),              // m >=1 0
        },
        {0, 2},
        4
    );
    auto graph = build_dep_graph(real_formula);
    identify_potential_variables(graph);

    vector<u64> expected_potent_vars = {0};
    CHECK(graph.potential_vars == expected_potent_vars);
}

TEST_CASE("Dep. analysis :: simplify (const var)") {
    // m(0) > x(1) > y(2) > z(3)
    Quantified_Atom_Conjunction real_formula(
        {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, -1, 0}),             // x <= y
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1, 0}, 10),   // m - y ~ 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, -1}),             // m <= y
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),             // m >= 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, 0}),              // m >= 0
        },
        {0, 2},
        4
    );
    auto graph = build_dep_graph(real_formula);
    identify_potential_variables(graph);

    Conjunction_State state({0, 0, 0, 1, 1});
    auto was_simplfied = simplify_graph(graph, state);

    CHECK(was_simplfied);

    // The graph should be simplified to
    // 0. x <= y, m ~ y, m <= z, m >= 1, m <= 1
    // 1. x <= y, 1 ~ y, 1 <= z
    // 2. 1 <= z
    CHECK(graph.atom_nodes[0].is_satisfied);
    CHECK(graph.atom_nodes[1].is_satisfied);
    CHECK(!graph.atom_nodes[2].is_satisfied);
    CHECK(graph.atom_nodes[3].is_satisfied);
    CHECK(graph.atom_nodes[4].is_satisfied);

    auto& actual_atom = graph.atom_nodes[2].atom;
    Presburger_Atom expected_atom(PR_ATOM_INEQ, {0, 0, 0, -1}); // 0 <= z
    CHECK(actual_atom == expected_atom);
    CHECK(state.constants[2] == -1);
}

TEST_CASE("Dep. analysis :: simplify (unbound vars)") {
    // x1(0) < x2 < x3 < x4
    Quantified_Atom_Conjunction conj(
        {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),   // 0  <= x1
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, -1, 0, 0}),   // 23 <= x2
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 5, 0, 0}),   // 5*x2 - x1 <= 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, -1, 1}),   // x4 <= x3
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 0, -1}),   // x4 >= 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 0,  1}),   // x4 <= 12
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, 0, -1}, 13),   // x1 ~ x4
        },
        {0, 1, 3},
        4
    );
    auto graph = build_dep_graph(conj);
    identify_potential_variables(graph);

    Conjunction_State state({0, -23, 0, 0, 0, 12, 0});

    auto was_simplified = simplify_graph(graph, state);
    CHECK(was_simplified);

    // 0)  0 <= x1, 23 <= x2, 5x2 <= x1, x4 <= x3, x4 >= 0, x4 <= 12, x1 ~ x4
    // 1)  23 <= x2, x4 <= x3, x4 >= 0, x4 <= 12
    // 2)  x4 <= x3, x4 >= 0, x4 <= 12
    // 3)  0 <= x3
    CHECK(graph.atom_nodes[0].is_satisfied);
    CHECK(graph.atom_nodes[1].is_satisfied);
    CHECK(graph.atom_nodes[2].is_satisfied);
    CHECK(!graph.atom_nodes[3].is_satisfied);
    CHECK(graph.atom_nodes[4].is_satisfied);
    CHECK(graph.atom_nodes[5].is_satisfied);
    CHECK(graph.atom_nodes[6].is_satisfied);

    CHECK(state.constants[3] == 0);
}

TEST_CASE("Dep. analysis :: simplify (presentation formula)") {
    // m(0) < x < y < z
    const u64 modulus = 299993;
    Quantified_Atom_Conjunction conj(
        {
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, -1, 0}),   // x - y <= -1
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, -1}),   // m - z <= -1
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 1, 0}),    // y <= -1
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {-1, 0, 0, 0}),   // -m <= 0
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_INEQ, {1, 0, 0, 0}),   // m <= 12
            Presburger_Atom(Presburger_Atom_Type::PR_ATOM_CONGRUENCE, {1, 0, -1, 0}, modulus),   // m - y ~ 303
        },
        {0, 2},
        4
    );
    auto graph = build_dep_graph(conj);
    identify_potential_variables(graph);

    Conjunction_State state({-1, -1, -1, 0, 0, 303});

    auto was_simplified = simplify_graph(graph, state);
    CHECK(was_simplified);
    // 0)  x - y <= -1, m - z <= -1, y <= -1, -m <= 0, m <= M, m - y ~ 303
    // 0)  x - y <= -1, -z <= -1, y <= -1, -y ~ 303   (instantiate y=-303)
    // 0)  x <= -304, -z <= -1
    CHECK(!graph.atom_nodes[0].is_satisfied);
    CHECK(!graph.atom_nodes[1].is_satisfied);
    CHECK(graph.atom_nodes[2].is_satisfied);
    CHECK(graph.atom_nodes[3].is_satisfied);
    CHECK(graph.atom_nodes[4].is_satisfied);
    CHECK(graph.atom_nodes[5].is_satisfied);

    std::cout << state.constants << std::endl;
    CHECK(state.constants[0] == -304);
    CHECK(state.constants[1] == -1);

    auto stateful_formula = convert_graph_into_formula(graph, state);
    auto& formula = stateful_formula.formula;

    CHECK(formula.atoms.size() == 2);
    CHECK(stateful_formula.state.constants.size() == 2);

    Presburger_Atom expected_atom0(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 1, 0, 0});
    Presburger_Atom expected_atom1(Presburger_Atom_Type::PR_ATOM_INEQ, {0, 0, 0, -1});
    CHECK(formula.atoms[0] == expected_atom0);
    CHECK(formula.atoms[1] == expected_atom1);
    CHECK(formula.bound_vars.empty());

    auto& new_state = stateful_formula.state;
    CHECK(new_state.constants[0] == -304);
    CHECK(new_state.constants[1] == -1);
}

TEST_CASE("Test extended Euclidean") {
    CHECK(compute_multiplicative_inverse(900, 37) == 73);
    CHECK(compute_multiplicative_inverse(13, 12) == 12);
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
