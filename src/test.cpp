#include "../include/hopcroft_leaf.hpp"
#include "../include/custom_leaf.hpp"
#include "../include/base.hpp"
#include "../include/wrapper.hpp"

#include <sylvan.h>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>

#include <iostream>
#include <set>
#include <vector>
#include <utility>
#include <assert.h>


using namespace sylvan;

MTBDD mk_normal_mtbdd(uint8_t cube[2], std::vector<State>&& dest_states) {
    sylvan::BDDSET variable_set = sylvan::mtbdd_set_empty();
    variable_set = sylvan::mtbdd_set_add(variable_set, 1);
    variable_set = sylvan::mtbdd_set_add(variable_set, 2);

    Transition_Destination_Set tds;
    tds.destination_set = new std::set<State>(dest_states.begin(), dest_states.end());

    sylvan::MTBDD leaf = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) &tds);

    sylvan::MTBDD mtbdd = sylvan::mtbdd_cube(variable_set, cube, leaf);

    return mtbdd;
}


int main() 
{
    init_machinery(); 
    LACE_ME;

    uint8_t cube[]{0, 1}; 
    auto normal_mtbdd0 = mk_normal_mtbdd(cube, std::vector<State>{1, 2, 3});
    
    auto normal_mtbdd1 = mk_normal_mtbdd(cube, std::vector<State>{3, 4, 5});

    sylvan::MTBDD hl_0 = mtbdd_uapply(normal_mtbdd0, TASK(hopcroft_leaf_from_normal), 100);
    sylvan::MTBDD hl_1 = mtbdd_uapply(normal_mtbdd1, TASK(hopcroft_leaf_from_normal), 101);

    MTBDD result = mtbdd_applyp(hl_0, hl_1, 0, TASK(hopcroft_leaf_union), HOPCROFT_UNION_OP_ID);

    mtbdd_fprintdot(stdout, result);
}
