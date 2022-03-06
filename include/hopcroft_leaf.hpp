#ifndef AMAYA_HOPCROFT_LEAF_H
#define AMAYA_HOPCROFT_LEAF_H

#include "../include/base.hpp"

#include <map>
#include <set>
#include <string>
#include <sylvan.h>

#define HOPCROFT_UNION_OP_ID 0x6000000

extern uint64_t mtbdd_leaf_type_hopcroft;

class Hopcroft_Leaf_Contents {
public:
    std::map<State, std::set<State>> destination_to_origin_states;

    Hopcroft_Leaf_Contents();
    Hopcroft_Leaf_Contents(Hopcroft_Leaf_Contents& original);
    ~Hopcroft_Leaf_Contents();

    std::string to_str();
};


// Callbacks for Sylvan to support this leaf type
void     hopcroft_leaf_create(uint64_t* hopcroft_leaf_contents_ptr_param);
uint64_t hopcroft_leaf_hash(const uint64_t contents_ptr, const uint64_t seed);
char*    hopcroft_leaf_to_str(int is_leaf_complemented, uint64_t leaf_contents_untyped_ptr, char* buffer, size_t buffer_size);
void     hopcroft_leaf_destroy(uint64_t hl_untyped_ptr);
int      hopcroft_leaf_equals(uint64_t hl_a_untyped_ptr, uint64_t hl_b_untyped_ptr);


TASK_DECL_2(sylvan::MTBDD, hopcroft_leaf_from_normal, sylvan::MTBDD, uint64_t);
TASK_DECL_3(sylvan::MTBDD, hopcroft_leaf_union, sylvan::MTBDD*, sylvan::MTBDD*, uint64_t);

#endif
