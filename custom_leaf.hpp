#ifndef AMAYA_CUSTOM_LEAF_H
#define AMAYA_CUSTOM_LEAF_H

#include <inttypes.h>
#include <stdlib.h>
#include <sylvan.h>
#include "base.hpp"

#ifndef rotl64
static inline uint64_t rotl64(uint64_t x, int8_t r)
{
    return ((x << r) | (x >> (64 - r)));
}
#endif

sylvan::MTBDD
make_set_leaf(Transition_Destination_Set* value);

void
mk_set_leaf(uint64_t *target_destination_set_ptr);

void
destroy_set_leaf(uint64_t leaf_value);

int set_leaf_equals(uint64_t a_ptr, uint64_t b_ptr);

uint64_t
set_leaf_hash(const uint64_t contents_ptr, const uint64_t seed);

char *
set_leaf_to_str(int comp, uint64_t leaf_val, char *buf, size_t buflen);

extern "C" {
	static uint64_t mtbdd_leaf_type_set;
}

#endif
