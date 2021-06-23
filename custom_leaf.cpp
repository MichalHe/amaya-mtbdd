#include "custom_leaf.hpp"
#include "base.hpp"

#include <inttypes.h>

#include <string>
#include <cstring>
#include <utility>
#include <stdlib.h>
#include <sstream>
#include <sylvan.h>

uint64_t mtbdd_leaf_type_set;

sylvan::MTBDD
make_set_leaf(Transition_Destination_Set* value) {
	sylvan::MTBDD leaf = sylvan::mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) value);
    return leaf;
}

/**
 * Function is called when a new node is being inserted to the hash table.
 *
 * @param state_set_ptr 	Pointer to an old valid leaf value, which is being
 * copied to the hash table. This means that its a pointer to a poiter to a
 * vector.
 */
void mk_set_leaf(uint64_t *target_destination_set_ptr)
{
    Transition_Destination_Set *original_tds = (Transition_Destination_Set *)*target_destination_set_ptr;
    Transition_Destination_Set *new_tds = new Transition_Destination_Set(*original_tds); // Copy construct

    // Accorting to the GMP leaf implementation, it should suffice to just write
    // the new value to the state_set_ptr.
    *(Transition_Destination_Set **)target_destination_set_ptr = new_tds;
}

void destroy_set_leaf(uint64_t leaf_value)
{
    Transition_Destination_Set *tds = (Transition_Destination_Set *)leaf_value;
    delete tds;
}

int set_leaf_equals(uint64_t a_ptr, uint64_t b_ptr)
{
    Transition_Destination_Set *a_tds = (Transition_Destination_Set *)a_ptr;
    Transition_Destination_Set *b_tds = (Transition_Destination_Set *)b_ptr;

	//if (a_tds->automaton_id == b_tds->automaton_id)
	if ((*a_tds->destination_set) == (*b_tds->destination_set))
	{
		return true;
	};
	return false;
}

uint64_t set_leaf_hash(const uint64_t contents_ptr, const uint64_t seed)
{
    Transition_Destination_Set *tds = (Transition_Destination_Set *)contents_ptr;

	unsigned long hash = seed;

	for (auto state : *tds->destination_set)
		hash = state + (hash << 6) + (hash << 16) - hash;

	return hash;

}

char *set_leaf_to_str(int comp, uint64_t leaf_val, char *buf, size_t buflen)
{
    (void)comp;
    std::stringstream ss;
    auto tds = (Transition_Destination_Set *)leaf_val;
    ss << "{";
    uint32_t cnt = 1;
    for (auto i : *tds->destination_set)
    {
        ss << i;
        if (cnt < tds->destination_set->size())
        {
            ss << ",";
        }
        cnt++;
    }
    ss << "}";

    const std::string str(ss.str());

    // Does the resulting string fits into the provided buffer?
    const size_t required_buf_size = str.size() + 1; // With nullbyte
    if (required_buf_size <= buflen)
    {
        const char *cstr = str.c_str();
        std::memcpy(buf, cstr, sizeof(char) * required_buf_size);
        return buf;
    }
    else
    {
        char *new_buf = (char *)malloc(sizeof(char) * required_buf_size);
        std::memcpy(new_buf, str.c_str(), sizeof(char) * required_buf_size);
        return new_buf;
    }
}
