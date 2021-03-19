#include <iterator>
#include <sylvan.h>
#include <sylvan_bdd.h>
#include <sylvan_common.h>
#include <sylvan_mt.h>
#include <sylvan_mtbdd.h>
#include <sylvan_stats.h>
#include <sylvan_gmp.h>
#include <assert.h>

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace sylvan;
using std::cout;
using std::endl;
using std::set;;
using std::vector;

static uint32_t mtbdd_leaf_type_set;

void print_states_set(set<int> *states) {
	int cnt = 1;
	cout << "{";
	for (auto state : *states) 
	{
		cout << state;
		if (cnt < states->size())
		{
			cout << ", ";
		}
		++cnt;
	}
	cout << "}" << endl;
}

#ifndef rotl64
static inline uint64_t
rotl64(uint64_t x, int8_t r)
{
    return ((x<<r) | (x>>(64-r)));
}
#endif

MTBDD 
make_set_leaf(std::set<int> *value)
{
	MTBDD leaf = mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) value);
	return leaf;
}

/**
 * Function is called when a new node is being inserted to the hash table
 *
 * @param state_set_ptr 	Pointer to an old valid leaf value, which is being copied to the hash table.
 * 							This means that its a pointer to a poiter to a vector.(which is our value).
 */
static void
mk_set_leaf(uint64_t *state_set_ptr)
{
	set<int>* old_states = (set<int> *) *state_set_ptr; // After dereferece we get the actual value - pointer to a set.
	set<int>* states_copy = new set<int>(*old_states); // Copy construct it

	// Accorting to the GMP leaf implementation, it should suffice to just write the new value to the state_set_ptr.
	*(set<int> **)state_set_ptr = states_copy;
}

static void
destroy_set_leaf(uint64_t leaf_value) 
{
	set<int>* states_set_ptr = (set<int> *) leaf_value;
	cout << "Delete called on: ";
	print_states_set(states_set_ptr);
	delete states_set_ptr;
}

int set_leaf_equals(uint64_t a_ptr, uint64_t b_ptr) 
{
	set<int> *states_a = (set<int> *) a_ptr;
	set<int> *states_b = (set<int> *) a_ptr;

	return (*states_a) == (*states_b);
}

static uint64_t
set_leaf_hash(const uint64_t contents_ptr, const uint64_t seed)
{
	cout << "Hashing!" << endl;
	set<int> *states = (set<int>*) contents_ptr;
    const uint64_t prime = 1099511628211;
    uint64_t hash = seed;

    for (auto i : *states) {
        hash = hash ^ i;
        hash = rotl64(hash, 47);
        hash = hash * prime;
    }

	return hash;
}



TASK_DECL_2(MTBDD, set_union, MTBDD*, MTBDD*);
TASK_IMPL_2(MTBDD, set_union, MTBDD*, pa, MTBDD*, pb) 
{
	MTBDD a = *pa, b = *pb;

	// When one leaf is empty set (false), 
	// the algorithm should return all reachable states for the other one
	if (a == mtbdd_false) 
	{
		return b; // If b is also mtbdd_false, nothing wrong is happening
	}
	if (b == mtbdd_false)
	{
		return a; // Same here
	}

	// If both are leaves, we calculate union
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
		std::set<int> set_a = *((std::set<int> *) mtbdd_getvalue(a));	
		std::set<int> set_b = *((std::set<int> *) mtbdd_getvalue(b));	

		std::set<int> *_union = new std::set<int>();
		std::set_union(
				set_a.begin(), set_a.end(), 
				set_b.begin(), set_b.end(),
				std::inserter(*_union, _union->begin())
				);

		
		print_states_set(_union);
		
		MTBDD union_leaf = make_set_leaf(_union); 

		return union_leaf;
    }

	// TODO: Perform pointer swap here, so the cache would be utilized (commutative).
	// TODO: Utilize operation cache.
    return mtbdd_invalid;
}

VOID_TASK_0(simple_fn)
{
	// Define custom leaf type
	mtbdd_leaf_type_set = sylvan_mt_create_type();

	sylvan_mt_set_hash(mtbdd_leaf_type_set , set_leaf_hash);
	sylvan_mt_set_equals(mtbdd_leaf_type_set, set_leaf_equals);
	sylvan_mt_set_create(mtbdd_leaf_type_set, mk_set_leaf);
	sylvan_mt_set_destroy(mtbdd_leaf_type_set, destroy_set_leaf);

    // sylvan_mt_set_to_str(gmp_type, gmp_to_str);
    // sylvan_mt_set_write_binary(gmp_type, gmp_write_binary);
    // sylvan_mt_set_read_binary(gmp_type, gmp_read_binary);

	BDD b1 = sylvan_ithvar(1);
	BDD b2 = sylvan_ithvar(2);

	std::set<int> set_A;
	set_A.insert(1);
	MTBDD leaf_A = make_set_leaf(&set_A);

	std::set<int> set_B;
	set_B.insert(2);
	MTBDD leaf_B = make_set_leaf(&set_B);

	std::set<int> set_C;
	set_C.insert(3);
	MTBDD leaf_C = make_set_leaf(&set_C);

	std::set<int> set_D;
	set_D.insert(4);
	MTBDD leaf_D = make_set_leaf(&set_D);

	// Transition t0:
	// Q0---(1, 0) --> A
	// 	 \--(1, 1) --> B
	MTBDD t0_subtree = mtbdd_makenode(b2, leaf_A, leaf_B);
	MTBDD t0_root = mtbdd_makenode(b1, mtbdd_false, t0_subtree);

	// Transition:
	// Q1---(0, *) --> C
	// 	 \--(1, 1) --> D
	MTBDD t1_subtree = mtbdd_makenode(b2, mtbdd_false, leaf_D);
	MTBDD t1_root = mtbdd_makenode(b1, leaf_C, t1_subtree);

	mtbdd_refs_push(t0_root);
	mtbdd_refs_push(t1_root);

	BDDSET vars = mtbdd_set_empty();
	vars = mtbdd_set_add(vars, 2);
	vars = mtbdd_set_add(vars, 1);

	MTBDD tfn_root = mtbdd_apply(t0_root, t1_root, TASK(set_union));
	mtbdd_refs_push(tfn_root);
	// Expected transition fn:
	// {Q0, Q1}----(0, *) --> {C}
	// 		   \---(1, 0) --> {A}
	// 		   \---(1, 1) --> {3, 4}

	{
		// Root low examination
		MTBDD root_low = mtbdd_getlow(tfn_root);
		assert(mtbdd_isleaf(root_low));
		set<int> *root_low_states = (set<int> *) mtbdd_getvalue(root_low);
		cout << "Root (0, *): ";
		print_states_set(root_low_states);
	}

	{
		// Root hight examination
		MTBDD root_high = mtbdd_gethigh(tfn_root);
		assert(!mtbdd_isleaf(root_high));
		{
			MTBDD root_high_high = mtbdd_gethigh(root_high);
			assert(mtbdd_isleaf(root_high_high));
			set<int> *root_low_states = (set<int> *) mtbdd_getvalue(root_high_high);
			cout << "Root (1, 1): ";
			print_states_set(root_low_states); // This should be {3, 4}
			cout << "Refs count: " << sylvan_count_refs() << endl;
		}
		{
			MTBDD root_high_low = mtbdd_getlow(root_high);
			assert(mtbdd_isleaf(root_high_low));
			set<int> *root_low_states = (set<int> *) mtbdd_getvalue(root_high_low);
			cout << "Root (1, 0): ";
			print_states_set(root_low_states);
		}
	}

	mtbdd_refs_pop(3);

	sylvan_gc();
}


VOID_TASK_1(_main, void*, arg)
{
	// args: Table min size, table max size, cache sizes
	sylvan_set_sizes(1LL<<22, 1LL<<26, 1LL<<22, 1LL<<26);
    sylvan_init_package();

	sylvan_init_mtbdd();

	CALL(simple_fn);
	
	sylvan_stats_report(stdout);

	sylvan_quit();
	(void) arg;
}

int main()
{
	int n_workers = 1;
	size_t dequeue_size = 0; 		// Auto select
	size_t program_stack_size = 0;	// Use default
	
	lace_init(n_workers, dequeue_size);

	lace_startup(program_stack_size, TASK(_main), NULL);
}
