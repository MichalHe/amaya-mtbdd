#include <assert.h>
#include <iterator>
#include <sylvan.h>
#include <sylvan_bdd.h>
#include <sylvan_common.h>
#include <sylvan_gmp.h>
#include <sylvan_mt.h>
#include <sylvan_mtbdd.h>
#include <sylvan_int.h>
#include <sylvan_stats.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

using namespace sylvan;
using std::cout;
using std::endl;
using std::set;
;
using std::vector;

static uint32_t mtbdd_leaf_type_set;

void print_states_set(set<int> *states) {
  int cnt = 1;
  cout << "{";
  for (auto state : *states) {
    cout << state;
    if (cnt < states->size()) {
      cout << ", ";
    }
    ++cnt;
  }
  cout << "}" << endl;
}

#ifndef rotl64
static inline uint64_t rotl64(uint64_t x, int8_t r) {
  return ((x << r) | (x >> (64 - r)));
}
#endif

MTBDD
make_set_leaf(std::set<int> *value) {
  MTBDD leaf = mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t)value);
  return leaf;
}

/**
 * Function is called when a new node is being inserted to the hash table
 *
 * @param state_set_ptr 	Pointer to an old valid leaf value, which is being
 * copied to the hash table. This means that its a pointer to a poiter to a
 * vector.(which is our value).
 */
static void mk_set_leaf(uint64_t *state_set_ptr) {
  set<int> *old_states =
      (set<int> *)*state_set_ptr; // After dereferece we get the actual value -
                                  // pointer to a set.
  set<int> *states_copy = new set<int>(*old_states); // Copy construct it

  // Accorting to the GMP leaf implementation, it should suffice to just write
  // the new value to the state_set_ptr.
  *(set<int> **)state_set_ptr = states_copy;
}

static void destroy_set_leaf(uint64_t leaf_value) {
  set<int> *states_set_ptr = (set<int> *)leaf_value;
  print_states_set(states_set_ptr);
  delete states_set_ptr;
}

int set_leaf_equals(uint64_t a_ptr, uint64_t b_ptr) {
  set<int> *states_a = (set<int> *)a_ptr;
  set<int> *states_b = (set<int> *)a_ptr;

  return (*states_a) == (*states_b);
}

static uint64_t set_leaf_hash(const uint64_t contents_ptr,
                              const uint64_t seed) {
  set<int> *states = (set<int> *)contents_ptr;
  const uint64_t prime = 1099511628211;
  uint64_t hash = seed;

  for (auto i : *states) {
    hash = hash ^ i;
    hash = rotl64(hash, 47);
    hash = hash * prime;
  }

  return hash;
}

TASK_DECL_2(MTBDD, set_union, MTBDD *, MTBDD *);
TASK_IMPL_2(MTBDD, set_union, MTBDD *, pa, MTBDD *, pb) {
  MTBDD a = *pa, b = *pb;

  // When one leaf is empty set (false),
  // the algorithm should return all reachable states for the other one
  if (a == mtbdd_false) {
    return b; // If b is also mtbdd_false, nothing wrong is happening
  }
  if (b == mtbdd_false) {
    return a; // Same here
  }

  // If both are leaves, we calculate union
  if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
    std::set<int> set_a = *((std::set<int> *)mtbdd_getvalue(a));
    std::set<int> set_b = *((std::set<int> *)mtbdd_getvalue(b));

    std::set<int> *_union = new std::set<int>();
    std::set_union(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                   std::inserter(*_union, _union->begin()));

    print_states_set(_union);

    MTBDD union_leaf = make_set_leaf(_union);

    return union_leaf;
  }

  // TODO: Perform pointer swap here, so the cache would be utilized
  // (commutative).
  // TODO: Utilize operation cache.
  return mtbdd_invalid;
}

TASK_DECL_2(MTBDD, my_exists_op, MTBDD *, MTBDD *);
TASK_IMPL_2(MTBDD, my_exists_op, MTBDD *, pa, MTBDD *, pb) {
    cout << "Calling my_exists_op" << endl;
    return CALL(set_union, pa, pb);
}

/**
 * The *abstraction* F definition.
 * Task returns a MTBDD. Task accepts two MTBDDs - left and right child (subtree) of a
 * node for variable <v> specified in variable set passed to abstract_apply with this operator.
 * The final int is a number of variables that are missing when reaching a leaf in MTBDD, but there is 
 * still <k> number of variables in the variable set passed to abstract_apply.
 */
TASK_DECL_3(MTBDD, my_abstract_exists_op, MTBDD, MTBDD, int);
TASK_IMPL_3(MTBDD, my_abstract_exists_op, MTBDD, a, MTBDD, b, int, k) {
    if (k == 0) {
        cout << "K is zero, returning my_exists_op" << endl;
        return CALL(set_union, &a, &b);
    } else {
        cout << "K is not zero, returning my_exists_op" << endl;
        return CALL(set_union, &a, &b);
    }
}

TASK_DECL_3(MTBDD, my_abstract_exists, MTBDD, MTBDD, MTBDD);
TASK_IMPL_3(MTBDD, my_abstract_exists, MTBDD, a, MTBDD, b, MTBDD, v) {
    cout << "Inside my_abstract_exists" << endl;
	// Check if there are any more variables in <v>
    if (v == mtbdd_true) {
		// We have exhausted the variables in <v>, no abstraction will be applied down the tree.
		// Just apply the ordinary operator.
        cout << "The variables <v> are exhausted" << endl;
		return mtbdd_apply(a, b, TASK(mtbdd_op_times));
	}

	// Maybe a or b are either terminals (true,false) or leaf nodes. Check
	// this case by trying to apply the operator, which will return mtbdd_invalid
	// if it cannot be applied directly to nodes, and therefor the nodes are not
	// leaves or terminals.
    cout << "Variables <v> are not empty, trying to apply my_exists_op" << endl;
    MTBDD result = CALL(my_exists_op, &a, &b);
    if (result != mtbdd_invalid) {
        cout << "my_exists_op did not return mtbdd_invalid" << endl;
        mtbdd_refs_push(result);
        result = mtbdd_abstract(result, v, TASK(my_abstract_exists_op));
        mtbdd_refs_pop(1);
        return result;
    }

    // At this point we know, that:
	// 	- v is not a constant (contains some variables)  and
    // 	- a or b is not a constant (are subtrees)

    int is_a_leaf = mtbdd_isleaf(a);
    int is_b_leaf = mtbdd_isleaf(b);
    // If the <a> or <b> are not leaves, then they represent a whole MTBDD, but we want only top node.
    // Get the top node.
    mtbddnode_t na = is_a_leaf ? 0 : MTBDD_GETNODE(a);
    mtbddnode_t nb = is_b_leaf ? 0 : MTBDD_GETNODE(b);

    uint32_t variable_a = is_a_leaf ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t variable_b = is_b_leaf ? 0xffffffff : mtbddnode_getvariable(nb);

    // The top nodes of the two subtrees might not be at the same level (might be different vars)
    // Get the one that is higher in the tree / variable order.
    uint32_t var = variable_a < variable_b ? variable_a : variable_b;

    // Get the topmost variable from the <v> variable set.
    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t variable_v = mtbddnode_getvariable(nv);

    if (variable_v < var) {
        cout << "The current variable from <v> is sooner than the one inside tree:" << endl;
        cout << "Variable from <v>:" << variable_v << ", tree var: " << var << endl;
        // We traversed deeper to the tree that the first variable specified in the variable set <v>.
        // What to do about it? We continue recursing even deeper, however after leaving the recursion
        // we must perform the abstraction -- with arguments (result, result) -- both subtrees of the 
        // imaginary missing variable `variable_v` are the same trees = `result`
        
        // Recurse deeper with the next variable from the variable set. 
        result = CALL(my_abstract_exists, a, b, node_gethigh(v, nv));
        mtbdd_refs_push(result);
        result = mtbdd_apply(result, result, TASK(my_exists_op)); // TODO: Change the operand here.
        mtbdd_refs_pop(1);
    } else  {
        MTBDD alow, ahigh, blow, bhigh;
        // TODO: Rewrite this to normal condition.

        // If <a> is a leaf then it does not have successors - both successors must be <a>
        alow  = (!is_a_leaf && variable_a == var) ? node_getlow(a, na)  : a;
        ahigh = (!is_a_leaf && variable_a == var) ? node_gethigh(a, na) : a;

        // If <b> is a leaf then it does not have successors - both successors must be <b>
        blow  = (!is_b_leaf && variable_b == var) ? node_getlow(b, nb)  : b;
        bhigh = (!is_b_leaf && variable_b == var) ? node_gethigh(b, nb) : b;

        // Check if we currently processing one of the desired variables.
        if (variable_v == var) {
            cout << "Variable from <v> and current tree var are on the same level" << endl;
            // Step down the recursion and since we are dealing with the desired variable to be
            // abstracted away we call the operation on it.
            mtbdd_refs_spawn(SPAWN(my_abstract_exists, ahigh, bhigh, node_gethigh(v, nv)));
            MTBDD low = mtbdd_refs_push(CALL(my_abstract_exists, alow, blow, node_gethigh(v, nv)));
            MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(my_abstract_exists)));
            result = CALL(mtbdd_apply, low, high, TASK(mtbdd_op_plus)); // TODO: Change the op here.
            mtbdd_refs_pop(2);
        } else  {
            cout << "Variable from <v> is deeper than the current one:" << endl;
            cout << "Variable from <v>:" << variable_v << ", tree var: " << var << endl;
            // variable_v > var (other possibilities were checked)
            // This means that we are somewhere higher in the tree, and we didn't reach the required variable yet.
            // The action is to go low via both subtrees and hight via both subtrees and just make 
            // new node out of resuts.
            mtbdd_refs_spawn(SPAWN(my_abstract_exists, ahigh, bhigh, v));
            MTBDD low = mtbdd_refs_push(CALL(my_abstract_exists, alow, blow, v));
            MTBDD high = mtbdd_refs_sync(SYNC(my_abstract_exists));
            mtbdd_refs_pop(1);
            result = mtbdd_makenode(var, low, high);
        }
    }

    return result;
}


VOID_TASK_0(simple_fn) {
  // Define custom leaf type
  mtbdd_leaf_type_set = sylvan_mt_create_type();

  sylvan_mt_set_hash(mtbdd_leaf_type_set, set_leaf_hash);
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
    set<int> *root_low_states = (set<int> *)mtbdd_getvalue(root_low);
    print_states_set(root_low_states);
  }

  {
    // Root hight examination
    MTBDD root_high = mtbdd_gethigh(tfn_root);
    assert(!mtbdd_isleaf(root_high));
    {
      MTBDD root_high_high = mtbdd_gethigh(root_high);
      assert(mtbdd_isleaf(root_high_high));
      set<int> *root_low_states = (set<int> *)mtbdd_getvalue(root_high_high);
      print_states_set(root_low_states); // This should be {3, 4}
    }
    {
      MTBDD root_high_low = mtbdd_getlow(root_high);
      assert(mtbdd_isleaf(root_high_low));
      set<int> *root_low_states = (set<int> *)mtbdd_getvalue(root_high_low);
      print_states_set(root_low_states);
    }
  }

  mtbdd_refs_pop(3);

  sylvan_gc();
}

VOID_TASK_1(_main, void *, arg) {
  // args: Table min size, table max size, cache sizes
  sylvan_set_sizes(1LL << 22, 1LL << 26, 1LL << 22, 1LL << 26);
  sylvan_init_package();

  sylvan_init_mtbdd();

  CALL(simple_fn);

  sylvan_stats_report(stdout);

  sylvan_quit();
  (void)arg;
}

int main() {
  int n_workers = 1;
  size_t dequeue_size = 0;       // Auto select
  size_t program_stack_size = 0; // Use default

  lace_init(n_workers, dequeue_size);

  lace_startup(program_stack_size, TASK(_main), NULL);
}
