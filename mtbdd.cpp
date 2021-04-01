#include "amaya_mtbdd.hpp"
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <sylvan_mtbdd_int.h>
#include <unistd.h>

using namespace sylvan;
using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;

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
  //print_states_set(states_set_ptr);
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

static char *set_leaf_to_str(int comp, uint64_t val, char *buf, size_t buflen){
	(void) comp;
	std::stringstream ss;
	auto list_contents = (std::set<int> *) val;
	ss << "{";
	uint32_t cnt = 1;
	for (auto i : *list_contents) {
		ss << i;
		if (cnt < list_contents->size()) {
			ss << ",";
		}
		cnt++;
	}
	ss << "}";
	const std::string str(ss.str());
	
	// Does the resulting string fits into the provided buffer?
	if (str.size() <= buflen) {
		const char *cstr = str.c_str();	
		std::memcpy(buf, cstr, str.size());
		return buf;
	} else {
		char *new_buf = (char *) malloc(sizeof(char) * str.size());
		std::memcpy(buf, new_buf, str.size());
		return new_buf;
	}
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
	// Even if we skipped <k> levels in the mtbdd the k-times applied union()
	//cout << "----- A ------ " << endl;
	//mtbdd_fprintdot(stdout, a);
	//cout << "----- B ------ " << endl;
	//mtbdd_fprintdot(stdout, b);
	MTBDD u = mtbdd_apply(a, b, TASK(set_union));
	//cout << "----- Result: ------ " << endl;
	//mtbdd_fprintdot(stdout, b);
	return u;
}


VOID_TASK_0(simple_fn) {
    // Define custom leaf type
    // sylvan_mt_set_write_binary(gmp_type, gmp_write_binary);
    // sylvan_mt_set_read_binary(gmp_type, gmp_read_binary);

    // We make a simple cube to test our algorithm
    // Transition
    // (1) -1-> (2) --> (3) -1-> {1, 2}

    // REMIND ME: makenode - Low, then high
    std::set<int> leaf_set_a {1, 2};
    std::set<int> leaf_set_b {1, 3};

    MTBDD leaf_a = make_set_leaf(&leaf_set_a); 
    MTBDD leaf_b = make_set_leaf(&leaf_set_b); 

    MTBDD subtree = mtbdd_makenode(3, leaf_b, leaf_a);
	MTBDD root = mtbdd_makenode(1, mtbdd_false, subtree);

    uint32_t arr[] = {1, 2, 3};
    BDDSET v = mtbdd_set_from_array(arr, 3);
   	
	// Project away the first variable.
	MTBDD result = mtbdd_abstract(root, v, TASK(my_abstract_exists_op));
	mtbdd_fprintdot(stdout, result);

    sylvan_gc();
}

VOID_TASK_1(_main, void *, arg) {
  // args: Table min size, table max size, cache sizes
  sylvan_set_sizes(1LL << 22, 1LL << 26, 1LL << 22, 1LL << 26);
  sylvan_init_package();

  sylvan_init_mtbdd();

  //sylvan_quit();
  (void)arg;
}

int main() {
  int n_workers = 1;
  size_t dequeue_size = 0;       // Auto select
  size_t program_stack_size = 0; // Use default

  lace_init(n_workers, dequeue_size);
  lace_startup(program_stack_size, TASK(_main), NULL);
}

void init_machinery() {
  int n_workers = 1;
  size_t dequeue_size = 0;       // Auto select
  size_t program_stack_size = 0; // Use default

  lace_init(n_workers, dequeue_size);
  //lace_startup(program_stack_size, TASK(_main), NULL);

  // THIS SEEMS TO BE THE SECRET
  // When the TASK parameter (the middle one) is set to be NULL, it does not spawn a new thread
  // and instead uses the current thread for all tasks (makes it possible to call from python)
  lace_startup(0, NULL, NULL); 
  sylvan_set_sizes(1LL << 22, 1LL << 26, 1LL << 22, 1LL << 26);
  sylvan_init_package();
  sylvan_init_mtbdd();

  // Initialize my own sylvan type.

    mtbdd_leaf_type_set = sylvan_mt_create_type();
    sylvan_mt_set_hash(mtbdd_leaf_type_set, set_leaf_hash);
    sylvan_mt_set_equals(mtbdd_leaf_type_set, set_leaf_equals);
    sylvan_mt_set_create(mtbdd_leaf_type_set, mk_set_leaf);
    sylvan_mt_set_destroy(mtbdd_leaf_type_set, destroy_set_leaf);
    sylvan_mt_set_to_str(mtbdd_leaf_type_set, set_leaf_to_str);
}

void shutdown_machinery() {
	sylvan_quit();
	lace_exit();
}

MTBDD w_mtbdd_makenode(uint32_t var, MTBDD low, MTBDD high) {
	return mtbdd_makenode(2, mtbdd_true, mtbdd_false);
}

MTBDD w_sylvan_ithvar(uint32_t var) {
	return mtbdd_ithvar(var);
}


MTBDD amaya_mtbdd_build_single_terminal(
		uint8_t *transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t transition_symbols_count,
		uint32_t variable_count,
		uint32_t *destination_set,
		uint32_t destination_set_size)
{
	// CALL(mtbdd_union_cube(mtbdd, variables, cube, terminal)	
	// Cube (transitions symbols) have to be in format -- suppose the transition_symbols are ready
	// x_0 = 0  ===> low transition
	// x_0 = 1  ===> high transition
	// x_0 = 2  ===> any
	// Construct variables set containing all the variables.
	BDDSET variables = mtbdd_set_empty();
	for (uint32_t i=1; i <= variable_count; i++) {
		variables = mtbdd_set_add(variables, i); // Variables are numbered from 1
	}

	// Construct the leaf
	std::set<int> *leaf_state_set = new std::set<int>();	
	for (uint32_t i = 0; i < destination_set_size; i++) {
		leaf_state_set->insert(destination_set[i]);
	}	
	MTBDD leaf = make_set_leaf(leaf_state_set);

	// Construct the initial MTBDD, then add the rest of the symbols
	// signature: mtbdd_cube(MTBDD variables, uint8_t *cube, MTBDD terminal)
	// initial_cube = transition_symbols[0*variable_count + (0..variable_count)] (size: variable_count)
	MTBDD mtbdd = mtbdd_cube(variables, transition_symbols, leaf);
	
	LACE_ME;
	for (uint32_t i = 1; i < transition_symbols_count; i++) {
		// Cube for this iteration:
		//   i*variable_count + (0..variable_count)
		mtbdd = CALL(mtbdd_union_cube, mtbdd, variables, &transition_symbols[i*variable_count], leaf);
	}
	//mtbdd_fprintdot(stdout, mtbdd);
	return mtbdd;
}


int* amaya_mtbdd_get_transition_target(MTBDD mtbdd, uint8_t* cube, uint32_t cube_size, uint32_t* result_size) 
{

	auto search_result = _get_transition_target(mtbdd, 1, cube, cube_size);
	if (search_result == NULL)
	{
		*result_size = 0;
		return NULL;
	} else 
	{
		const uint32_t rs = search_result->size();
		*result_size = rs;
		int* result_arr = (int*) malloc(sizeof(uint32_t) * rs);
		uint32_t i = 0;
		for (auto v : *search_result) 
		{
			result_arr[i++] = v;
		}
		return result_arr;
	}
}

void amaya_do_free(void *ptr)
{
	free(ptr);
}

std::set<int>* _get_transition_target(
		MTBDD root, 
		uint32_t current_variable,
		uint8_t* variable_assigments, 
		uint32_t var_count)
{
	int is_leaf = mtbdd_isleaf(root);
	if (var_count == 0) {
		if (is_leaf) {
			if (root == mtbdd_false) return NULL;
			auto leaf_value = (std::set<int>*) mtbdd_getvalue(root);	
			std::set<int>* result = new std::set<int>(*leaf_value); // Return a copy (might modify elsewhere)
			return result;
		} else {
			// We have exhausted all variables, yet we did not hit a leaf. 
			// That means the cube did not contain all the variables.
			return NULL;  // Not enough information to decide.
		}
	}
	
	if (is_leaf) {
		if (root == mtbdd_false) return NULL;
		else {
			// Found a solution
			auto leaf_value = (std::set<int>*) mtbdd_getvalue(root);	
			std::set<int>* result = new std::set<int>(*leaf_value); // Return a copy (might modify elsewhere)
			return result;
		}
	}
	
	// If we got here then root is not a leaf and we still have still some variables
	uint32_t root_var = mtbdd_getvar(root);

	if (root_var == current_variable) {
		if (*variable_assigments == 0) 
		{
			// Go low
			return _get_transition_target(
					mtbdd_getlow(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);
		} 
		else if (*variable_assigments == 1) 
		{
			// Go high
			return _get_transition_target(
					mtbdd_gethigh(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);
		}
		else 
		{
			// That means that the current assigment is don't care (2)
			// First go low
			std::set<int>* low_result = _get_transition_target(
					mtbdd_getlow(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);

			// Then go high
			std::set<int>* high_result = _get_transition_target(
					mtbdd_gethigh(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);

			// Now decide what to do with the results
			if (low_result == NULL) return high_result; // Might return NULL
			if (high_result == NULL) return low_result; // The low result is not NULL

			low_result->insert(high_result->begin(), high_result->end());
			delete high_result; // Only low result will be propagated upwards
			return low_result;
		}
	} else 
	{
		// The current variable in the tree skipped a variable.
		// That means the current variable assignment did not matter.
		return _get_transition_target(
				root, 
				current_variable + 1,   		// Move to the next variable
				variable_assigments + 1,		// --^
				var_count - 1);
	}
}


void amaya_print_dot(MTBDD m, int32_t fd) {
	FILE* file = fdopen(fd, "w");
	cout << file << endl;
	mtbdd_fprintdot(file, m);
	fflush(file);
	// fclose(file);  The file will be closed in python
}

MTBDD amaya_project_variables_away(MTBDD m, uint32_t *variables, uint32_t var_count) {
	// Construct the variables set.
	BDDSET var_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < var_count; i++) {
		var_set = mtbdd_set_add(var_set, variables[i]);
	}
	
	// Do the projection itself.
	LACE_ME;
	MTBDD result = mtbdd_abstract(m, var_set, TASK(my_abstract_exists_op));
	
	//mtbdd_fprintdot(stdout, result);
	return result;
}

MTBDD amaya_unite_mtbdds(MTBDD m1, MTBDD m2) {
	LACE_ME;
	MTBDD u = mtbdd_apply(m1, m2, TASK(set_union));
	//mtbdd_fprintdot(stdout, u);

	return u;
}
