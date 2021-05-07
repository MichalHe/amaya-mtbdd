#include "amaya_mtbdd.hpp"
#include <algorithm>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <sylvan_mtbdd_int.h>
#include <unistd.h>
#include <utility>

using namespace sylvan;
using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::map;
using std::pair;

static map<pair<int, int>, int>* intersection_state = NULL;
static bool DEBUG_ON = false;

Transition_Destination_Set::Transition_Destination_Set() 
{
	this->destination_set = NULL;
	this->automaton_id = 0;
}

Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other) 
{
    this->destination_set = new std::set<int>(*other.destination_set); // Do copy
    this->automaton_id = other.automaton_id;
}

Transition_Destination_Set::Transition_Destination_Set(uint32_t automaton_id, std::set<int>* destination_set) 
{
	this->automaton_id = automaton_id;
    this->destination_set = destination_set;
}

Transition_Destination_Set::~Transition_Destination_Set() 
{
	if (this->destination_set != NULL) {
		delete this->destination_set;
	}
}

void Transition_Destination_Set::print_dest_states() 
{
    int cnt = 1;
    cout << "{";
    for (auto state : *this->destination_set) {
      cout << state;
      if (cnt < this->destination_set->size()) {
        cout << ", ";
      }
      ++cnt;
    }
    cout << "}" << endl;
}


static uint32_t mtbdd_leaf_type_set;


#ifndef rotl64
static inline uint64_t rotl64(uint64_t x, int8_t r) 
{
  return ((x << r) | (x >> (64 - r)));
}
#endif

MTBDD
make_set_leaf(Transition_Destination_Set* value) {
    MTBDD leaf = mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) value);
    return leaf;
}

/**
 * Function is called when a new node is being inserted to the hash table
 *
 * @param state_set_ptr 	Pointer to an old valid leaf value, which is being
 * copied to the hash table. This means that its a pointer to a poiter to a
 * vector.
 */
static void mk_set_leaf(uint64_t *target_destination_set_ptr) 
{
    Transition_Destination_Set* original_tds = (Transition_Destination_Set *) *target_destination_set_ptr; 
    Transition_Destination_Set* new_tds = new Transition_Destination_Set(*original_tds); // Copy construct

    // Accorting to the GMP leaf implementation, it should suffice to just write
    // the new value to the state_set_ptr.
    *(Transition_Destination_Set **)target_destination_set_ptr = new_tds;
}

static void destroy_set_leaf(uint64_t leaf_value) 
{
    Transition_Destination_Set *tds = (Transition_Destination_Set *)leaf_value;
    delete tds;
}

int set_leaf_equals(uint64_t a_ptr, uint64_t b_ptr) 
{
    Transition_Destination_Set* a_tds = (Transition_Destination_Set*) a_ptr;
    Transition_Destination_Set* b_tds = (Transition_Destination_Set*) b_ptr;

    if (a_tds->automaton_id == b_tds->automaton_id) {
        if ((*a_tds->destination_set) == (*b_tds->destination_set)) {
            return true;
        };
    }
    return false;
}

static uint64_t set_leaf_hash(const uint64_t contents_ptr,
                              const uint64_t seed) 
{
    Transition_Destination_Set* tds = (Transition_Destination_Set*) contents_ptr;
    const uint64_t prime = 1099511628211;
    uint64_t hash = seed;

    for (auto i : *tds->destination_set) {
        hash = hash ^ i;
        hash = rotl64(hash, 47);
        hash = hash * prime;
    }

	// Hash also the automaton_id, TODO: Remove this?
	//hash = hash ^ tds->automaton_id;
	//hash = rotl64(hash, 31);
	//hash = hash * prime;

    return hash;
}

static char *set_leaf_to_str(int comp, uint64_t leaf_val, char *buf, size_t buflen){
	(void) comp;
	std::stringstream ss;
	auto tds = (Transition_Destination_Set*) leaf_val;
    ss << "#" <<tds->automaton_id << " ";
	ss << "{";
	uint32_t cnt = 1;
	for (auto i : *tds->destination_set) {
		ss << i;
		if (cnt < tds->destination_set->size()) {
			ss << ",";
		}
		cnt++;
	}
	ss << "}";

	const std::string str(ss.str());
	
	// Does the resulting string fits into the provided buffer?
	const size_t required_buf_size = str.size() + 1; // With nullbyte
	if (required_buf_size <= buflen) {
		const char *cstr = str.c_str();	
		std::memcpy(buf, cstr, sizeof(char) * required_buf_size);
		return buf;
	} else {
		char *new_buf = (char *) malloc(sizeof(char) * required_buf_size);
		std::memcpy(new_buf, str.c_str(), sizeof(char) * required_buf_size);
		return new_buf;
	}
}


TASK_DECL_3(MTBDD, set_union, MTBDD*, MTBDD*, uint64_t);
TASK_IMPL_3(MTBDD, set_union, MTBDD*, pa, MTBDD*, pb, uint64_t, param) 
{
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
    auto& tds_a = *((Transition_Destination_Set*) mtbdd_getvalue(a));
    auto& tds_b = *((Transition_Destination_Set*) mtbdd_getvalue(b));

    std::set<int> *union_set = new std::set<int>();
    std::set_union(
            tds_a.destination_set->begin(), tds_a.destination_set->end(),
            tds_b.destination_set->begin(), tds_b.destination_set->end(),
            std::inserter(*union_set, union_set->begin()));

    Transition_Destination_Set* union_tds = new Transition_Destination_Set((uint32_t) param, union_set);

    MTBDD union_leaf = make_set_leaf(union_tds);  // Wrap the TDS with a MTBDD leaf.
    return union_leaf;
  }

  // TODO: Perform pointer swap here, so the cache would be utilized
  // (commutative).
  // TODO: Utilize operation cache.
  return mtbdd_invalid;
}

TASK_DECL_3(MTBDD, set_intersection_op, MTBDD *, MTBDD *, uint64_t);
TASK_IMPL_3(MTBDD, set_intersection_op, MTBDD *, pa, MTBDD *, pb, uint64_t, param) 
{
    MTBDD a = *pa, b = *pb;
	// Intersection with an empty set (mtbdd_false) is empty set (mtbdd_false)
    if (a == mtbdd_false || b == mtbdd_false) {
        return mtbdd_false;
    }

    // If both are leaves calculate intersection
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
		auto intersect_info = (intersect_info_t*) param;

        auto& tds_a = *((Transition_Destination_Set*) mtbdd_getvalue(a));
        auto& tds_b = *((Transition_Destination_Set*) mtbdd_getvalue(b));

        std::set<int> *left_states  = tds_a.destination_set;  
        std::set<int> *right_states = tds_b.destination_set;  

        // Calculate cross product
        pair<int, int> metastate; 
        int state;

        std::set<int>* intersection_leaf_states = new std::set<int>();

        for (auto left_state : *left_states) {
            for (auto right_state : *right_states) {
                metastate = std::make_pair(left_state, right_state);
                auto pos = intersection_state->find(metastate);
                bool contains_metastate = pos != intersection_state->end();
                if (contains_metastate) {
                    state = pos->second;
                } else {
                    // We have discovered a new state.
                    // Update the global intersection state, so in the future every such
                    // state will get the same integer
                    state = intersection_state->size();
                    intersection_state->insert(std::make_pair(metastate, state));

                    // Update the discoveries local for this intersection, so that the Python
                    // side knows whats up.
                    intersect_info->discoveries->push_back(left_state);
                    intersect_info->discoveries->push_back(right_state);
                    intersect_info->discoveries->push_back(state);
                }

                intersection_leaf_states->insert(state); 
            }
        }

        if (intersection_leaf_states->empty()) {
            delete intersection_leaf_states;
            return mtbdd_false;
        }
         
        auto intersection_tds = new Transition_Destination_Set(intersect_info->automaton_id, intersection_leaf_states);
        MTBDD intersection_leaf = make_set_leaf(intersection_tds);
        return intersection_leaf;
    }

    // TODO: Perform pointer swap here, so the cache would be utilized
    // (commutative).
    // TODO: Utilize operation cache.
    return mtbdd_invalid;
}

/**
 * The *abstraction* F definition.
 * Task returns a MTBDD. Task accepts two MTBDDs - left and right child (subtree) of a
 * node for variable <v> specified in variable set passed to abstract_apply with this operator.
 * The final int is a number of variables that are missing when reaching a leaf in MTBDD, but there is 
 * still <k> number of variables in the variable set passed to abstract_apply.
 */
TASK_DECL_3(MTBDD, my_abstract_exists_op, MTBDD, MTBDD, int);
TASK_IMPL_3(MTBDD, my_abstract_exists_op, MTBDD, a, MTBDD, b, int, k) 
{
	// Even if we skipped <k> levels in the mtbdd the k-times applied union()
	//cout << "----- A ------ " << endl;
	//mtbdd_fprintdot(stdout, a);
	//cout << "----- B ------ " << endl;
	//mtbdd_fprintdot(stdout, b);
	MTBDD u = mtbdd_applyp(a, b, (uint64_t) -1, TASK(set_union), AMAYA_EXISTS_OPERATION_ID);
	//cout << "----- Result: ------ " << endl;
	//mtbdd_fprintdot(stdout, b);
	return u;
}

void init_machinery() 
{
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

void shutdown_machinery() 
{
	sylvan_quit();
	lace_exit();
}

MTBDD w_mtbdd_makenode(uint32_t var, MTBDD low, MTBDD high) 
{
	return mtbdd_makenode(2, mtbdd_true, mtbdd_false);
}

MTBDD w_sylvan_ithvar(uint32_t var) 
{
	return mtbdd_ithvar(var);
}


MTBDD amaya_mtbdd_build_single_terminal(
		uint32_t  automaton_id,
		uint8_t*  transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t  transition_symbols_count,
		uint32_t  variable_count,
		int* 	  destination_set,
		uint32_t  destination_set_size)
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
    
	// Construct the destination set
	std::set<int> *leaf_state_set = new std::set<int>();	
	for (uint32_t i = 0; i < destination_set_size; i++) {
		leaf_state_set->insert(destination_set[i]);
	}
    
    Transition_Destination_Set* tds = new Transition_Destination_Set(automaton_id, leaf_state_set);
	MTBDD leaf = make_set_leaf(tds);

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

	return mtbdd;
}


int* amaya_mtbdd_get_transition_target(
        MTBDD mtbdd, 
        uint8_t* cube, 
        uint32_t cube_size, 
        uint32_t* result_size) 
{

	auto search_result_tds = _get_transition_target(mtbdd, 1, cube, cube_size);
	if (search_result_tds == NULL)
	{
		*result_size = 0;
		return NULL;
	} else 
	{
		const uint32_t rs = search_result_tds->destination_set->size();
		*result_size = rs;
		int* result_arr = (int*) malloc(sizeof(uint32_t) * rs);
		uint32_t i = 0;
		for (auto v : *search_result_tds->destination_set) 
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

Transition_Destination_Set* _get_transition_target(
		MTBDD root, 
		uint32_t current_variable,
		uint8_t* variable_assigments, 
		uint32_t var_count)
{
	int is_leaf = mtbdd_isleaf(root);
	if (var_count == 0) {
		if (is_leaf) {
			if (root == mtbdd_false) return NULL;
			auto tds = (Transition_Destination_Set*) mtbdd_getvalue(root);	
			auto tds_copy = new Transition_Destination_Set(*tds);
			return tds_copy;
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
			auto tds = (Transition_Destination_Set*) mtbdd_getvalue(root);	
			auto tds_copy = new Transition_Destination_Set(*tds); // Return a copy (might modify elsewhere)
			return tds_copy;
		}
	}
	
	// If we got here then root is not a leaf and we still have some variables
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
            // Then we need to explore both branches and return the results
			auto low_tds = _get_transition_target(
					mtbdd_getlow(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);

			auto high_tds = _get_transition_target(
					mtbdd_gethigh(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);

			// Now decide what to do with the results
			if (low_tds == NULL) return high_tds; // Might return NULL
			if (high_tds == NULL) return low_tds; // The low result is not NULL

            // Perform state merging.
			low_tds->destination_set->insert(
                    high_tds->destination_set->begin(),
                    high_tds->destination_set->end());
			delete high_tds; // Only low result will be propagated upwards
			return low_tds;
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

void amaya_mtbdd_rename_states(
		MTBDD* mtbdd_roots, 
		uint32_t root_count,
		int* names, // [(old, new), (old, new), (old, new)] 
		uint32_t name_count)
{
	int old_name, new_name;

	 // No MTBDD interference can happen, since the leaves contain automaton_id, and therefore
	 // each MTBDD has unique leaves.
	std::set<MTBDD> leaves {};
	for (uint32_t i = 0; i < root_count; i++) {
		// This only inserts them into leaves, does not clear the container
		collect_mtbdd_leaves(mtbdd_roots[i], leaves); 
	}

	if (leaves.empty()) {
		return;
	}

	for (auto leaf : leaves) {
		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
	}

    for (uint32_t i = 0; i < name_count; i++) 
    {
        old_name = names[2*i];
        new_name = names[2*i + 1];
    }

	
	for (auto leaf : leaves) {

		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);

		auto* new_leaf_contents = new std::set<int>();
		for (uint32_t i = 0; i < name_count; i++) 
		{
			old_name = names[2*i];
			new_name = names[2*i + 1];
			bool is_in_leaf = (tds->destination_set->find(old_name) != tds->destination_set->end());

			if (is_in_leaf) {
				new_leaf_contents->insert(new_name);
			}
		}
        
		delete tds->destination_set;
		tds->destination_set = new_leaf_contents;
	}
}

void collect_mtbdd_leaves(MTBDD root, std::set<MTBDD>& dest)
{
	LACE_ME;
	MTBDD support = mtbdd_support(root);
	uint32_t support_size = mtbdd_set_count(support);
		
	// Stores the tree path to the leaf
	uint8_t* arr = (uint8_t*) malloc(sizeof(uint8_t) * support_size);

	MTBDD leaf = mtbdd_enum_first(root, support, arr, NULL);
	
 	while (leaf != mtbdd_false)
	{
		dest.insert(leaf);
 	    leaf = mtbdd_enum_next(root, support, arr, NULL);
 	}
}

int* amaya_mtbdd_get_leaves(
		MTBDD root, 
		uint32_t** leaf_sizes,	// OUT, Array containing the sizes of leaves inside dest
		uint32_t*  leaf_cnt,	// OUT, Number of leaves in the tree
        void***    leaf_ptrs)   // Pointer to an array of pointers to MTBDD leaves
{
	std::set<MTBDD> leaves {};
	collect_mtbdd_leaves(root, leaves);	

	// This probably is not the most efficient way how to do it.	
	// First compute the destination size (for malloc)
	uint32_t size_cnt = 0; // Counter for the total number of states.
	for (MTBDD leaf : leaves) 
	{
		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		size_cnt += tds->destination_set->size();
	}

	// Do the allocations
	uint32_t* _leaf_sizes = (uint32_t*) malloc(sizeof(uint32_t) * leaves.size());  // One slot for each leaf 
	int* _leaf_states = (int *) malloc(sizeof(int) * size_cnt);
	
	// Populate
	uint32_t state_i = 0;
	uint32_t size_i = 0;
	for (MTBDD leaf : leaves)
	{
		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (int state : *tds->destination_set) 
		{
			_leaf_states[state_i] = state;
			state_i++;
		}
		_leaf_sizes[size_i] = tds->destination_set->size();
		size_i++;
	}
    if (leaf_ptrs != NULL) {
        void** leaf_ptr_arr = (void**) malloc(sizeof(void*) * leaves.size());    
        
        uint32_t i = 0;
        for (auto leaf : leaves) {
            Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
            leaf_ptr_arr[i] = (void*) tds;
            i++;
        }
        *leaf_ptrs = leaf_ptr_arr;
    }
	
	*leaf_sizes = _leaf_sizes;
	*leaf_cnt = (uint32_t) leaves.size();
	return _leaf_states;
}

void amaya_replace_leaf_contents_with(void *leaf_tds, int* new_contents, uint32_t contents_size)
{
    auto tds = (Transition_Destination_Set*) leaf_tds;
    tds->destination_set->clear();
    for (uint32_t i = 0; i < contents_size; i++) {
        tds->destination_set->insert(new_contents[i]);
    }
}


void amaya_print_dot(MTBDD m, int32_t fd) 
{
	FILE* file = fdopen(fd, "w");
	mtbdd_fprintdot(file, m);
	fflush(file);
	// fclose(file);  The file will be closed in python
}

MTBDD amaya_project_variables_away(MTBDD m, uint32_t *variables, uint32_t var_count) 
{
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

int* amaya_mtbdd_get_state_post(MTBDD m, uint32_t *post_size)
{
	std::set<MTBDD> leaves;
	collect_mtbdd_leaves(m, leaves);
	
	std::set<int> state_post_set;

	for (auto leaf: leaves)
	{
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		state_post_set.insert(
                tds->destination_set->begin(),
                tds->destination_set->end());
	}
	
	int* result = (int*) malloc(sizeof(int) * state_post_set.size());
	uint32_t i = 0;
	for (int state: state_post_set)
	{
		result[i++] = state;
	}

	*post_size = state_post_set.size();
	return result;
}


TASK_DECL_3(MTBDD, pad_closure_op, MTBDD *, MTBDD *, uint64_t);
TASK_IMPL_3(MTBDD, pad_closure_op, MTBDD *, p_left, MTBDD *, p_right, uint64_t, op_param) {
	MTBDD left = *p_left, right = *p_right;
    if (left == mtbdd_false) {
        return right; // Padding closure has no effect, when one mtbdd path leads to nothing
    }
    if (right == mtbdd_false) {
        return left; // Same here
    }
    // If both are non-false leaves, we can try doing the actual padding closure
    if (mtbdd_isleaf(left) && mtbdd_isleaf(right)) {
        auto left_tds  = (Transition_Destination_Set*) mtbdd_getvalue(left);
        auto right_tds = (Transition_Destination_Set*) mtbdd_getvalue(right);

	    pad_closure_info_t* pci = (pad_closure_info_t*) op_param;

		// Check whether the transition destination state even leads the the right state (might not)
		if (left_tds->destination_set->find(pci->right_state) == left_tds->destination_set->end()) {
			// Does not lead to the right state, therefore cannot propagate padding.
			return left; 
		}

	    bool is_final;
	    for (int rs: *right_tds->destination_set)
	    {
			is_final = false;

			// Is the current right state final?
			for (uint32_t i = 0; i < pci->final_states_cnt; i++) {
				if (rs == pci->final_states[i]) {
					is_final = true;
					break; // We have the information we required, we don't need to iterate further.
				}
			}
			
			if (is_final) {
				const bool is_missing_from_left = left_tds->destination_set->find(rs) == left_tds->destination_set->end();
				if (is_missing_from_left) {
					// We've located a state that is final, and is not in left, pad closure adds it to the left states.
					left_tds->destination_set->insert(rs); // The actual pad closure
					pci->had_effect = true; // To utilize the python cache
				}
			}
	    }

		// The pad closure is in place modification 
		// (no new MTBDD is created, only leaf states are modified)
        return left;
    }

    return mtbdd_invalid;
}


bool amaya_mtbdd_do_pad_closure(int left_state, MTBDD left, int right_state, MTBDD right, int* final_states, uint32_t final_states_cnt)
{
	// We need to invalidate MTBDD cache here - for the same pair of mtbdds can the 
	// padding closure be called multiple times as the final states bubble up the state-pre chain
	// TODO: This this is kindof inefficient - maybe try populating the cache with purposely
	// different value before the pad closure, so that the `applyp` will not recognize it???
	cache_clear();
	pad_closure_info_t pci = {0};
	pci.final_states = final_states;
	pci.final_states_cnt = final_states_cnt;
	pci.had_effect = false;

	pci.right_state = right_state;
	pci.left_state = left_state;
	
	LACE_ME;
	mtbdd_applyp(
			left, 
			right, 
			(uint64_t) &pci, 
			TASK(pad_closure_op), 
			AMAYA_PAD_CLOSURE_OPERATION_ID);

	return pci.had_effect;
}

uint8_t* amaya_mtbdd_get_transitions(
		MTBDD root, 	// The mtbdd
		uint32_t* vars,	
		uint32_t var_count,
		uint32_t* symbols_cnt,
		int** dest_states,
		uint32_t** dest_states_sizes)
{
	LACE_ME;
	
	// Construct MTBDD set from given vars
	MTBDD variable_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < var_count; i++) {
		variable_set= mtbdd_set_add(variable_set, vars[i]);
	}

	// Stores the tree path to the leaf
	uint8_t* arr = (uint8_t*) malloc(sizeof(uint8_t) * var_count);

	MTBDD leaf = mtbdd_enum_first(root, variable_set, arr, NULL);
	
	std::vector<int> dest_states_vec;
	std::vector<uint8_t> symbols;
	std::vector<uint32_t> state_sizes;
 	while (leaf != mtbdd_false)
	{
		// Add the destination states to the oveall vector
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (int state: *tds->destination_set) {
			dest_states_vec.push_back(state);
		}

		for (uint32_t i = 0; i < var_count; i++) {
			symbols.push_back(arr[i]);
		}

		state_sizes.push_back(tds->destination_set->size());

 	    leaf = mtbdd_enum_next(root, variable_set, arr, NULL);
 	}

	// Create output arrays
	uint8_t* symbols_out = (uint8_t*) malloc(sizeof(uint8_t) * symbols.size());
	for (uint32_t i = 0; i < symbols.size(); i++) symbols_out[i] = symbols.at(i);

	uint32_t* dest_states_sizes_out = (uint32_t*) malloc(sizeof(uint32_t) * state_sizes.size());
	for (uint32_t i = 0; i < state_sizes.size(); i++) dest_states_sizes_out[i] = state_sizes.at(i);

	int* dest_states_out = (int*) malloc(sizeof(int) * dest_states_vec.size());
	for (uint32_t i = 0; i < dest_states_vec.size(); i++) dest_states_out [i] = dest_states_vec.at(i);

	// Do the return arrays assignment
	*dest_states = dest_states_out;
	*dest_states_sizes = dest_states_sizes_out;

	*symbols_cnt = state_sizes.size(); // For every destination size, there exists 1 symbol leading to it.
	return symbols_out;
}

void amaya_mtbdd_change_automaton_id_for_leaves(
        MTBDD* roots,
        uint32_t root_cnt,
        uint32_t new_id)
{
	std::set<MTBDD> leaves {};
    for (uint32_t i = 0; i < root_cnt; i++) {
        collect_mtbdd_leaves(roots[i], leaves);
    }

    for (auto leaf : leaves) {
        auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf); 
        tds->automaton_id = new_id;
    }
}

MTBDD amaya_mtbdd_intersection(
        MTBDD a, 
        MTBDD b,
        uint32_t result_automaton_id, 
        int** discovered_states,         // OUT
        int*  discovered_states_cnt)     // OUT
{
    auto intersect_info = (intersect_info_t*) malloc(sizeof(intersect_info_t)); // TODO: replace malloc with new
    intersect_info->automaton_id = result_automaton_id;
    intersect_info->discoveries = new std::vector<int>();

	LACE_ME;
	MTBDD result = mtbdd_applyp(a, b, (uint64_t) intersect_info, TASK(set_intersection_op), AMAYA_INTERSECTION_OP_ID);
    
    
    if (!intersect_info->discoveries->empty()) {
        auto discovered_states_arr = (int*) malloc(sizeof(int) * intersect_info->discoveries->size());
        for (uint32_t i = 0; i < intersect_info->discoveries->size(); i++) {
            discovered_states_arr[i] = intersect_info->discoveries->at(i);
        }

        // Send the new discoveries to the python side.
        *discovered_states = discovered_states_arr;
        *discovered_states_cnt = intersect_info->discoveries->size() / 3;
    } else {
        // Nothing new was discovered, do not malloc
        *discovered_states = NULL;
        *discovered_states_cnt = 0;
    }

    delete intersect_info->discoveries;
    free(intersect_info);

    return result;
}

MTBDD amaya_unite_mtbdds(MTBDD m1, MTBDD m2, uint32_t automaton_id) {
	LACE_ME;
	MTBDD u = mtbdd_applyp(m1, m2, (uint64_t) automaton_id, TASK(set_union), AMAYA_UNION_OP_ID);
	return u;
}

void amaya_begin_intersection() 
{
    intersection_state = new map<pair<int, int>, int>();
}

void amaya_update_intersection_state(int* metastates, int* renamed_metastates, uint32_t cnt)
{
    for (uint32_t i = 0; i < cnt; i++) {
        const auto metastate = std::make_pair(metastates[2*i], metastates[2*i + 1]);
        intersection_state->insert(std::make_pair(metastate, renamed_metastates[i])); 
    }
}

void amaya_end_intersection() 
{
    delete intersection_state;
    intersection_state = NULL;
}

void amaya_set_debugging(bool debug) {
	DEBUG_ON = debug;
}
