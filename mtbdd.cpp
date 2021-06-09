#include "amaya_mtbdd.hpp"
#include <algorithm>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <sylvan_mtbdd_int.h>
#include <unistd.h>
#include <assert.h>
#include <unordered_set>
#include <lace.h>
#include <utility>

using namespace sylvan;
using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::stringstream;
using std::map;
using std::pair;
using std::unordered_set;


static Intersection_State* intersection_state = NULL;
static bool DEBUG_ON = false;
static Amaya_Stats STATS = {0};
static uint32_t mtbdd_leaf_type_set;

Transition_Destination_Set::Transition_Destination_Set() 
{
	this->destination_set = NULL;
	this->automaton_id = 0;
}

Transition_Destination_Set::Transition_Destination_Set(const Transition_Destination_Set &other) 
{
    this->destination_set = new std::set<State>(*other.destination_set); // Do copy
    this->automaton_id = other.automaton_id;
}

Transition_Destination_Set::Transition_Destination_Set(uint32_t automaton_id, std::set<State>* destination_set) 
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


MTBDD
make_set_leaf(Transition_Destination_Set* value) {
    MTBDD leaf = mtbdd_makeleaf(mtbdd_leaf_type_set, (uint64_t) value);
    return leaf;
}

/**
 * Function is called when a new node is being inserted to the hash table.
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
	
	// Hash the automaton id too
	hash = hash ^ tds->automaton_id;
	hash = rotl64(hash, 31);
	hash = hash * prime;

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
  if (a == mtbdd_false && b == mtbdd_false) return mtbdd_false;

  // When one leaf is empty set (false),
  // the algorithm should return all reachable states for the other one
  if (a == mtbdd_false) {
    return b; 
  }

  if (b == mtbdd_false) {
	return a; 
  }
  // If both are leaves, we calculate union
  if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
    auto& tds_a = *((Transition_Destination_Set*) mtbdd_getvalue(a));
    auto& tds_b = *((Transition_Destination_Set*) mtbdd_getvalue(b));
		
	// If the passed parameter is missing (-1) that means that we should derive it 
	// from the leaves (TDS). This can happen when padding closure is performed, as
	// the leaf automaton id is not modified.
	if (param == -1) {
		assert(tds_a.automaton_id == tds_b.automaton_id);
		param = (uint32_t) tds_a.automaton_id;
	}

    std::set<State> *union_set = new std::set<State>();
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
		auto intersect_info = (Intersection_Op_Info*) param;

        auto& tds_a = *((Transition_Destination_Set*) mtbdd_getvalue(a));
        auto& tds_b = *((Transition_Destination_Set*) mtbdd_getvalue(b));

        auto left_states  = tds_a.destination_set;  
        auto right_states = tds_b.destination_set;  

		if (left_states->empty() || right_states->empty()) {
			return mtbdd_false;
		}

        // Calculate cross product
        pair<State, State> metastate; 
        State state;

        auto intersection_leaf_states = new std::set<State>();
		
		auto already_discovered_intersection_states = intersection_state->intersection_state_pairs_numbers;

        for (auto left_state : *left_states) {
            for (auto right_state : *right_states) {
                metastate = std::make_pair(left_state, right_state);
                auto pos = already_discovered_intersection_states->find(metastate);
                bool contains_metastate = (pos != already_discovered_intersection_states->end());
                if (contains_metastate) {
                    state = pos->second;
                } else {
                    // We have discovered a new state.
					// Check (if early pruning is on) whether the state should be pruned.  
					if (intersection_state->should_do_early_prunining) {
						const bool is_left_in_pruned = (intersection_state->prune_final_states->find(left_state) != intersection_state->prune_final_states->end());
						const bool is_right_in_pruned = (intersection_state->prune_final_states->find(right_state) != intersection_state->prune_final_states->end());
						
						// Pruning is performed only when exactly one of the states is in the pruned set
						if (is_left_in_pruned != is_right_in_pruned) {
							// The state should be pruned
							STATS_INC(STATS.intersection_states_pruned);
							continue;
						}
					}

                    // Update the global intersection state, so in the future every such
                    // state will get the same integer
                    state = already_discovered_intersection_states->size();
                    already_discovered_intersection_states->insert(std::make_pair(metastate, state));

                    // Update the discoveries local for this intersection, so that the Python
                    // side knows whats up.
                    intersect_info->discoveries->push_back(left_state);
                    intersect_info->discoveries->push_back(right_state);
                    intersect_info->discoveries->push_back(state);
                }

                intersection_leaf_states->insert(state); 
            }
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

TASK_DECL_2(MTBDD, remove_states_op, MTBDD, uint64_t);
TASK_IMPL_2(MTBDD, remove_states_op, MTBDD, dd, uint64_t, param) {
	if (dd == mtbdd_false) return mtbdd_false;

	if (mtbdd_isleaf(dd)) {
		// Param is ignored, instead use the global variable REMOVE_STATES_OP_PARAM
		(void) param;

		auto states_to_remove = (set<State>*) (REMOVE_STATES_OP_PARAM);
		auto tds = (Transition_Destination_Set *) mtbdd_getvalue(dd);

		auto new_tds = new Transition_Destination_Set(*tds); // Make leaf value copy.
		
		for (auto state : *tds->destination_set) {
			bool should_be_removed = states_to_remove->find(state) != states_to_remove->end();	
			if (should_be_removed) new_tds->destination_set->erase(state);
		}

		if (new_tds->destination_set->empty()) {
			delete new_tds;
			return mtbdd_false;
		}
		
		return make_set_leaf(new_tds);
	}

	return mtbdd_invalid;
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
    sylvan_set_sizes(1LL << 23, 1LL << 27, 1LL << 23, 1LL << 27);
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

MTBDD amaya_mtbdd_build_single_terminal(
		uint32_t  automaton_id,
		uint8_t*  transition_symbols,  // 2D array of size (variable_count) * transition_symbols_count
		uint32_t  transition_symbols_count,
		uint32_t  variable_count,
		State* 	  destination_set,
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
	auto leaf_state_set = new std::set<State>();
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

TASK_DECL_2(MTBDD, complete_mtbdd_with_trapstate_op, MTBDD, uint64_t);
TASK_IMPL_2(MTBDD, complete_mtbdd_with_trapstate_op, MTBDD, dd, uint64_t, param)
{
	// param is used just to avoid sylvan miss-cachces
	if (dd == mtbdd_false) {
		(void) param;
		auto op_info = (Complete_With_Trapstate_Op_Info*) ADD_TRAPSTATE_OP_PARAM;

		Transition_Destination_Set* tds = new Transition_Destination_Set();
		tds->automaton_id = op_info->automaton_id;
		tds->destination_set = new set<State>();
		tds->destination_set->insert(op_info->trapstate);

		op_info->had_effect = true;
		return make_set_leaf(tds);
	} else if (mtbdd_isleaf(dd)) {
		// @TODO: Check that the complete with trapstate does have the same automaton ID as the node 
		// being returned.
		return dd;
	}
	return mtbdd_invalid;
}


MTBDD amaya_complete_mtbdd_with_trapstate(
		MTBDD 	 dd, 
		uint32_t automaton_id, 
		State 	 trapstate,
		bool* 	 had_effect) 
{

	Complete_With_Trapstate_Op_Info op_info = {};
	op_info.trapstate = trapstate;
	op_info.automaton_id = automaton_id;
	op_info.had_effect = false;
	
	ADD_TRAPSTATE_OP_PARAM = &op_info;
	LACE_ME;
	MTBDD result = mtbdd_uapply(dd, TASK(complete_mtbdd_with_trapstate_op), ADD_TRAPSTATE_OP_COUNTER);
	ADD_TRAPSTATE_OP_COUNTER++;

	if (DEBUG_ON) printf("Had the complete with trapstate effect? %s\n", (op_info.had_effect ? "Yes": "No"));
	*had_effect = op_info.had_effect;
	return result;
}

State* amaya_mtbdd_get_transition_target(
        MTBDD 		mtbdd, 
        uint8_t* 	cube, 
        uint32_t 	cube_size, 
        uint32_t* 	result_size) 
{

	auto search_result_tds = _get_transition_target(mtbdd, 1, cube, cube_size);
	if (search_result_tds == NULL)
	{
		*result_size = 0;
		return NULL;
	} else {
		const uint32_t rs = search_result_tds->destination_set->size();
		*result_size = rs;
		auto result_arr = (State*) malloc(sizeof(uint32_t) * rs);
		uint32_t i = 0;
		for (auto state : *search_result_tds->destination_set) 
		{
			result_arr[i++] = state;
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
		if (*variable_assigments == 0) {
			// Go low
			return _get_transition_target(
					mtbdd_getlow(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);
		} 
		else if (*variable_assigments == 1) {
			// Go high
			return _get_transition_target(
					mtbdd_gethigh(root), 
					current_variable + 1, 
					variable_assigments + 1, 
					var_count - 1);
		}
		else {
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

	} else {
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
		MTBDD* 		mtbdd_roots, 
		uint32_t 	root_count,
		State* 		names, // [(old, new), (old, new), (old, new)] 
		uint32_t 	name_count)
{
	State old_state_name, new_state_name;

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

		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		auto new_leaf_contents = new std::set<State>();

		for (uint32_t i = 0; i < name_count; i++) {

			old_state_name = names[2*i];
			new_state_name = names[2*i + 1];

			bool is_in_leaf = (tds->destination_set->find(old_state_name) != tds->destination_set->end());

			if (is_in_leaf) {
				new_leaf_contents->insert(new_state_name);
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

	free(arr);
}

State* amaya_mtbdd_get_leaves(
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
	for (MTBDD leaf : leaves) {
		Transition_Destination_Set* tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		size_cnt += tds->destination_set->size();
	}

	// Do the allocations
	auto _leaf_sizes = (uint32_t*) malloc(sizeof(uint32_t) * leaves.size());  // One slot for each leaf 
	auto _leaf_states = (State*) malloc(sizeof(State) * size_cnt);
	
	// Populate
	uint32_t state_i = 0;
	uint32_t size_i = 0;

	for (MTBDD leaf : leaves) {
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (int state : *tds->destination_set) {
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

void amaya_replace_leaf_contents_with(void *leaf_tds, State* new_contents, uint32_t contents_size)
{
    auto tds = (Transition_Destination_Set*) leaf_tds;
    tds->destination_set->clear();
    for (uint32_t i = 0; i < contents_size; i++) {
        tds->destination_set->insert(new_contents[i]);
    }
}

State* amaya_rename_metastates_to_int(
		MTBDD* 		roots, 							// MTBDDs resulting from determinization
		uint32_t 	root_cnt,						// Root count
		State 		start_numbering_metastates_from,
		uint32_t 	resulting_automaton_id,
		uint32_t**	out_metastates_sizes,
		uint32_t*	out_metastates_cnt)
{

	set<MTBDD> leaves;
	for (uint32_t root_i = 0; root_i < root_cnt; root_i++) {
		collect_mtbdd_leaves(roots[root_i], leaves);
	}

	uint32_t metastates_cnt = leaves.size();
	uint32_t metastates_sizes_first_empty_i = 0;
	auto metastates_sizes = (uint32_t *) malloc(metastates_cnt * sizeof(uint32_t));
	assert(metastates_sizes != NULL);

	vector<State> metastates_serialized;

	for (auto leaf: leaves) {
		auto leaf_contents = (Transition_Destination_Set*) mtbdd_getvalue(leaf);

		// Serialize the mapping so that the Python side will know what was mapped to what.
		for (auto state : *leaf_contents->destination_set) {
			metastates_serialized.push_back(state);
		}

		metastates_sizes[metastates_sizes_first_empty_i++] = leaf_contents->destination_set->size();
		
		// We need to convert the internal transition destination set into a sorted vector.
		// @Warn: Currently we rely on the fact that we use std::set that stores the destination states
		// 		  in a sorted manner (a balanced binary tree)
		
		// Do the actual renaming.
		int metastate_num = start_numbering_metastates_from++;
		leaf_contents->automaton_id = resulting_automaton_id;
		leaf_contents->destination_set->clear();
		leaf_contents->destination_set->insert(metastate_num);
	} 

	*out_metastates_sizes = metastates_sizes;
	*out_metastates_cnt = leaves.size();
	
	auto metastates_serialized_arr = (State*) malloc(metastates_serialized.size() * sizeof(State));
	assert(metastates_serialized_arr != NULL);

	for (uint32_t i = 0; i < metastates_serialized.size(); i++) {
		metastates_serialized_arr[i] = metastates_serialized.at(i); 
	}

	return metastates_serialized_arr;
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

State* amaya_mtbdd_get_state_post(MTBDD dd, uint32_t *post_size)
{
	set<MTBDD> leaves;
	collect_mtbdd_leaves(dd, leaves);
	
	set<State> state_post_set;

	for (auto leaf: leaves)	{
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		state_post_set.insert(
                tds->destination_set->begin(),
                tds->destination_set->end());
	}
	
	auto result = (State*) malloc(sizeof(State) * state_post_set.size());
	uint32_t i = 0;

	for (int state: state_post_set)	{
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

	    auto pci = (Pad_Closure_Info*) op_param;

		// Check whether the transition destination state even leads the the right state (might not)
		if (left_tds->destination_set->find(pci->right_state) == left_tds->destination_set->end()) {
			// Does not lead to the right state, therefore cannot propagate padding.
			return left; 
		}

	    bool is_final;
	    for (State rs: *right_tds->destination_set)
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


bool amaya_mtbdd_do_pad_closure(
		State 		left_state, 
		MTBDD 		left_dd,
		State 		right_state, 
		MTBDD 		right_dd, 
		State* 		final_states, 
		uint32_t 	final_states_cnt)
{
	Pad_Closure_Info pci = {0};
	pci.final_states = final_states;
	pci.final_states_cnt = final_states_cnt;
	pci.had_effect = false;

	pci.right_state = right_state;
	pci.left_state = left_state;
	

	LACE_ME;
	//sylvan_clear_cache();
	mtbdd_applyp(left_dd, 
				 right_dd,
				 (uint64_t) &pci, 
				 TASK(pad_closure_op), 
				 CUR_PADDING_CLOSURE_ID);

	// Every padding closure that is performed (even repeated) should be treated
	// as a unique operation --> so that the caching problems will not occur
	CUR_PADDING_CLOSURE_ID += 4;

	return pci.had_effect;
}

uint8_t* amaya_mtbdd_get_transitions(
		MTBDD 		root, 	// The mtbdd
		uint32_t* 	vars,	
		uint32_t 	var_count,
		uint32_t* 	symbols_cnt,
		State** 	dest_states,
		uint32_t** 	dest_states_sizes)
{
	LACE_ME;
	
	// Construct MTBDD set from given vars
	MTBDD variable_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < var_count; i++) {
		variable_set= mtbdd_set_add(variable_set, vars[i]);
	}

	// Stores the tree path to the leaf
	auto arr = (uint8_t*) malloc(sizeof(uint8_t) * var_count);

	MTBDD leaf = mtbdd_enum_first(root, variable_set, arr, NULL);
	
	vector<State> dest_states_vec;
	vector<uint8_t> symbols;
	vector<uint32_t> state_sizes;
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
	auto symbols_out = (uint8_t*) malloc(sizeof(uint8_t) * symbols.size());
	for (uint32_t i = 0; i < symbols.size(); i++) symbols_out[i] = symbols.at(i);

	auto dest_states_sizes_out = (uint32_t*) malloc(sizeof(uint32_t) * state_sizes.size());
	for (uint32_t i = 0; i < state_sizes.size(); i++) dest_states_sizes_out[i] = state_sizes.at(i);

	auto  dest_states_out = (State*) malloc(sizeof(State) * dest_states_vec.size());
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
	set<MTBDD> leaves {};
    for (uint32_t i = 0; i < root_cnt; i++) {
        collect_mtbdd_leaves(roots[i], leaves);
    }

    for (auto leaf : leaves) {
        auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf); 
        tds->automaton_id = new_id;
    }
}

MTBDD amaya_mtbdd_intersection(
        MTBDD 		a, 
        MTBDD 		b,
        uint32_t 	result_automaton_id, 
        State** 	discovered_states,         // OUT
        uint32_t*  	discovered_states_cnt)     // OUT
{
    Intersection_Op_Info intersect_info = {0};
    intersect_info.automaton_id = result_automaton_id;
    intersect_info.discoveries = new vector<State>();

	LACE_ME;
	MTBDD result = mtbdd_applyp(a, b, (uint64_t) &intersect_info, TASK(set_intersection_op), AMAYA_INTERSECTION_OP_ID);
    
    
    if (!intersect_info.discoveries->empty()) {
        auto discovered_states_arr = (State*) malloc(sizeof(State) * intersect_info.discoveries->size());
        for (uint32_t i = 0; i < intersect_info.discoveries->size(); i++) {
            discovered_states_arr[i] = intersect_info.discoveries->at(i);
        }

        // Send the new discoveries to the python side.
        *discovered_states = discovered_states_arr;
        *discovered_states_cnt = intersect_info.discoveries->size() / 3;
    } else {
        // Nothing new was discovered, do not malloc
        *discovered_states = NULL;
        *discovered_states_cnt = 0;
    }

    delete intersect_info.discoveries;

    return result;
}

MTBDD amaya_unite_mtbdds(MTBDD m1, MTBDD m2, uint32_t automaton_id) {
	LACE_ME;
	MTBDD u = mtbdd_applyp(m1, m2, (uint64_t) automaton_id, TASK(set_union), AMAYA_UNION_OP_ID);
	return u;
}

void amaya_begin_intersection(
		bool prune_state_pairs_with_one_final, 
		State* prune_final_states, 
		uint32_t final_states_cnt) 
{
    intersection_state = (Intersection_State*) malloc(sizeof(Intersection_State));
	intersection_state->intersection_state_pairs_numbers = new map<pair<State, State>, State>();

	if (prune_state_pairs_with_one_final) {
		intersection_state->should_do_early_prunining = true;
		intersection_state->prune_final_states = new unordered_set<State>();
		for (uint32_t i = 0; i < final_states_cnt; i++) {
			intersection_state->prune_final_states->insert(prune_final_states[i]);
		}
	} else {
		intersection_state->should_do_early_prunining = false;
		intersection_state->prune_final_states = NULL;
	}
}

void amaya_update_intersection_state(State* metastates, State* renamed_metastates, uint32_t cnt)
{
    for (uint32_t i = 0; i < cnt; i++) {
        const auto metastate = std::make_pair(metastates[2*i], metastates[2*i + 1]);
        intersection_state
			->intersection_state_pairs_numbers
			->insert(std::make_pair(metastate, renamed_metastates[i])); 
    }
}

void amaya_end_intersection() 
{
	delete intersection_state->intersection_state_pairs_numbers;
	if (intersection_state->should_do_early_prunining) {
		delete intersection_state->prune_final_states;
	}
    free(intersection_state);
    intersection_state = NULL;
}

void amaya_set_debugging(bool debug) {
	DEBUG_ON = debug;
}

State* amaya_get_state_post_with_some_transition(
		MTBDD mtbdd,
		uint32_t* variables,
		uint32_t variable_cnt,
		uint8_t** out_symbols,
		uint32_t* transition_cnt)
{
	LACE_ME;

	// Create a MTBDD Cube containing all requred variables.
	MTBDD variable_set = mtbdd_set_empty();
	for (uint32_t i = 0; i < variable_cnt; i++) {
		variable_set = mtbdd_set_add(variable_set, variables[i]);
	}
		
	// Stores the tree path to the leaf
	uint8_t* arr = (uint8_t*) malloc(sizeof(uint8_t) * variable_cnt);

	MTBDD leaf = mtbdd_enum_first(mtbdd, variable_set, arr, NULL);

	vector<State> 	reachable_states;
	vector<uint8_t> transition_symbols;

 	while (leaf != mtbdd_false) {
		auto tds = (Transition_Destination_Set*) mtbdd_getvalue(leaf);
		for (auto state : *tds->destination_set) {

			// @Optimize: We use linear search over vector - there should be a relatively small number of states 
			// 			  reachable from everystate, therefore it should be faster to use ordinary vector (chache lines)
			bool state_already_located = false;
			for (auto already_located_state : reachable_states) {
				if (already_located_state == state) {
					state_already_located = true;
					break;
				}
			}

			if (!state_already_located) {
				reachable_states.push_back(state);
				// Copy the transition symbol so that the Python side will have the information available 
				// (this method is used in DFS) 
				
				for (uint32_t i = 0; i < variable_cnt; i++) {
					transition_symbols.push_back(arr[i]);
				}   	
			}
		}

 	    leaf = mtbdd_enum_next(mtbdd, variable_set, arr, NULL);
 	}

	free(arr);

	auto _out_symbols = (uint8_t*) malloc(sizeof(uint8_t) * transition_symbols.size());
	for (uint32_t i = 0; i < transition_symbols.size(); i++) {
		_out_symbols[i] = transition_symbols.at(i);
	} 

	*out_symbols = _out_symbols;
	*transition_cnt = (uint32_t) reachable_states.size();
	
	auto _reachable_states = (State*) malloc(sizeof(State) * reachable_states.size());
	for (uint32_t i = 0; i < reachable_states.size(); i++) _reachable_states[i] = reachable_states.at(i); 

	return _reachable_states;
}


MTBDD* amaya_remove_states_from_transitions(
	MTBDD* 		transition_roots,
	uint32_t 	transition_cnt,
	State* 		states_to_remove,
	uint32_t 	states_to_remove_cnt)
{

	auto mtbdds_after_removal = (MTBDD*) malloc(sizeof(MTBDD) * transition_cnt);
	LACE_ME;

	set<State> states_to_remove_set;
	for (uint32_t i = 0; i < states_to_remove_cnt; i++) {
		states_to_remove_set.insert(states_to_remove[i]);
	}

	REMOVE_STATES_OP_PARAM = &states_to_remove_set;
	for (uint32_t i = 0; i < transition_cnt; i++) {
		MTBDD result_mtbdd = mtbdd_uapply(transition_roots[i], TASK(remove_states_op), REMOVE_STATES_OP_COUNTER);
		mtbdds_after_removal[i] = result_mtbdd;
	}

	REMOVE_STATES_OP_COUNTER++;

	return mtbdds_after_removal;
}
