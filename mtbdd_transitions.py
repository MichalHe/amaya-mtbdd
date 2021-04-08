import ctypes as ct
from typing import Dict, Any, Tuple, Union, List

mtbdd_wrapper = ct.CDLL('./amaya-mtbdd.so', mode=1)
mtbdd_wrapper.init_machinery()
mtbdd_wrapper.amaya_mtbdd_get_transition_target.argtypes = (
    ct.c_ulong,
    ct.POINTER(ct.c_uint8),
    ct.c_uint32,
    ct.POINTER(ct.c_uint32)
)
mtbdd_wrapper.amaya_mtbdd_get_transition_target.restype = ct.POINTER(ct.c_int)
mtbdd_wrapper.amaya_mtbdd_rename_states.argtypes = (
    ct.c_long,
    ct.POINTER(ct.c_int),
    ct.c_uint32,
)

mtbdd_wrapper.amaya_project_variables_away.argtypes = (
    ct.c_ulong,
    ct.POINTER(ct.c_uint32),
    ct.c_uint32
)

mtbdd_wrapper.amaya_mtbdd_get_leaves.argtypes = (
    ct.c_ulong,
    ct.POINTER(ct.POINTER((ct.c_uint32))),  # Leaf sizes
    ct.POINTER(ct.c_uint32)  # Leaf count
)
mtbdd_wrapper.amaya_mtbdd_get_leaves.restype = ct.POINTER(ct.c_int)  # Pointer to array containing the states

mtbdd_wrapper.amaya_unite_mtbdds.argtypes = (
    ct.c_ulong,  # MTBDD a
    ct.c_ulong,  # MTBDD b
)
mtbdd_wrapper.amaya_unite_mtbdds.restype = ct.c_ulong

mtbdd_wrapper.amaya_mtbdd_get_state_post.argtypes = (
    ct.c_ulong,  # MTBDD m
    ct.POINTER(ct.c_uint32),  # MTBDD result array size
)
mtbdd_wrapper.amaya_mtbdd_get_state_post.restype = ct.POINTER(ct.c_int)

mtbdd_wrapper.amaya_mtbdd_do_pad_closure.argtypes = (
    ct.c_ulong,  # MTBDD left
    ct.c_ulong,  # MTBDD right
    ct.POINTER(ct.c_int),  # Array with final states.
    ct.c_uint32  # Number of final states
)
mtbdd_wrapper.amaya_mtbdd_do_pad_closure.restype = ct.c_bool

mtbdd_wrapper.amaya_mtbdd_get_transitions.argtypes = (
    ct.c_ulong,                           # MTBDD root
    ct.POINTER(ct.c_uint32),              # Array with variables.
    ct.c_uint32,                          # Size of array with variables
    ct.POINTER(ct.c_uint32),              # OUT, number of transitions
    ct.POINTER(ct.POINTER(ct.c_int)),     # OUT Pointer to array containing destinations states
    ct.POINTER(ct.POINTER(ct.c_uint32)),  # OUT Pointer to array containing sizes of serialized destinations states
)
mtbdd_wrapper.amaya_mtbdd_get_transitions.restype = ct.POINTER(ct.c_uint8)


mtbdd_false = ct.c_ulong.in_dll(mtbdd_wrapper, 'w_mtbdd_false')
MTBDD = ct.c_ulong
Symbol = Tuple[Union[str, int], ...]


class MTBDDTransitionFn():
    def __init__(self):
        self.mtbdds: Dict[Any, MTBDD] = {}
        self.post_cache = dict()
        self.is_post_cache_valid = False

    def insert_transition(self,
                          source: Any,
                          symbol: Symbol,
                          dest: int):
        assert type(dest) == int
        assert type(source) == int

        if source not in self.mtbdds:
            current_mtbdd = mtbdd_false
            self.mtbdds[source] = current_mtbdd
        else:
            current_mtbdd = self.mtbdds[source]

        cube_size, cube = self.convert_symbol_to_mtbdd_cube(symbol)

        # Destination:
        dest_state = (ct.c_uint32 * 1)()
        dest_state[0] = dest

        # Construct cube from destination state
        mtbdd_new = mtbdd_wrapper.amaya_mtbdd_build_single_terminal(
            cube,            # Transition symbols, 2d array
            ct.c_uint32(1),  # Transitions symbols count
            cube_size,       # Variables count
            dest_state,      # Set of the terminal states
            ct.c_uint32(1)   # Destination set size
        )

        resulting_mtbdd = mtbdd_wrapper.amaya_unite_mtbdds(mtbdd_new, current_mtbdd)
        self.mtbdds[source] = resulting_mtbdd

    def write_mtbdd_dot_to_file(self, m, filename):
        '''Writes the dot representation of given MTBDD m to the file.'''
        output_f = open(filename, 'w')
        fd = output_f.fileno()

        mtbdd_wrapper.amaya_print_dot(
            m,  # The mtbdd to be printed out
            fd  # File descriptor
        )

        output_f.close()  # The C side does not close the file descriptor

    def convert_symbol_to_mtbdd_cube(self, symbol: Symbol):
        cube = (ct.c_uint8 * len(symbol) * 1)()
        for i, bit in enumerate(symbol):
            if bit == '*':
                cube[0][i] = ct.c_uint8(2)
            else:
                cube[0][i] = ct.c_uint8(int(bit))

        cube_size = ct.c_uint32(len(symbol))

        return cube_size, cube

    def get_transition_target(self, source, symbol):
        '''Retrieve the set of states that lead from `source` via `symbol`.'''
        if source not in self.mtbdds:
            return set()
        mtbdd = self.mtbdds[source]

        cs, c = self.convert_symbol_to_mtbdd_cube(symbol)
        result_size = ct.c_uint32(0)
        result = mtbdd_wrapper.amaya_mtbdd_get_transition_target(
            mtbdd,
            ct.cast(c, ct.POINTER(ct.c_uint8)),      # Cube
            cs,     # Cube size
            ct.byref(result_size)
        )

        dest_ = []
        for i in range(int(result_size.value)):
            dest_.append(result[i])

        mtbdd_wrapper.amaya_do_free(result)
        return set(dest_)

    def rename_states(self, mappings: Dict[int, int]):
        '''Renames all states referenced within stored mtbdds with the
        provided mapping.

        Requires all referenced states to be present in the mapping.
        '''
        flat_mapping_size = 2*len(mappings)
        arr = (ct.c_int * flat_mapping_size)()

        for i, mapping in enumerate(mappings.items()):
            old, new = mapping
            arr[2*i] = ct.c_int(old)
            arr[2*i + 1] = ct.c_int(new)

        mapping_ptr = ct.cast(arr, ct.POINTER(ct.c_int))
        mapping_size = ct.c_uint32(len(mappings))

        for mtbdd in self.mtbdds.values():
            mtbdd_wrapper.amaya_mtbdd_rename_states(
                mtbdd,
                mapping_ptr,
                mapping_size
            )

    def project_variable_away(self, variable: int):
        '''Not sure what happens when trying to project a variable that is not
        present.'''
        assert variable > 0, 'MTBDD variables are numbered via ints from 1 up'

        variables = (ct.c_uint32 * 1)(variable)
        for state in self.mtbdds:
            mtbdd = self.mtbdds[state]
            new_mtbdd = mtbdd_wrapper.amaya_project_variables_away(
                mtbdd,
                variables,
                ct.c_uint32(1)
            )
            self.mtbdds[state] = new_mtbdd

    def get_mtbdd_leaves(self, mtbdd: MTBDD) -> List[List[int]]:
        ''' Internal procedure. '''
        leaf_sizes = ct.POINTER(ct.c_uint32)()
        leaf_cnt = ct.c_uint32()
        leaf_contents = mtbdd_wrapper.amaya_mtbdd_get_leaves(
            mtbdd,
            ct.byref(leaf_sizes),
            ct.byref(leaf_cnt)
        )

        state_i = 0
        leaves = []
        for i in range(leaf_cnt.value):
            leaf = []
            for j in range(leaf_sizes[i]):
                state = leaf_contents[state_i]
                state_i += 1
                leaf.append(state)
            leaves.append(leaf)

        mtbdd_wrapper.amaya_do_free(leaf_sizes)
        mtbdd_wrapper.amaya_do_free(leaf_contents)

        return leaves

    def get_union_mtbdd_for_states(self, states: List[int]) -> MTBDD:
        '''Does what name suggests.'''
        resulting_mtbdd = mtbdd_false
        for state in states:
            if state not in self.mtbdds:
                continue
            mtbdd = self.mtbdds[state]
            resulting_mtbdd = mtbdd_wrapper.amaya_unite_mtbdds(
                mtbdd,
                resulting_mtbdd)
        return resulting_mtbdd

    def get_state_post(self, state: int) -> List[int]:
        mtbdd = self.mtbdds.get(state, None)
        if mtbdd is None:
            return []
        result_size = ct.c_uint32()
        state_post_arr = mtbdd_wrapper.amaya_mtbdd_get_state_post(mtbdd, ct.byref(result_size))

        result = []
        for i in range(result_size.value):
            result.append(state_post_arr[i])
        return result

    def get_state_pre(self, state: int, initial_state: int) -> List[int]:
        if not self.is_post_cache_valid:
            self.post_cache = self._build_morph_map(initial_state)

        # Use the post cache to calculate state pre.
        state_pre = set()
        for s in self.post_cache:
            s_post = self.post_cache[s]
            if state in s_post:
                state_pre.add(s)
        return list(state_pre)

    def _build_morph_map(self, initial_state: int) -> Dict[int, List[int]]:
        '''Builds the image of the transition function with the information
        about transition symbols left out.'''

        morph_map = {}
        work_queue = [initial_state]
        while work_queue:
            state = work_queue.pop(-1)
            state_post = self.get_state_post(state)
            if not state_post:
                # state is terminal
                continue
            morph_map[state] = state_post
            for new_state in state_post:
                # Check whether the state is not scheduled/processed already
                if new_state not in work_queue and new_state not in morph_map:
                    work_queue.append(new_state)
        return morph_map

    def do_pad_closure(self, initial_state: int, final_states: List[int]):
        '''Performs padding closure on the underlying automatic structure.'''
        # Initialize the working queue with all states, that have some final
        # state in their Post
        final_states_pre_set = set()
        for fs in final_states:
            final_states_pre_set.update(self.get_state_pre(fs, initial_state=initial_state))
        work_queue = list(final_states_pre_set)

        while work_queue:
            state = work_queue.pop()
            state_pre_list = self.get_state_pre(state, initial_state=initial_state)
            for pre_state in state_pre_list:
                print(f'Applying PC on: {pre_state} ---> {state}')
                had_pc_effect = self._do_pad_closure_single(pre_state, state, final_states)
                if had_pc_effect:
                    if pre_state not in work_queue:
                        work_queue.append(pre_state)

    def _do_pad_closure_single(self, left_state: int, right_state: int, final_states: List[int]) -> bool:
        '''(left_state) --A--> (right_state) --B--> final_states'''
        left_mtbdd = self.mtbdds.get(left_state, None)
        right_mtbdd = self.mtbdds.get(right_state, None)

        if left_mtbdd is None or right_mtbdd is None:
            return  # We need a pair of mtbdds to perform padding closure.

        # Convert the set into C array
        final_states_arr = (ct.c_int * len(final_states))(*list(final_states))
        final_states_cnt = ct.c_uint32(len(final_states))

        was_modified = mtbdd_wrapper.amaya_mtbdd_do_pad_closure(
            left_mtbdd,
            right_mtbdd,
            final_states_arr,
            final_states_cnt
        )

        return bool(was_modified)

    def iter_transitions(self, state: int, variables: List[int]):
        mtbdd = self.mtbdds.get(state, None)
        if mtbdd is None:
            return

        _vars = (ct.c_uint32 * len(variables))(*variables)

        transition_dest_states = ct.POINTER(ct.c_int)()
        transition_dest_states_sizes = ct.POINTER(ct.c_uint32)()
        transition_count = ct.c_uint32()

        symbols = mtbdd_wrapper.amaya_mtbdd_get_transitions(
            mtbdd,
            ct.cast(_vars, ct.POINTER(ct.c_uint32)),
            ct.c_uint32(len(variables)),
            ct.byref(transition_count),
            ct.byref(transition_dest_states),
            ct.byref(transition_dest_states_sizes)
        )

        s_length = len(variables)
        i = 0
        # ti = transition index
        for ti in range(transition_count.value):
            symbol = []
            for s in range(s_length):
                symbol.append(symbols[s_length*ti + s])

            dest = []
            dest_size = transition_dest_states_sizes[ti]
            for _ in range(dest_size):
                dest.append(transition_dest_states[i])
                i += 1

            for dest_state in dest:
                yield (state, symbol, dest_state)

        mtbdd_wrapper.amaya_do_free(transition_dest_states)
        mtbdd_wrapper.amaya_do_free(transition_dest_states_sizes)
        mtbdd_wrapper.amaya_do_free(symbols)


def determinize_mtbdd(tfn: MTBDDTransitionFn, initial_states: List[int]):
    work_queue = [tuple(initial_states)]
    states = set()
    while work_queue:
        c_metastate = work_queue.pop(-1)
        states.add(c_metastate)
        print(c_metastate)
        transition = tfn.get_union_mtbdd_for_states(c_metastate)

        reachable_states = tfn.get_mtbdd_leaves(transition)
        for rs in reachable_states:
            rs = tuple(rs)
            if rs not in states:
                # @Optimize: this is a linear search.
                if rs not in work_queue:
                    work_queue.append(rs)
    return states


if __name__ == '__main__':
    tfn = MTBDDTransitionFn()
    zeta0 = (0, 0, 1)
    zeta1 = (1, 0, 0)

    tfn.insert_transition(0, zeta0, 1)
    tfn.insert_transition(0, zeta1, 1)

    tfn.insert_transition(1, zeta1, 2)
    tfn.insert_transition(1, zeta0, 3)

    tfn.insert_transition(2, zeta1, 3)

    final_states = [3]
    initial_state = 0

    # Initialstate, final states
    # tfn.do_pad_closure(initial_state, final_states)
    print(list(tfn.iter_transitions(0, [1, 2, 3])))
