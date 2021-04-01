import ctypes as ct
from typing import Dict, Any, Tuple, Union

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


mtbdd_false = ct.c_ulong.in_dll(mtbdd_wrapper, 'w_mtbdd_false')
MTBDD = ct.c_ulong
Symbol = Tuple[Union[str, int], ...]


class MTBDDTransitionFn():
    def __init__(self):
        self.mtbdds: Dict[Any, MTBDD] = {}

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


M = MTBDDTransitionFn()

M.insert_transition(0, (0, 0, 1), 1)
M.insert_transition(0, (0, 0, 1), 2)
M.insert_transition(0, (0, 1, 1), 3)
M.insert_transition(0, (0, 1, 1), 4)
M.rename_states({1: 10, 2: 20, 3: 30, 4: 40})
M.write_mtbdd_dot_to_file(M.mtbdds[0], '/tmp/AYAY.dot')
print(M.get_transition_target(0, (0, 0, 1)))
print(M.get_transition_target(0, (0, 1, 1)))
