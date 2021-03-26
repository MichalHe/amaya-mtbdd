#!/bin/python3
import ctypes

h = ctypes.CDLL('./amaya-mtbdd.so', mode=1)
h.init_machinery()
mtbdd_false = ctypes.c_ulong.in_dll(h, 'w_mtbdd_false')
mtbdd_true = ctypes.c_ulong.in_dll(h, 'w_mtbdd_true')

v = ctypes.c_uint32(2)


def write_mtbdd_dot_to_file(m, filename):
    '''Writes the dot representation of given MTBDD m to the file.'''
    output_f = open(filename, 'w')
    fd = output_f.fileno()

    h.amaya_print_dot(
        m,  # The mtbdd to be printed out
        fd  # File descriptor
    )

    output_f.close()  # The C side does not close the file descriptor


# Fil in the transitions:
t = [(1, 1, 1, 1)]
transitions = (ctypes.c_uint8 * 4 * len(t))()
for i in range(len(t)):
    for j in range(4):
        transitions[i][j] = t[i][j]


t_count = ctypes.c_uint32(len(t))
t_size = ctypes.c_uint32(4)

destination_set = (ctypes.c_uint32 * 2)(2, 3)


# Construct two MTBDDs
x = h.amaya_mtbdd_build_single_terminal(
    transitions,        # Transition symbols, 2d array
    t_count,            # Transitions symbols count
    t_size,             # Variables count
    destination_set,    # Set of the terminal states
    ctypes.c_uint32(2)  # Destination set size
)

write_mtbdd_dot_to_file(x, '/tmp/amaya_mtbdd_first.dot')

# Second mtbdd
destination_set = (ctypes.c_uint32 * 2)(1, 3)
transitions[0][2] = 0  # Transition via symbol (1, 1, 0, 0) to {1, 3}
transitions[0][3] = 0
y = h.amaya_mtbdd_build_single_terminal(
    transitions,        # Transition symbols, 2d array
    t_count,            # Transitions symbols count
    t_size,             # Variables count
    destination_set,    # Set of the terminal states
    ctypes.c_uint32(2)  # Destination set size
)

write_mtbdd_dot_to_file(y, '/tmp/amaya_mtbdd_second.dot')

# Conceptually the two created MTBDDs can be viewed as the result of the
# automaton building procedure that converts inequation to automaton.
# Say that we are in a state (1) and we know that via symbols (A, B) we can get
# to (1, 2), and via symbol (C) to (3). Then we build two mtbdds:
# (1) ---- (A, B) ---> (1, 2)
# (1) ---- (C) ---> (3)
# The union of these mtbdds is the transition fn for the state (1)

u = h.amaya_unite_mtbdds(x, y)

write_mtbdd_dot_to_file(u, '/tmp/amaya_mtbdd_after_union.dot')

# We would like to project the 3rd and 4th variables away (to force terminal
# merge)
variables = (ctypes.c_uint32 * 2)(3, 4)
f = h.amaya_project_variables_away(
   u,                  # The MTBDD
   variables,          # The variables to be projected away
   ctypes.c_uint32(2)  # The number of variables in the variables array
)

write_mtbdd_dot_to_file(f, '/tmp/amaya_mtbdd_after_projection.dot')

h.shutdown_machinery()  # Cleanup
