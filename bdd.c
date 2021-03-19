#include <assert.h>
#include <stdio.h>
#include <sylvan.h>
#include <sylvan_bdd.h>
#include <sylvan_common.h>
#include <sylvan_mtbdd.h>
#include <sylvan_stats.h>

VOID_TASK_4(_my_callback, void *, ctx, BDDVAR *, var, uint8_t *, arr, int,
            arr_length) {
  (void)ctx;
}

TASK_DECL_2(MTBDD, my_op, MTBDD *, MTBDD *);
TASK_IMPL_2(MTBDD, my_op, MTBDD *, pa, MTBDD *, pb) {
  MTBDD a = *pa, b = *pb;

  /* Check for partial functions */
  if (a == mtbdd_false)
    return b;
  if (b == mtbdd_false)
    return a;

  if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
    printf("Both are leaves with values: %ld %ld\n", mtbdd_getvalue(a),
           mtbdd_getvalue(b));
    return mtbdd_int64(mtbdd_getvalue(a) + mtbdd_getvalue(b)); // Do not compute
  }

  /* Commutative, so swap a,b for better cache performance */
  if (a < b) {
    *pa = b;
    *pb = a;
  }

  return mtbdd_invalid;
}

VOID_TASK_0(simple_fn) {
  BDD b1 = sylvan_ithvar(1);
  BDD b2 = sylvan_ithvar(2);

  MTBDD leaf = mtbdd_int64(42);

  BDDSET vars = mtbdd_set_empty();
  vars = mtbdd_set_add(vars, 2);
  vars = mtbdd_set_add(vars, 1);

  uint8_t cube_vars[] = {1, 1};

  MTBDD my_mtbdd = mtbdd_cube(vars, cube_vars, leaf);
  MTBDD h1 = mtbdd_gethigh(my_mtbdd);
  MTBDD p_leaf = mtbdd_gethigh(h1);
  printf("Is leaf after 2 highs?: %s\n",
         mtbdd_isleaf(p_leaf) ? "TRUE" : "FALSE");
  if (mtbdd_isleaf(p_leaf)) {
    printf("Leaf value: %ld\n", mtbdd_getvalue(p_leaf));
  }

  uint8_t cube_vars2[] = {1, 1};
  MTBDD my_mtbdd2 = mtbdd_cube(vars, cube_vars2, leaf);

  MTBDD mtbdd_apply = mtbdd_apply(my_mtbdd, my_mtbdd2, TASK(my_op));

  MTBDD tn;
  tn = mtbdd_gethigh(mtbdd_apply);
  assert(!mtbdd_isleaf(tn));

  tn = mtbdd_gethigh(tn);
  assert(mtbdd_isleaf(tn));

  if (mtbdd_isleaf(tn)) {
    printf("Apply plus value: %ld\n", mtbdd_getvalue(tn));
  }

  int ctx = 0;

  puts("Done");
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
