#include "../include/tfa_leaf.h"
#include "../include/base.hpp"
#include "../include/operations.hpp"

#include <inttypes.h>

#include <string>
#include <cstring>
#include <utility>
#include <stdlib.h>
#include <sstream>
#include <sylvan.h>

using namespace sylvan;


u64 mtbdd_tfa_leaf_type_id;
u64 mtbdd_tfa_intersection_leaf_type_id;
u64 mtbdd_tfa_pareto_leaf_type_id;

void init_tfa_leaf() {
    mtbdd_tfa_leaf_type_id = sylvan::sylvan_mt_create_type();
    sylvan::sylvan_mt_set_hash(mtbdd_tfa_leaf_type_id, tfa_leaf_hash);
    sylvan::sylvan_mt_set_equals(mtbdd_tfa_leaf_type_id, tfa_leaf_equals);
    sylvan::sylvan_mt_set_create(mtbdd_tfa_leaf_type_id, mk_tfa_leaf_from_contents);
    sylvan::sylvan_mt_set_destroy(mtbdd_tfa_leaf_type_id, destroy_tfa_leaf);
    sylvan::sylvan_mt_set_to_str(mtbdd_tfa_leaf_type_id, tfa_leaf_to_str);
}

void init_tfa_intersection_leaf() {
    mtbdd_tfa_intersection_leaf_type_id = sylvan::sylvan_mt_create_type();
    sylvan::sylvan_mt_set_hash(mtbdd_tfa_intersection_leaf_type_id, tfa_intersection_leaf_hash);
    sylvan::sylvan_mt_set_equals(mtbdd_tfa_intersection_leaf_type_id, tfa_intersection_leaf_equals);
    sylvan::sylvan_mt_set_create(mtbdd_tfa_intersection_leaf_type_id, mk_tfa_intersection_leaf_from_contents);
    sylvan::sylvan_mt_set_destroy(mtbdd_tfa_intersection_leaf_type_id, destroy_tfa_intersection_leaf);
    sylvan::sylvan_mt_set_to_str(mtbdd_tfa_intersection_leaf_type_id, tfa_leaf_intersection_to_str);
}

void init_tfa_pareto_leaf() {
    mtbdd_tfa_pareto_leaf_type_id = sylvan::sylvan_mt_create_type();
    sylvan::sylvan_mt_set_hash(mtbdd_tfa_pareto_leaf_type_id, tfa_pareto_leaf_hash);
    sylvan::sylvan_mt_set_equals(mtbdd_tfa_pareto_leaf_type_id, tfa_pareto_leaf_equals);
    sylvan::sylvan_mt_set_create(mtbdd_tfa_pareto_leaf_type_id, mk_tfa_pareto_leaf_from_contents);
    sylvan::sylvan_mt_set_destroy(mtbdd_tfa_pareto_leaf_type_id, destroy_tfa_pareto_leaf);
    sylvan::sylvan_mt_set_to_str(mtbdd_tfa_pareto_leaf_type_id, tfa_pareto_leaf_to_str);
}

void init_tfa_leaves() {
   init_tfa_leaf(); 
   init_tfa_intersection_leaf(); 
   init_tfa_pareto_leaf(); 
}


sylvan::MTBDD make_tfa_leaf(TFA_Leaf_Contents* contents) {
    sylvan::MTBDD leaf = sylvan::mtbdd_makeleaf(mtbdd_tfa_leaf_type_id, (uint64_t) contents);
    return leaf;
}

template <typename Leaf>
void mk_leaf_from_contents(u64* contents_ptr) {
    auto leaf_contents_ptr = reinterpret_cast<Leaf**>(contents_ptr);
    Leaf* contents = *leaf_contents_ptr;
    Leaf* contents_copy = new Leaf(*contents);
    *contents_ptr = reinterpret_cast<u64>(contents_copy);
}

template <typename Leaf>
void destroy_leaf(u64 leaf_value) {
    auto leaf_contents_ptr = reinterpret_cast<Leaf*>(leaf_value);
    delete leaf_contents_ptr;
}

template <typename Leaf>
int leaf_equals(uint64_t a_ptr, uint64_t b_ptr) {
    auto left_leaf_contents  = reinterpret_cast<Leaf*>(a_ptr);
    auto right_leaf_contents = reinterpret_cast<Leaf*>(b_ptr);
    return (left_leaf_contents == right_leaf_contents);
}

void mk_tfa_leaf_from_contents(u64* contents_ptr) {
    mk_leaf_from_contents<TFA_Leaf_Contents>(contents_ptr);
}

void destroy_tfa_leaf(u64 leaf_value) {
    destroy_leaf<TFA_Leaf_Contents>(leaf_value);
}

int tfa_leaf_equals(uint64_t a_ptr, uint64_t b_ptr) {
    return leaf_equals<TFA_Leaf_Contents>(a_ptr, b_ptr);
}

uint64_t tfa_leaf_hash(const uint64_t contents_ptr, const uint64_t seed) {
    auto leaf_contents = reinterpret_cast<TFA_Leaf_Contents*>(contents_ptr);
    return (leaf_contents->post + (1ull << 63) * leaf_contents->is_accepting);
}

char* tfa_leaf_to_str(int comp, uint64_t leaf_val, char *buf, size_t buflen) {
   (void)comp;

    std::stringstream ss;
    auto leaf_contents = reinterpret_cast<TFA_Leaf_Contents*>(leaf_val);

    ss << "{.accepting=" << (leaf_contents->is_accepting ? "True" : "False")
       << ", .post=" << leaf_contents->post << "}";

    const std::string str(ss.str());

    const size_t required_buf_size = str.size() + 1; // With '\0' at the end
    if (required_buf_size <= buflen) {
        const char *cstr = str.c_str();
        std::memcpy(buf, cstr, sizeof(char) * required_buf_size);
        return buf;
    } else {
        char *new_buf = (char *)malloc(sizeof(char) * required_buf_size);
        std::memcpy(new_buf, str.c_str(), sizeof(char) * required_buf_size);
        return new_buf;
    }
}

MTBDD make_tfa_intersection_leaf(TFA_Leaf_Intersection_Contents* contents) {
    sylvan::MTBDD leaf = sylvan::mtbdd_makeleaf(mtbdd_tfa_intersection_leaf_type_id, (u64) contents);
    return leaf;
}

void mk_tfa_intersection_leaf_from_contents(u64* dest_contents) {
    mk_leaf_from_contents<TFA_Leaf_Intersection_Contents>(dest_contents);
}

void destroy_tfa_intersection_leaf(u64 leaf_value) {
    destroy_leaf<TFA_Leaf_Intersection_Contents>(leaf_value);
}

int tfa_intersection_leaf_equals(u64 a_ptr, u64 b_ptr) {
    return leaf_equals<TFA_Leaf_Intersection_Contents>(a_ptr, b_ptr);
}

u64 tfa_intersection_leaf_hash(const u64 contents_ptr, const u64 seed) {
    auto contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(contents_ptr);

    u64 hash = hash_vector(contents->post, seed);
    hash *= (contents->is_accepting * 33);

    contents->hash = hash;
    contents->hash_valid = true;

    return hash;
}

char* tfa_leaf_intersection_to_str(int comp, u64 leaf_val, char *buf, size_t buflen) {
    (void)comp;

    std::stringstream ss;
    auto leaf_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(leaf_val);

    ss << "{.accepting=" << (leaf_contents->is_accepting ? "True" : "False")
       << ", .post=" << leaf_contents->post << "}";

    const std::string str(ss.str());

    const size_t required_buf_size = str.size() + 1; // With '\0' at the end
    if (required_buf_size <= buflen) {
        const char *cstr = str.c_str();
        std::memcpy(buf, cstr, sizeof(char) * required_buf_size);
        return buf;
    } else {
        char *new_buf = (char *)malloc(sizeof(char) * required_buf_size);
        std::memcpy(new_buf, str.c_str(), sizeof(char) * required_buf_size);
        return new_buf;
    }
}

TASK_IMPL_3(MTBDD, tfa_mtbdd_union, MTBDD *, left_op_ptr, MTBDD *, right_op_ptr, u64, param)
{
    MTBDD left_mtbdd = *left_op_ptr, right_mtbdd = *right_op_ptr;

    if (left_mtbdd == mtbdd_false)  return right_mtbdd;
    if (right_mtbdd == mtbdd_false) return left_mtbdd;

    if (mtbdd_isleaf(left_mtbdd) && mtbdd_isleaf(right_mtbdd)) {
        assert(false);
    }

    return sylvan::mtbdd_invalid;
}

MTBDD perform_tfa_mtbdd_union(MTBDD left, MTBDD right) {
    LACE_ME;
    return mtbdd_applyp(left, right, 0u, TASK(tfa_mtbdd_union), AMAYA_TFA_UNION_OP_ID);
}

TASK_IMPL_3(MTBDD, tfa_mtbdd_intersection, MTBDD *, left_op_ptr, MTBDD *, right_op_ptr, u64, param)
{
    MTBDD left_mtbdd = *left_op_ptr, right_mtbdd = *right_op_ptr;

    // If left or right are false, then there is no post since we are doing an intersection intersection
    if (left_mtbdd == mtbdd_false)  return sylvan::mtbdd_false;

    if (mtbdd_isleaf(left_mtbdd) && mtbdd_isleaf(right_mtbdd)) {
        auto left_raw_val = sylvan::mtbdd_getvalue(left_mtbdd);
        auto left_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(left_raw_val);

        auto right_raw_val = sylvan::mtbdd_getvalue(left_mtbdd);
        auto right_contents = reinterpret_cast<TFA_Leaf_Contents*>(left_raw_val);

        TFA_Leaf_Intersection_Contents result = {
           .post = left_contents->post,
           .is_accepting = left_contents->is_accepting || right_contents->is_accepting,
        };
        result.post.push_back(right_contents->post);

        MTBDD result_leaf = make_tfa_intersection_leaf(&result);
        return result_leaf;
    }

    return sylvan::mtbdd_invalid;
}

MTBDD perform_tfa_mtbdd_intersection(MTBDD left, MTBDD right) {
    LACE_ME;
    return mtbdd_applyp(left, right, 0u, TASK(tfa_mtbdd_intersection), AMAYA_TFA_INTERSECTION_OP_ID);
}

MTBDD make_tfa_intersection_top() {
    TFA_Leaf_Intersection_Contents c = {.post = {}, .is_accepting = false};
    return make_tfa_intersection_leaf(&c);
}

MTBDD make_tfa_pareto_leaf(TFA_Pareto_Leaf* contents) {
    sylvan::MTBDD leaf = sylvan::mtbdd_makeleaf(mtbdd_tfa_pareto_leaf_type_id, reinterpret_cast<u64>(contents));
    return leaf;
}

void mk_tfa_pareto_leaf_from_contents(u64* contents) {
    mk_leaf_from_contents<TFA_Pareto_Leaf>(contents);
}

void destroy_tfa_pareto_leaf(u64 leaf_value) {
    destroy_leaf<TFA_Pareto_Leaf>(leaf_value);
}

int tfa_pareto_leaf_equals(u64 a_ptr, u64 b_ptr) {
    return leaf_equals<TFA_Pareto_Leaf>(a_ptr, b_ptr);
}

u64 tfa_pareto_leaf_hash(const u64 contents_ptr, const u64 seed) {
    u64 hash = seed;
    auto* leaf_contents = reinterpret_cast<TFA_Pareto_Leaf*>(contents_ptr);
    hash = hash_combine(hash, leaf_contents->elements.hash(0));
    hash *= leaf_contents->is_accepting ? 33 : 1;
    return hash;
}

MTBDD make_union_singleton_pareto_leaf(MTBDD not_bot, u64 prefix_size) {
    auto right_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(mtbdd_getvalue(not_bot));
    TFA_Pareto_Leaf new_leaf_contents = {.elements = Pareto_Set(prefix_size)};
    new_leaf_contents.elements.insert(right_contents->post, 0);
    return make_tfa_pareto_leaf(&new_leaf_contents);
}

MTBDD make_union_with_singleton(MTBDD pareto_leaf, MTBDD singleton) {
    auto pareto_contents = reinterpret_cast<TFA_Pareto_Leaf*>(sylvan::mtbdd_getvalue(pareto_leaf));
    auto singleton_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(sylvan::mtbdd_getvalue(singleton));

    TFA_Pareto_Leaf leaf = {.elements = pareto_contents->elements};
    leaf.elements.insert(singleton_contents->post, 0);
    leaf.is_accepting = pareto_contents->is_accepting || leaf.is_accepting;

    return make_tfa_pareto_leaf(&leaf);
}

TASK_IMPL_3(MTBDD, tfa_mtbdd_pareto_union, MTBDD*, left_ptr, MTBDD*, right_ptr, u64, prefix_size) {
    MTBDD left = *left_ptr;
    MTBDD right = *right_ptr;

    if (!sylvan::mtbdd_isleaf(left) || !sylvan::mtbdd_isleaf(right)) return mtbdd_invalid;

    if (left == sylvan::mtbdd_false && right == sylvan::mtbdd_false) return mtbdd_false;

    if (left != mtbdd_false && right != mtbdd_false) {
        u64 left_type  = sylvan::mtbdd_gettype(left);
        u64 right_type = sylvan::mtbdd_gettype(right);
        if (left_type == mtbdd_tfa_pareto_leaf_type_id && right_type == mtbdd_tfa_pareto_leaf_type_id) {
            auto left_contents  = reinterpret_cast<TFA_Pareto_Leaf*>(sylvan::mtbdd_getvalue(left));
            auto right_contents = reinterpret_cast<TFA_Pareto_Leaf*>(sylvan::mtbdd_getvalue(right));
            TFA_Pareto_Leaf new_leaf = {
                .elements = merge_pareto_sets(left_contents->elements, right_contents->elements),
                .is_accepting = left_contents->is_accepting || right_contents->is_accepting
            };
            MTBDD result = make_tfa_pareto_leaf(&new_leaf);
            return result;
        } else if (left_type == mtbdd_tfa_intersection_leaf_type_id && right_type == mtbdd_tfa_intersection_leaf_type_id) {
            auto left_contents  = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(sylvan::mtbdd_getvalue(left));
            auto right_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(sylvan::mtbdd_getvalue(right));
            TFA_Pareto_Leaf new_leaf = {
                .elements = Pareto_Set(prefix_size),
                .is_accepting = left_contents->is_accepting || right_contents->is_accepting
            };
            new_leaf.elements.insert(left_contents->post, 0);
            new_leaf.elements.insert(right_contents->post, 0);
            MTBDD result = make_tfa_pareto_leaf(&new_leaf);
            return result;
        }

        // One is pareto, other one is intersection result
        MTBDD intersection_leaf = left;
        MTBDD pareto_leaf = right;
        if (left_type == mtbdd_tfa_pareto_leaf_type_id) std::swap(intersection_leaf, pareto_leaf);

        auto intersection_contents = reinterpret_cast<TFA_Leaf_Intersection_Contents*>(sylvan::mtbdd_getvalue(intersection_leaf));
        auto pareto_contents       = reinterpret_cast<TFA_Pareto_Leaf*>(sylvan::mtbdd_getvalue(pareto_leaf));

        TFA_Pareto_Leaf new_leaf = {
            .elements = pareto_contents->elements,
            .is_accepting = pareto_contents->is_accepting || intersection_contents->is_accepting
        };
        new_leaf.elements.insert(intersection_contents->post, 0);
        MTBDD result = make_tfa_pareto_leaf(&new_leaf);
        return result;
    }

    // One of them is false
    MTBDD not_bot = left;
    MTBDD bot     = right;
    if (not_bot == mtbdd_false) std::swap(not_bot, bot);

    u64 not_bot_type = mtbdd_gettype(not_bot);
    if (not_bot_type == mtbdd_tfa_pareto_leaf_type_id) return not_bot;

    assert(mtbdd_isleaf(not_bot));
    assert(mtbdd_gettype(not_bot) == mtbdd_tfa_intersection_leaf_type_id);

    auto result = make_union_singleton_pareto_leaf(not_bot, prefix_size);
    return result;
    
    return mtbdd_invalid;
}

MTBDD perform_tfa_pareto_union(MTBDD left, MTBDD right, u64 prefix_size) {
    LACE_ME;
    return mtbdd_applyp(left, right, prefix_size, TASK(tfa_mtbdd_pareto_union), AMAYA_TFA_PARETO_UNION_OP_ID);
}

u64 projection_prefix_size = 0;

TASK_IMPL_3(MTBDD, tfa_mtbdd_pareto_projection_op, MTBDD, left_mtbdd, MTBDD, right_mtbdd, int, k) {
    MTBDD u = perform_tfa_pareto_union(left_mtbdd, right_mtbdd, projection_prefix_size);
    return u;
}

MTBDD perform_tfa_pareto_projection(MTBDD bdd, BDDSET var_to_project_away, u64 prefix_size) {
    LACE_ME;
    projection_prefix_size = prefix_size;
    MTBDD result = mtbdd_abstract(bdd, var_to_project_away, TASK(tfa_mtbdd_pareto_projection_op));
    return result;
}

char* tfa_pareto_leaf_to_str(int comp, u64 leaf_val, char *buf, size_t buflen) {
   (void) comp;

    std::stringstream ss;
    auto leaf_contents = reinterpret_cast<TFA_Pareto_Leaf*>(leaf_val);

    ss << "{.accepting=" << (leaf_contents->is_accepting ? "True" : "False")
       << ", .post=" << leaf_contents->elements<< "}";

    const std::string str(ss.str());

    const size_t required_buf_size = str.size() + 1; // With '\0' at the end
    if (required_buf_size <= buflen) {
        const char *cstr = str.c_str();
        std::memcpy(buf, cstr, sizeof(char) * required_buf_size);
        return buf;
    } else {
        char *new_buf = (char*) malloc(sizeof(char) * required_buf_size);
        std::memcpy(new_buf, str.c_str(), sizeof(char) * required_buf_size);
        return new_buf;
    }
}

TASK_IMPL_3(MTBDD, are_mtbdds_identic_op, sylvan::MTBDD*, left_ptr, sylvan::MTBDD*, right_ptr, u64, result_raw_ptr) {
    MTBDD left  = *left_ptr;
    MTBDD right = *right_ptr;

    if (!mtbdd_isleaf(left) || !mtbdd_isleaf(right)) return mtbdd_invalid;

    if (left == mtbdd_false) return right == mtbdd_false ? mtbdd_true : mtbdd_false;
    if (right == mtbdd_false) return left == mtbdd_false ? mtbdd_true : mtbdd_false;

    u64 left_type = mtbdd_gettype(left);
    u64 right_type = mtbdd_gettype(right);

    if (left_type != right_type) return mtbdd_false; 

    u64 left_contents_ptr  = mtbdd_getvalue(left);
    u64 right_contents_ptr = mtbdd_getvalue(right);

    if (left_type == mtbdd_tfa_pareto_leaf_type_id) {
        return tfa_pareto_leaf_equals(left_contents_ptr, right_contents_ptr) ? mtbdd_true : mtbdd_false;
    }

    assert(false);
    return mtbdd_invalid;
}

bool check_mtbdds_are_identical(MTBDD left, MTBDD right) {
    LACE_ME;
    return mtbdd_applyp(sylvan::mtbdd_false, sylvan::mtbdd_false, 0, TASK(are_mtbdds_identic_op), AMAYA_ARE_MTBDDS_IDENTIC_OP);
}
