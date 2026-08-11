#include <quaternion.h>
#include "internal.h"
#include <stdlib.h>
// get shortest equivalent ideal

#define DEBUG_PRINTS 1

void
quat_lideal_shortest_equivalent(quat_left_ideal_t *equiv,
                                quat_alg_elem_t *elem,
                                const quat_left_ideal_t *lideal,
                                const quat_alg_t *alg)
{

    ibz_mat_4x4_t gram, red;
    quat_alg_elem_t new_alpha;
    quat_alg_elem_init(&new_alpha);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);

    // computing the reduced basis
    quat_lideal_reduce_basis(&red, &gram, lideal, alg);

    ibz_set(&(new_alpha.coord[0]), 1);
    ibz_mat_4x4_eval(&new_alpha.coord, &red, &new_alpha.coord);
    ibz_copy(&new_alpha.denom, &lideal->lattice.denom);
    assert(quat_lattice_contains(NULL, &lideal->lattice, &new_alpha));
#ifndef NDEBUG
    ibz_t n, d;
    ibz_init(&n);
    ibz_init(&d);
    quat_alg_norm(&n, &d, &new_alpha, alg);
    assert(ibz_is_one(&d));
    ibz_div(&n, &d, &n, &lideal->norm);
    assert(ibz_is_zero(&d));
    ibz_finalize(&n);
    ibz_finalize(&d);
#endif
    quat_alg_elem_copy(elem, &new_alpha);
    ibz_mul(&new_alpha.denom, &new_alpha.denom, &lideal->norm);
    ibz_neg(&new_alpha.coord[0], &new_alpha.coord[0]);
    quat_lideal_mul(equiv, lideal, &new_alpha, alg);

    quat_alg_elem_finalize(&new_alpha);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);
}

void
quat_lideal_generator_small_coprime(quat_alg_elem_t *gen,
                                    const quat_left_ideal_t *lideal,
                                    const quat_alg_t *alg,
                                    int sampling_bound_bits)
{
    int found = 0;
    ibz_t n, d, n2, gcd;
    ibz_vec_4_t coeffs;
    ibz_init(&n);
    ibz_init(&d);
    ibz_init(&n2);
    ibz_init(&gcd);
    ibz_vec_4_init(&coeffs);
    ibz_copy(&gen->denom, &lideal->lattice.denom);
    ibz_mul(&n2, &lideal->norm, &lideal->norm);
    while (!found) {
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval_bits(&(coeffs[i]), sampling_bound_bits);
        }
        ibz_mat_4x4_eval(&(gen->coord), &(lideal->lattice.basis), &coeffs);

        // check a_alpha invertible
        if (ibz_is_one(&(gen->denom))) {
            ibz_mul(&gcd, &(gen->coord[0]), &ibz_const_two);
            ibz_gcd(&gcd, &lideal->norm, &gcd);
            found = ibz_is_one(&gcd);
        } else {
            found = 1;
            if (0 == ibz_cmp(&ibz_const_two, &(gen->denom))) {
                ibz_gcd(&gcd, &lideal->norm, &(gen->coord[0]));
                found = ibz_is_one(&gcd);
            }
        }
        // check generator
        if (found) {
            quat_alg_norm(&n, &d, gen, alg);
            assert(ibz_is_one(&d));
            ibz_gcd(&gcd, &n2, &n);
            found = (ibz_cmp(&gcd, &lideal->norm) == 0);
        }
    }

    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&n);
    ibz_finalize(&d);
    ibz_finalize(&n2);
    ibz_finalize(&gcd);
}

// enumerate in dim 2

// helper for cvp
void
ibz_rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t r, sign_q, abs_b;
    ibz_init(&r);
    ibz_init(&sign_q);
    ibz_init(&abs_b);

    // assumed to round towards 0
    ibz_abs(&abs_b, b);
    // q is of same sign as a*b (and 0 if a is 0)
    ibz_mul(&sign_q, a, b);
    ibz_div(q, &r, a, b);
    ibz_abs(&r, &r);
    ibz_add(&r, &r, &r);
    ibz_set(&sign_q, (1 - 2 * (ibz_cmp(&sign_q, &ibz_const_zero) < 0)) * (ibz_cmp(&r, &abs_b) > 0));
    ibz_add(q, q, &sign_q);
    ibz_finalize(&r);
    ibz_finalize(&sign_q);
    ibz_finalize(&abs_b);
}

int
quat_dim2_lattice_contains(const ibz_mat_2x2_t *basis, const ibz_t *coord1, const ibz_t *coord2)
{
    int res = 1;
    ibz_t prod, sum, det, r;
    ibz_init(&det);
    ibz_init(&r);
    ibz_init(&sum);
    ibz_init(&prod);
    // compute det, then both coordinates (inverse*det)*vec, where vec is (coord1, coord2) and check wthether det
    // divides both results
    ibz_mat_2x2_det_from_ibz(&det, &((*basis)[0][0]), &((*basis)[0][1]), &((*basis)[1][0]), &((*basis)[1][1]));
    ibz_mul(&sum, coord1, &((*basis)[1][1]));
    ibz_mul(&prod, coord2, &((*basis)[0][1]));
    ibz_sub(&sum, &sum, &prod);
    ibz_div(&prod, &r, &sum, &det);
    res = res & ibz_is_zero(&r);
    ibz_mul(&sum, coord2, &((*basis)[0][0]));
    ibz_mul(&prod, coord1, &((*basis)[1][0]));
    ibz_sub(&sum, &sum, &prod);
    ibz_div(&prod, &r, &sum, &det);
    res = res & ibz_is_zero(&r);
    ibz_finalize(&det);
    ibz_finalize(&r);
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    return (res);
}

void
quat_dim2_lattice_norm(ibz_t *norm, const ibz_t *coord1, const ibz_t *coord2, const ibz_t *norm_q)
{
    ibz_t prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_mul(&sum, coord1, coord1);
    ibz_mul(&prod, coord2, coord2);
    ibz_mul(&prod, &prod, norm_q);
    ibz_add(norm, &sum, &prod);
    ibz_finalize(&prod);
    ibz_finalize(&sum);
}

void
quat_dim2_lattice_bilinear(ibz_t *res,
                           const ibz_t *v11,
                           const ibz_t *v12,
                           const ibz_t *v21,
                           const ibz_t *v22,
                           const ibz_t *norm_q)
{
    ibz_t prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_mul(&sum, v11, v21);
    ibz_mul(&prod, v12, v22);
    ibz_mul(&prod, &prod, norm_q);
    ibz_add(res, &sum, &prod);
    ibz_finalize(&prod);
    ibz_finalize(&sum);
}

// algo 3.1.14 Cohen (exact solution for shortest vector in dimension 2, than take a second, orthogonal vector)
void
quat_dim2_lattice_short_basis(ibz_mat_2x2_t *reduced, const ibz_mat_2x2_t *basis, const ibz_t *norm_q)
{
    ibz_vec_2_t a, b, t;
    ibz_t prod, sum, norm_a, norm_b, r, norm_t, n;
    ibz_vec_2_init(&a);
    ibz_vec_2_init(&b);
    ibz_vec_2_init(&t);
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_init(&r);
    ibz_init(&n);
    ibz_init(&norm_t);
    ibz_init(&norm_a);
    ibz_init(&norm_b);
    // init a,b
    ibz_copy(&(a[0]), &((*basis)[0][0]));
    ibz_copy(&(a[1]), &((*basis)[1][0]));
    ibz_copy(&(b[0]), &((*basis)[0][1]));
    ibz_copy(&(b[1]), &((*basis)[1][1]));
    // compute initial norms
    quat_dim2_lattice_norm(&norm_a, &(a[0]), &(a[1]), norm_q);
    quat_dim2_lattice_norm(&norm_b, &(b[0]), &(b[1]), norm_q);
    // exchange if needed
    if (ibz_cmp(&norm_a, &norm_b) < 0) {
        ibz_copy(&sum, &(a[0]));
        ibz_copy(&(a[0]), &(b[0]));
        ibz_copy(&(b[0]), &sum);
        ibz_copy(&sum, &(a[1]));
        ibz_copy(&(a[1]), &(b[1]));
        ibz_copy(&(b[1]), &sum);
        ibz_copy(&sum, &norm_a);
        ibz_copy(&norm_a, &norm_b);
        ibz_copy(&norm_b, &sum);
    }
    int test = 1;
    while (test) {
        // compute n
        quat_dim2_lattice_bilinear(&n, &(a[0]), &(a[1]), &(b[0]), &(b[1]), norm_q);
        // set r
        ibz_rounded_div(&r, &n, &norm_b);
        // compute t_norm
        ibz_set(&prod, 2);
        ibz_mul(&prod, &prod, &n);
        ibz_mul(&prod, &prod, &r);
        ibz_sub(&sum, &norm_a, &prod);
        ibz_mul(&prod, &r, &r);
        ibz_mul(&prod, &prod, &norm_b);
        ibz_add(&norm_t, &sum, &prod);
        // test:
        if (ibz_cmp(&norm_b, &norm_t) > 0) {
            // compute t, a, b
            ibz_copy(&norm_a, &norm_b);
            ibz_copy(&norm_b, &norm_t);
            // t is a -rb, a is b, b is t
            ibz_mul(&prod, &r, &(b[0]));
            ibz_sub(&(t[0]), &(a[0]), &prod);
            ibz_mul(&prod, &r, &(b[1]));
            ibz_sub(&(t[1]), &(a[1]), &prod);
            ibz_copy(&(a[0]), &(b[0]));
            ibz_copy(&(a[1]), &(b[1]));
            ibz_copy(&(b[0]), &(t[0]));
            ibz_copy(&(b[1]), &(t[1]));
        } else {
            test = 0;
        }
    }
    // output : now b is short: need to get 2nd short vector: idea: take shortest among t and a
    if (ibz_cmp(&norm_t, &norm_a) < 0) {
        ibz_mul(&prod, &r, &(b[0]));
        ibz_sub(&(a[0]), &(a[0]), &prod);
        ibz_mul(&prod, &r, &(b[1]));
        ibz_sub(&(a[1]), &(a[1]), &prod);
    }
    ibz_copy(&((*reduced)[0][0]), &(b[0]));
    ibz_copy(&((*reduced)[1][0]), &(b[1]));
    ibz_copy(&((*reduced)[0][1]), &(a[0]));
    ibz_copy(&((*reduced)[1][1]), &(a[1]));

    ibz_finalize(&prod);
    ibz_finalize(&sum);
    ibz_finalize(&norm_a);
    ibz_finalize(&norm_b);
    ibz_finalize(&norm_t);
    ibz_vec_2_finalize(&a);
    ibz_vec_2_finalize(&b);
    ibz_vec_2_finalize(&t);
    ibz_finalize(&r);
    ibz_finalize(&n);
}

// qlapoti
int
quat_qlapoti_check_mod_condition(const ibz_t *m, const ibz_t *a, const ibz_t *b)
{
    int ok = 1;
    int A = ibz_get(a) & 1;
    int B = ibz_get(b) & 1;
    int M8 = ibz_get(m) & 7;
    int M = M8 & 3;
#ifndef NDEBUG
    if (DEBUG_PRINTS > 3)
        ibz_printf("m %d a %d b %d\n", M, A, B);
#endif
    if (M8 == 0)
        ok = 0;
    if (ok) {
        if ((A == B) && (A == 0)) {
            if (M != 0) {
                ok = 0;
#ifndef NDEBUG
                if (DEBUG_PRINTS > 2)
                    ibz_printf("case 0\n");
#endif
            }
        } else {
            if (((A == B) && (A == 1))) {
                if (M != 2) {
                    ok = 0;
#ifndef NDEBUG
                    if (DEBUG_PRINTS > 2)
                        ibz_printf("case 1\n");
#endif
                }
            } else {
                if (M != 1) {
#ifndef NDEBUG
                    if (DEBUG_PRINTS > 2)
                        ibz_printf("mod 4\n");
#endif
                    ok = 0;
                }
            }
        }
    }

    return (ok);
}

int
quat_dim2_lattice_qlapoti_cvp_condition(quat_alg_elem_t *elem, const ibz_vec_2_t *vec, const void *params)
{
    // remember a_alpha and b_alpha are integers
    int found = 1;
    qlapoti_enumeration_parameters_t *q_params = ((qlapoti_enumeration_parameters_t *)params);
    ibz_t m2, tmp, A, B;
    ibz_vec_2_t a, b;
    ibz_init(&m2);
    ibz_init(&tmp);
    ibz_init(&A);
    ibz_init(&B);
    ibz_vec_2_init(&a);
    ibz_vec_2_init(&b);

    // compute A,B
    // Assumes vec is target-close (up to one sign)
    ibz_neg(&A, &((*vec)[0]));
    ibz_neg(&B, &((*vec)[1]));
#ifndef NDEBUG
    if (DEBUG_PRINTS > 3)
        ibz_printf("A %Zd\nB %Zd\n", &A, &B);
    // assert (2*(a_alpha*A + b_alpha*B)) % N == M % N
    ibz_mul(&tmp, &A, q_params->a_alpha);
    ibz_mul(&m2, &B, q_params->b_alpha);
    ibz_add(&tmp, &tmp, &m2);
    ibz_mod(&tmp, &tmp, q_params->n);
    ibz_mod(&m2, q_params->m, q_params->n);
    assert(0 == ibz_cmp(&m2, &tmp));
#endif
    // compute m2
    // M2 = M - 2 *a_alpha *A - 2 *b_alpha *B; M2 = ZZ(M2 / N)
    ibz_mul(&tmp, &A, q_params->a_alpha);
    ibz_mul(&m2, &B, q_params->b_alpha);
    ibz_add(&tmp, &tmp, &m2);
    ibz_sub(&m2, q_params->m, &tmp);
    ibz_div(&m2, &tmp, &m2, q_params->n);
    assert(ibz_is_zero(&tmp));
    // Complete the square
    // Comp: M3 = M2 - A ** 2 - B ** 2
    ibz_mul(&tmp, &A, &A);
    ibz_sub(&m2, &m2, &tmp);
    ibz_mul(&tmp, &B, &B);
    ibz_sub(&m2, &m2, &tmp);

    // Comp; M4 = 2 * M3 + A * *2 + B * *2
    ibz_add(&m2, &m2, &m2);
    ibz_add(&m2, &m2, &tmp);
    ibz_mul(&tmp, &A, &A);
    ibz_add(&m2, &m2, &tmp);

#ifndef NDEBUG
    if ((DEBUG_PRINTS > 2))
        ibz_printf("M4 %Zd\n", &m2);
    // If the first one is too small, there is no point in trying others...
    // if M4 < 0:if first_vec : break else : continue
    // Must be communicated to exterior loop (enforce change of alpha)? Might be fine already since better enum (if
    // correct bound)
    // Test for unsolvable cases
    if ((DEBUG_PRINTS > 2) && found && !(ibz_cmp(&ibz_const_zero, &m2) < 0))
        ibz_printf("size\n");
#endif
    if (found)
        found = quat_qlapoti_check_mod_condition(&m2, &A, &B);
    if (found) {
        // cornacchia
        found = ibz_cornacchia_extended(&(a[0]), &(b[0]), &m2, q_params->cornacchia_params);
#ifndef NDEBUG
        if ((DEBUG_PRINTS > 2) && !found)
            ibz_printf("cor \n");
#endif
    }
    if (found) {
// treat output
#ifndef NDEBUG
        if (DEBUG_PRINTS > 3)
            ibz_printf("a_0 %Zd, b_0 %Zd\n", &(a[0]), &(b[0]));
#endif
        if (!((ibz_get(&A) & 1) == (ibz_get(&(a[0])) & 1)))
            ibz_swap(&(a[0]), &(b[0]));
        assert(((ibz_get(&A) & 1) == (ibz_get(&(a[0])) & 1)));
        assert(((ibz_get(&B) & 1) == (ibz_get(&(b[0])) & 1)));
        // a1 = ZZ((ad1 + A)/2);b1 = ZZ((bd1 + B)/2)
        ibz_add(&(a[0]), &(a[0]), &A);
        ibz_div(&(a[0]), &tmp, &(a[0]), &ibz_const_two);
        assert(ibz_is_zero(&tmp));
        ibz_add(&(b[0]), &(b[0]), &B);
        ibz_div(&(b[0]), &tmp, &(b[0]), &ibz_const_two);
        assert(ibz_is_zero(&tmp));
        // a2 = A - a1;b2 = B - b1
        ibz_sub(&(a[1]), &A, &(a[0]));
        ibz_sub(&(b[1]), &B, &(b[0]));
        ibz_copy(&(elem->coord[0]), &(a[0]));
        ibz_copy(&(elem->coord[1]), &(a[1]));
        ibz_copy(&(elem->coord[2]), &(b[0]));
        ibz_copy(&(elem->coord[3]), &(b[1]));
    }
    ibz_finalize(&m2);
    ibz_finalize(&tmp);
    ibz_finalize(&A);
    ibz_finalize(&B);
    ibz_vec_2_finalize(&a);
    ibz_vec_2_finalize(&b);
    return (found);
}

int
quat_elem_is_odd_norm(const quat_alg_elem_t *elem)
{
    int found = 0;
    if ((ibz_get(&elem->coord[0]) & 1) == (ibz_get(&elem->coord[2]) & 1)) {
        return found;
    }
    for (int i = 0; i < 4; i++) {
        if (found) {
            if ((ibz_get(&elem->coord[i]) & 3) == 2) {
                return 0;
            }
        } else {
            if ((ibz_get(&elem->coord[i]) & 3) == 2) {
                found = 1;
            }
        }
    }
    return found;
}

int
get_endtype(const quat_alg_elem_t *elem)
{
    int t1, t2, t3, t4;
    // transformlist = {1: [[2, 0, 1, 2], [2, 2, 3, 0], [2, 2, 1, 0], [2, 0, 3, 2]], 2: [[0, 2, 2, 1], [0, 2, 2, 3], [2,
    // 2, 0, 1], [2, 2, 0, 3]]}
    t1 = ibz_get(&elem->coord[0]) & 3;
    t2 = ibz_get(&elem->coord[1]) & 3;
    // NB! These are swapped because I computed this using j and k swapped
    t3 = ibz_get(&elem->coord[3]) & 3;
    t4 = ibz_get(&elem->coord[2]) & 3;
    // End NB!
    if (t1 == 2) {
        if (t2 == 2) {
            if (t3 == 1 && t4 == 0) {
                return 1;
            } else if (t3 == 3 && t4 == 2) {
                return 1;
            } else if (t3 == 0) {
                if (t4 == 1) {
                    return 2;
                } else if (t4 == 3) {
                    return 2;
                } else {
                    return 0;
                }
            } else {
                return 0;
            }
        } else if (t2 == 0) {
            if (t3 == 1 && t4 == 2) {
                return 1;
            } else if (t3 == 3 && t4 == 2) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    } else if (t1 == 0) {
        if (t2 == 2 && t3 == 2) {
            if (t4 == 1) {
                return 2;
            } else if (t4 == 3) {
                return 2;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}
// document all functions
int
quat_qlapoti(quat_alg_elem_t *mu1,
             quat_alg_elem_t *mu2,
             quat_alg_elem_t *theta,
             quat_alg_elem_t *smallest,
             const quat_left_ideal_t *lideal,
             const quat_alg_t *alg,
             int max_counter_alpha,
             int gen_sampling_bound_bits,
             int two_power,
             const ibz_cornacchia_extended_params_t *cornacchia_params)
{
    int found = 0;
    int transformed = 0;
    int keep_alpha = 0;
    quat_left_ideal_t small;
    ibz_mat_4x4_t gram;
    quat_alg_elem_t alpha, encoding, alpha_0;
    ibz_mat_2x2_t L, L_red, L_inv;
    ibz_vec_2_t v_target, v_close, v_diff;
    ibz_t n, m, norm_n, norm_d, a_alpha, b_alpha, x, T, alpha0norm, psqrt;
    ibz_t tmp, lam, two_e, L_det;
    qlapoti_enumeration_parameters_t params;
    quat_alg_elem_t gamma1, gamma2, temp;
    quat_alg_elem_init(&gamma1);
    quat_alg_elem_init(&gamma2);
    ibz_init(&T);
    ibz_init(&L_det);
    ibz_vec_2_init(&v_target);
    ibz_vec_2_init(&v_close);
    ibz_vec_2_init(&v_diff);
    ibz_mat_2x2_init(&L);
    ibz_mat_2x2_init(&L_red);
    ibz_mat_2x2_init(&L_inv);
    ibz_init(&tmp);
    ibz_init(&x);
    ibz_init(&n);
    ibz_init(&m);
    ibz_init(&two_e);
    ibz_init(&lam);
    ibz_init(&norm_n);
    ibz_init(&norm_d);
    ibz_init(&a_alpha);
    ibz_init(&b_alpha);
    ibz_init(&alpha0norm);
    ibz_init(&psqrt);
    ibz_mat_4x4_init(&gram);
    quat_alg_elem_init(&alpha);
    quat_alg_elem_init(&alpha_0);
    quat_alg_elem_init(&encoding);
    quat_alg_elem_init(&temp);
    quat_left_ideal_init(&small);
    ibz_sqrt_floor(&psqrt, &alg->p);
    quat_lideal_shortest_equivalent(&small, smallest, lideal, alg);
    params.cornacchia_params = cornacchia_params;
    quat_lideal_reduce_basis(&small.lattice.basis, &gram, &small, alg);
    ibz_copy(&n, &small.norm);
    ibz_pow(&two_e, &ibz_const_two, two_power);

    // int num_alphas_tried = 0;
    // int num_loops_total = 0;
    for (int counter = 0; counter < max_counter_alpha * 20; counter++) {
        found = 1;

        // compute alpha, lam, alpha_0, alpha0norm
        {
            if (counter % max_counter_alpha == 0) {
                keep_alpha = 0;
                ibz_set(&lam, 1);
            }
            if (!keep_alpha) {
                // num_alphas_tried++;
                quat_lideal_generator_small_coprime(&alpha_0, &small, alg, gen_sampling_bound_bits);
                quat_alg_elem_copy(&alpha, &alpha_0);
                quat_alg_norm(&alpha0norm, &b_alpha, &alpha_0, alg);
                assert(ibz_is_one(&b_alpha));
                ibz_set(&lam, 1);
            } else {
                ibz_add(&lam, &lam, &ibz_const_one);
                quat_alg_add(&alpha, &alpha, &alpha_0);
                ibz_gcd(&b_alpha, &(small.norm), &lam);
                while (!ibz_is_one(&b_alpha)) {
                    ibz_add(&lam, &lam, &ibz_const_one);
                    quat_alg_add(&alpha, &alpha, &alpha_0);
                    ibz_gcd(&b_alpha, &(small.norm), &lam);
                }
            }
            assert(quat_lattice_contains(NULL, &small.lattice, &alpha));
            quat_alg_normalize(&alpha);
        }

        // compute m
        // 2^e-(2norm(alpha0)lam^2)/n
        {
            ibz_add(&norm_n, &alpha0norm, &alpha0norm);
            ibz_mul(&norm_d, &lam, &lam);
            ibz_mul(&norm_n, &norm_n, &norm_d);
            ibz_div_floor(&norm_n, &norm_d, &norm_n, &n);
            assert(ibz_is_zero(&norm_d));
            ibz_sub(&m, &two_e, &norm_n);
            // Discard bad alphas directly: this case should never happen
            if (ibz_cmp(&m, &ibz_const_zero) < 0) {
                found = 0;
                continue;
            }
            assert(ibz_cmp(&m, &ibz_const_zero) > 0);
        }
        // set a_alpha, b_alpha (to double of what they are in sage)
        {
            ibz_copy(&a_alpha, &(alpha.coord[0]));
            ibz_copy(&b_alpha, &(alpha.coord[1]));
            if (ibz_is_one(&alpha.denom)) {
                ibz_add(&a_alpha, &a_alpha, &a_alpha);
                ibz_add(&b_alpha, &b_alpha, &b_alpha);
            }
#ifndef NDEBUG
            else {
                assert(ibz_cmp(&alpha.denom, &ibz_const_two) == 0);
            }
#endif
        }
        // Prepare target vector: T, v_target
        {
            // #A + (2*b_alpha/2*a_alpha)B = M/(2*a_alpha) (mod N)
            // x = ZZ(Z_N(2*b_alpha)*(Z_N(2*a_alpha)**-1))
            ibz_invmod(&tmp, &a_alpha, &n);
            // T = ZZ(Z_N(M)*Z_N(2*a_alpha)**-1)
            ibz_mul(&T, &tmp, &m);
            ibz_mod(&T, &T, &n);
            // v = vector(ZZ, [-T, 0])
            ibz_set(&(v_target[1]), 0);
            ibz_neg(&(v_target[0]), &T);
        }
        // Prepare 2d lattice: L, L_red, L_inv. Set keep_alpha accordingly
        // tmp must contain invmod(a_alpha,n) here
        if (!keep_alpha) {
            // build lattice
            ibz_mul(&x, &b_alpha, &tmp);
            ibz_mod(&x, &x, &n);
            // L = Matrix(ZZ, [[N-x, 1], [N, 0]])
            ibz_sub(&(L[0][0]), &n, &x);
            ibz_copy(&(L[0][1]), &n);
            ibz_set(&(L[1][1]), 0);
            ibz_set(&(L[1][0]), 1);
            quat_dim2_lattice_short_basis(&L_red, &L, &ibz_const_one);
            ibz_mat_2x2_inv_with_det_as_denom(&L_inv, &L_det, &L_red);
            // set keep_alpha: test short basis is suitable
            {
                ibz_gcd(&tmp, &a_alpha, &b_alpha);
                keep_alpha = ibz_is_one(&tmp);
                ibz_mul(&norm_d, &L_red[0][1], &L_red[0][1]);
                ibz_mul(&norm_n, &L_red[1][1], &L_red[1][1]);
                ibz_add(&norm_n, &norm_d, &norm_n);
                keep_alpha = keep_alpha && (ibz_cmp(&norm_n, &psqrt) < 0);
            }
        }
        // find close vector
        {
            ibz_mat_2x2_eval(&v_close, &L_inv, &v_target);
            ibz_rounded_div(&(v_close[0]), &(v_close[0]), &L_det);
            ibz_rounded_div(&(v_close[1]), &(v_close[1]), &L_det);
            ibz_mat_2x2_eval(&v_close, &L_red, &v_close);
            ibz_sub(&(v_diff[0]), &(v_target[0]), &(v_close[0]));
            ibz_sub(&(v_diff[1]), &(v_target[1]), &(v_close[1]));
        }
        // set parameters and call condition
        params.a_alpha = &a_alpha;
        params.b_alpha = &b_alpha;
        params.m = &m;
        params.n = &n;
        found = quat_dim2_lattice_qlapoti_cvp_condition(&encoding, &v_diff, &params);

        if (!found) {
            continue;
        } else {
            // extract ideals from encoding of output
            ibz_mul(&(gamma1.coord[0]), &(encoding.coord[0]), &n);
            ibz_mul(&(gamma2.coord[0]), &(encoding.coord[1]), &n);
            ibz_mul(&(gamma1.coord[1]), &(encoding.coord[2]), &n);
            ibz_mul(&(gamma2.coord[1]), &(encoding.coord[3]), &n);
            quat_alg_add(&gamma1, &gamma1, &alpha);
            quat_alg_add(&gamma2, &gamma2, &alpha);
#ifndef NDEBUG
            ibz_t n1, n2, d;
            ibz_init(&n1);
            ibz_init(&n2);
            ibz_init(&d);
            assert(quat_lattice_contains(NULL, &small.lattice, &gamma1));
            assert(quat_lattice_contains(NULL, &small.lattice, &gamma2));
            quat_alg_norm(&n1, &d, &gamma1, alg);
            assert(ibz_is_one(&d));
            quat_alg_norm(&n2, &d, &gamma2, alg);
            assert(ibz_is_one(&d));
            ibz_div(&n1, &d, &n1, &(small.norm));
            assert(ibz_is_zero(&d));
            ibz_div(&n2, &d, &n2, &(small.norm));
            assert(ibz_is_zero(&d));
            ibz_add(&n2, &n1, &n2);
            assert(0 == ibz_cmp(&two_e, &n2));
            ibz_finalize(&n1);
            ibz_finalize(&n2);
            ibz_finalize(&d);
#endif
            quat_alg_elem_copy(mu1, &gamma1);
            quat_alg_elem_copy(mu2, &gamma2);
            quat_alg_conj(&(gamma1), &(gamma1));
            ibz_mul(&gamma1.denom, &gamma1.denom, &n);
            quat_alg_mul(theta, mu2, &gamma1, alg);
#ifndef NDEBUG
            quat_left_ideal_t I1, I2;
            quat_left_ideal_init(&I1);
            quat_left_ideal_init(&I2);
            ibz_neg(&(gamma2.coord[0]), &(gamma2.coord[0]));
            ibz_mul(&gamma2.denom, &gamma2.denom, &n);
            quat_lideal_mul(&I1, &small, &gamma1, alg);
            quat_lideal_mul(&I2, &small, &gamma2, alg);
            ibz_add(&tmp, &I1.norm, &I2.norm);
            assert(ibz_cmp(&tmp, &two_e) == 0);
            quat_left_ideal_finalize(&I1);
            quat_left_ideal_finalize(&I2);
#endif

            // enforces suitable output
            quat_alg_normalize(theta);
            if (ibz_is_one(&theta->denom)) {
                int endtype = get_endtype(theta);
                if (endtype == 1) {
                    quat_alg_add(&temp, mu1, mu2);
                    ibz_add(&temp.denom, &temp.denom, &temp.denom);
                    quat_alg_sub(mu1, mu1, mu2);
                    ibz_add(&mu1->denom, &mu1->denom, &mu1->denom);
                    quat_alg_elem_copy(mu2, &temp);
                } else if (endtype == 2) {
                    quat_alg_elem_set(&temp, 1, 0, 1, 0, 0);
                    quat_alg_mul(mu2, &temp, mu2, alg);
                    quat_alg_add(&temp, mu1, mu2);
                    ibz_add(&temp.denom, &temp.denom, &temp.denom);
                    quat_alg_sub(mu1, mu1, mu2);
                    ibz_add(&mu1->denom, &mu1->denom, &mu1->denom);
                    quat_alg_elem_copy(mu2, &temp);
                } else {
                    keep_alpha = 0;
                    ibz_set(&lam, 1);
                    ibz_vec_4_set(&gamma1.coord, 0, 0, 0, 0);
                    ibz_vec_4_set(&gamma2.coord, 0, 0, 0, 0);
                    ibz_set(&gamma1.denom, 1);
                    ibz_set(&gamma2.denom, 1);
                    found = 0;
                    transformed = 0;
                    continue;
                }

                assert(quat_lattice_contains(NULL, &small.lattice, mu1));
                assert(quat_lattice_contains(NULL, &small.lattice, mu2));
                quat_alg_conj(&(gamma1), mu1);
                ibz_mul(&gamma1.denom, &gamma1.denom, &n);
                quat_alg_mul(theta, mu2, &gamma1, alg);
                quat_alg_normalize(theta);
                transformed = 1;
            } else {
                if (!(quat_elem_is_odd_norm(theta))) {
                    keep_alpha = 0;
                    ibz_set(&lam, 1);
                    ibz_vec_4_set(&gamma1.coord, 0, 0, 0, 0);
                    ibz_vec_4_set(&gamma2.coord, 0, 0, 0, 0);
                    ibz_set(&gamma1.denom, 1);
                    ibz_set(&gamma2.denom, 1);
                    found = 0;
                    transformed = 0;
                    continue;
                }
            }
            found = 1;
            // num_loops_total = counter;
            break;
        }
    }

    quat_alg_elem_finalize(&gamma1);
    quat_alg_elem_finalize(&gamma2);
    ibz_finalize(&n);
    ibz_finalize(&m);
    ibz_finalize(&tmp);
    ibz_finalize(&x);
    ibz_finalize(&lam);
    ibz_finalize(&two_e);
    ibz_finalize(&norm_n);
    ibz_finalize(&norm_d);
    ibz_finalize(&a_alpha);
    ibz_finalize(&b_alpha);
    ibz_finalize(&alpha0norm);
    ibz_finalize(&psqrt);
    ibz_finalize(&L_det);
    quat_left_ideal_finalize(&small);
    quat_alg_elem_finalize(&alpha);
    quat_alg_elem_finalize(&alpha_0);
    quat_alg_elem_finalize(&encoding);
    quat_alg_elem_finalize(&temp);
    ibz_finalize(&T);
    ibz_vec_2_finalize(&v_target);
    ibz_vec_2_finalize(&v_close);
    ibz_vec_2_finalize(&v_diff);
    ibz_mat_2x2_finalize(&L);
    ibz_mat_2x2_finalize(&L_red);
    ibz_mat_2x2_finalize(&L_inv);
    ibz_mat_4x4_finalize(&gram);

    return (found + found * transformed * 2);
}


////// qlapoty
int
ibz_mod8(const ibz_t *n)
{
    return (ibz_get(n) % 8);
}

void
ibz_mul_2exp(ibz_t *prod, const ibz_t *ipt, int exp)
{
    mpz_mul_2exp(*prod, *ipt, exp);
}

void
ibz_sum_two_squares(ibz_t *sum, const ibz_t *a, const ibz_t *b)
{
    ibz_t prod;
    ibz_init(&prod);
    ibz_mul(&prod, b, b);
    ibz_mul(sum, a, a);
    ibz_add(sum, sum, &prod);
    ibz_finalize(&prod);
}

// gaussian helpers
void
ibz_vec_2_gaussian_mul(ibz_vec_2_t *prod, const ibz_vec_2_t *a, const ibz_vec_2_t *b)
{
    ibz_t r, s, t;
    ibz_init(&r);
    ibz_init(&s);
    ibz_init(&t);
    ibz_add(&r, &(*a)[0], &(*a)[1]);
    ibz_sub(&s, &(*a)[0], &(*a)[1]);
    ibz_add(&t, &(*b)[0], &(*b)[1]);
    ibz_mul(&r, &r, &(*b)[0]);
    ibz_mul(&s, &s, &(*b)[1]);
    ibz_mul(&t, &t, &(*a)[1]);
    ibz_sub(&(*prod)[0], &r, &t);
    ibz_add(&(*prod)[1], &t, &s);
    ibz_finalize(&t);
    ibz_finalize(&s);
    ibz_finalize(&r);
}
void
ibz_vec_2_add(ibz_vec_2_t *sum, const ibz_vec_2_t *a, const ibz_vec_2_t *b)
{
    ibz_add(&(*sum)[0], &(*a)[0], &(*b)[0]);
    ibz_add(&(*sum)[1], &(*a)[1], &(*b)[1]);
}
void
ibz_vec_2_sub(ibz_vec_2_t *diff, const ibz_vec_2_t *a, const ibz_vec_2_t *b)
{
    ibz_sub(&(*diff)[0], &(*a)[0], &(*b)[0]);
    ibz_sub(&(*diff)[1], &(*a)[1], &(*b)[1]);
}
void
ibz_vec_2_copy(ibz_vec_2_t *copy, const ibz_vec_2_t *copied)
{
    ibz_copy(&(*copy)[0], &(*copied)[0]);
    ibz_copy(&(*copy)[1], &(*copied)[1]);
}

int
ibz_vec_2_is_zero(const ibz_vec_2_t *x)
{
    return (ibz_is_zero(&(*x)[0]) && ibz_is_zero(&(*x)[1]));
}

int
ibz_vec_2_gaussian_is_unit(const ibz_vec_2_t *x)
{
    int res;
    ibz_t n;
    ibz_init(&n);
    ibz_sum_two_squares(&n, &(*x)[0], &(*x)[1]);
    res = ibz_is_one(&n);
    ibz_finalize(&n);
    return (res);
}
void
ibz_vec_2_gaussian_euclidean_division(ibz_vec_2_t *q, ibz_vec_2_t *r, const ibz_vec_2_t *a, const ibz_vec_2_t *b)
{
    ibz_t n, tmp, sum;
    ibz_vec_2_t tmp_q, prod, tmp_r;
    ibz_init(&tmp);
    ibz_init(&sum);
    ibz_init(&n);
    ibz_vec_2_init(&tmp_q);
    ibz_vec_2_init(&tmp_r);
    ibz_vec_2_init(&prod);
    ibz_sum_two_squares(&n, &(*b)[0], &(*b)[1]);
    ibz_vec_2_copy(&tmp_q, b);
    ibz_neg(&tmp_q[1], &tmp_q[1]);
    ibz_vec_2_gaussian_mul(&tmp_q, &tmp_q, a);
    ibz_rounded_div(&tmp_q[1], &tmp_q[1], &n);
    ibz_rounded_div(&tmp_q[0], &tmp_q[0], &n);
    ibz_vec_2_gaussian_mul(&prod, &tmp_q, b);
    ibz_vec_2_sub(&tmp_r, a, &prod);
#ifndef NDEBUG
    ibz_sum_two_squares(&tmp, &(tmp_r)[0], &(tmp_r)[1]);
    assert(ibz_cmp(&tmp, &n) < 0);
    ibz_vec_2_gaussian_mul(&prod, &tmp_q, b);
    ibz_vec_2_add(&prod, &prod, &tmp_r);
    assert(ibz_cmp(&prod[0], &(*a)[0]) == 0);
    assert(ibz_cmp(&prod[1], &(*a)[1]) == 0);
#endif
    if (r != NULL)
        ibz_vec_2_copy(r, &tmp_r);
    if (q != NULL)
        ibz_vec_2_copy(q, &tmp_q);
    ibz_finalize(&tmp);
    ibz_finalize(&sum);
    ibz_finalize(&n);
    ibz_vec_2_finalize(&tmp_q);
    ibz_vec_2_finalize(&tmp_r);
    ibz_vec_2_finalize(&prod);
}
void
ibz_vec_2_gaussian_gcd(ibz_vec_2_t *gcd, const ibz_vec_2_t *a, const ibz_vec_2_t *b)
{
    ibz_vec_2_t q, r;
    ibz_vec_2_init(&q);
    ibz_vec_2_init(&r);
    ibz_vec_2_copy(&r, b);
    ibz_vec_2_copy(&q, a);
    if (ibz_vec_2_is_zero(a)) {
        ibz_vec_2_copy(&r, a);
        ibz_vec_2_copy(&q, b);
    }
    while (!ibz_vec_2_is_zero(&r)) {
        ibz_vec_2_gaussian_euclidean_division(NULL, &q, &q, &r);
        ibz_swap(&q[0], &r[0]);
        ibz_swap(&q[1], &r[1]);
    }
    ibz_vec_2_copy(gcd, &q);
#ifndef NDEBUG
    ibz_vec_2_gaussian_euclidean_division(NULL, &q, a, gcd);
    assert(ibz_vec_2_is_zero(&q));
    ibz_vec_2_gaussian_euclidean_division(NULL, &r, b, gcd);
    assert(ibz_vec_2_is_zero(&r));
#endif
    ibz_vec_2_finalize(&q);
    ibz_vec_2_finalize(&r);
}

void
quat_qlapoty_gen_to_dim2_lattice(ibz_mat_2x2_t *lat, ibz_vec_2_t *target, const quat_alg_elem_t *gen, const ibz_t *N)
{
    ibz_vec_2_t uv, ab, kl;
    ibz_t c, gcd, N4div, N4, r;
    ibz_vec_2_init(&ab);
    ibz_vec_2_init(&uv);
    ibz_vec_2_init(&kl);
    ibz_init(&c);
    ibz_init(&gcd);
    ibz_init(&N4div);
    ibz_init(&N4);
    ibz_init(&r);
    assert(ibz_cmp(&gen->denom, &ibz_const_two) == 0);
    ibz_copy(&ab[0], &gen->coord[0]);
    ibz_copy(&ab[1], &gen->coord[1]);
    ibz_mul_2exp(&N4, N, 2);
#ifndef NDEBUG
    ibz_gcd(&gcd, &ab[0], &ab[1]);
    ibz_gcd(&gcd, &gcd, &N4);
    assert(ibz_is_one(&gcd));
#endif
    // set gcd and Bézout coefficients
    ibz_gcd(&gcd, &ab[0], &N4);
    ibz_div(&N4div, &r, &N4, &gcd);
    assert(ibz_is_zero(&r));
    ibz_div(&c, &r, &ab[0], &gcd);
    assert(ibz_is_zero(&r));
    // compute c
    ibz_invmod(&c, &c, &N4div);
    ibz_mul(&c, &c, &ab[1]);
    ibz_neg(&c, &c);
    ibz_mod(&c, &c, &N4div);
    // copy into basis
    ibz_copy(&(*lat)[0][0], &c);
    ibz_copy(&(*lat)[1][0], &gcd);
    ibz_copy(&(*lat)[0][1], &N4div);
    ibz_copy(&(*lat)[1][1], &ibz_const_zero);
    // set target
    if (target != NULL) {
        ibz_xgcd(&gcd, &uv[0], &uv[1], &ab[0], &ab[1]);
        ibz_xgcd(&gcd, &kl[0], &kl[1], &gcd, &N4);
        assert(ibz_is_one(&gcd));
        ibz_xgcd(&gcd, &uv[0], &uv[1], &ab[0], &ab[1]);
        ibz_mul(&(*target)[0], &kl[0], &uv[0]);
        ibz_mul(&(*target)[1], &kl[0], &uv[1]);
        ibz_mod(&(*target)[0], &(*target)[0], &N4);
        ibz_mod(&(*target)[1], &(*target)[1], &N4);
    }
    ibz_vec_2_finalize(&ab);
    ibz_vec_2_finalize(&uv);
    ibz_vec_2_finalize(&kl);
    ibz_finalize(&c);
    ibz_finalize(&gcd);
    ibz_finalize(&N4div);
    ibz_finalize(&N4);
    ibz_finalize(&r);
}

void
quat_qlapoty_get_short_basis(ibz_mat_2x2_t *basis,
                             ibz_vec_2_t *target,
                             const ibz_vec_2_t *st,
                             const quat_alg_elem_t *alpha0,
                             const ibz_t *N)
{
    ibz_t gcd, N4, tmp, r, Nmod, c;
    ibz_vec_2_t uv, kl, omega;
    ibz_mat_2x2_t xy;
    ibz_vec_2_init(&uv);
    ibz_vec_2_init(&kl);
    ibz_vec_2_init(&omega);
    ibz_init(&gcd);
    ibz_init(&N4);
    ibz_init(&tmp);
    ibz_init(&r);
    ibz_init(&Nmod);
    ibz_init(&c);
    ibz_mat_2x2_init(&xy);
    assert(ibz_cmp(&alpha0->denom, &ibz_const_two) == 0);

    // setup
    ibz_mul_2exp(&N4, N, 2);
    ibz_vec_2_copy(&omega, st);
    quat_qlapoty_gen_to_dim2_lattice(basis, target, alpha0, N);
    assert(quat_dim2_lattice_contains(basis, &omega[0], &omega[1]));
    // copies
    ibz_copy(&c, &(*basis)[0][0]);
    ibz_copy(&gcd, &(*basis)[1][0]);
    ibz_copy(&Nmod, &(*basis)[0][1]);
#ifndef NDEBUG
    ibz_mul(&tmp, &gcd, &Nmod);
    assert(ibz_cmp(&tmp, &N4) == 0);
#endif
    // a1
    ibz_vec_2_copy(&omega, st);
    ibz_div(&xy[0][0], &r, &omega[1], &gcd);
    assert(ibz_is_zero(&r));
    // a2
    ibz_mul(&xy[1][0], &c, &omega[1]);
    ibz_mul(&tmp, &gcd, &omega[0]);
    ibz_sub(&xy[1][0], &tmp, &xy[1][0]);
    ibz_div(&xy[1][0], &r, &xy[1][0], &N4);
    assert(ibz_is_zero(&r));
    // a3 a4
    ibz_xgcd(&r, &xy[1][1], &xy[0][1], &xy[0][0], &xy[1][0]);
    assert(ibz_is_one(&r));
    ibz_neg(&xy[0][1], &xy[0][1]);
    // basis
    ibz_2x2_mul(basis, basis, &xy);
#ifndef NDEBUG
    // check proposed algo
    assert(ibz_cmp(&(*basis)[0][0], &omega[0]) == 0);
    assert(ibz_cmp(&(*basis)[1][0], &omega[1]) == 0);
    ibz_mul(&tmp, &xy[0][1], &c);
    ibz_mul(&r, &xy[1][1], &Nmod);
    ibz_add(&tmp, &tmp, &r);
    assert(ibz_cmp(&(*basis)[0][1], &tmp) == 0);
    ibz_mul(&tmp, &xy[0][1], &gcd);
    assert(ibz_cmp(&(*basis)[1][1], &tmp) == 0);
#endif
    quat_dim2_lattice_short_basis(basis, basis, &ibz_const_one);
    ibz_finalize(&gcd);
    ibz_finalize(&N4);
    ibz_finalize(&tmp);
    ibz_finalize(&r);
    ibz_finalize(&Nmod);
    ibz_finalize(&c);
    ibz_vec_2_finalize(&uv);
    ibz_vec_2_finalize(&kl);
    ibz_vec_2_finalize(&omega);
    ibz_mat_2x2_finalize(&xy);
}

// f and big interval needs to be of length 2
// 2d lattice Li should already be reduced
// C=1
// gen canonical generator (x+yi+j)/2
int
quat_qlapoty_initialize_bounds(ibz_vec_2_t *ab,
                               ibz_vec_2_t *omega,
                               quat_alg_elem_t *gen,
                               uint32_t *f_and_big_interval,
                               const ibz_t *N,
                               const ibz_mat_2x2_t *Li,
                               const quat_alg_t *alg,
                               const ibz_t *twoe)
{
    ibz_vec_2_t B;
    quat_alg_elem_t fivevec;
    ibz_t li, Ba, tmp, gcd, fivef, r, cst, cst2;
    int ok = 0;
    ibz_vec_2_init(&B);
    quat_alg_elem_init(&fivevec);
    ibz_init(&li);
    ibz_init(&Ba);
    ibz_init(&gcd);
    ibz_init(&cst);
    ibz_init(&cst2);
    ibz_init(&tmp);
    ibz_init(&r);
    ibz_init(&fivef);
    // get min
    ibz_copy(&(*omega)[0], &(*Li)[0][0]);
    ibz_copy(&(*omega)[1], &(*Li)[1][0]);
    ibz_sum_two_squares(&li, &(*omega)[0], &(*omega)[1]);
    // B2 (N ) ← sqrt(8N) , B1 (N ) = 0.8sqrtN
    // but square everything
    ibz_mul_2exp(&B[1], N, 3);
    ibz_set(&tmp, 64);
    ibz_mul(&B[0], N, &tmp);
    ibz_set(&tmp, 100);
    ibz_div(&B[0], &r, &B[0], &tmp);
    // Ba =  sqrt((2 ^(e+3) N ) - 16N 2 )/p), here also compute square
    ibz_mul_2exp(&Ba, twoe, 3);
    ibz_mul_2exp(&tmp, N, 4);
    ibz_sub(&Ba, &Ba, &tmp);
    ibz_mul(&Ba, &Ba, N);
    ibz_div(&Ba, &r, &Ba, &alg->p);
    // if li ≤ 2B[0] /Ba: Return ⊥
    ibz_mul(&r, &li, &Ba);
    ibz_mul_2exp(&tmp, &B[0], 1);
    if (ibz_cmp(&r, &tmp) <= 0) {
        ok = 0;
        goto fin;
    }
    // a = B[0]/li
    ibz_div(&(*ab)[0], &r, &B[0], &li);
    if (!ibz_is_zero(&r))
        ibz_add(&(*ab)[0], &(*ab)[0], &ibz_const_one);
    // b = min Ba/2, B[1]/li
    ibz_div(&(*ab)[1], &r, &B[1], &li);
    ibz_div_2exp(&tmp, &Ba, 1);
    if (ibz_cmp(&(*ab)[1], &tmp) > 0)
        ibz_copy(&(*ab)[1], &tmp);
    // if
    // assume log in base 2, and C=1
    ibz_sub(&tmp, &(*ab)[1], &(*ab)[0]);
    ibz_set(&cst, ibz_bitsize(N));
    ibz_mul(&cst, &cst, &cst);
    f_and_big_interval[1] = (ibz_cmp(&tmp, &cst) > 0);
    if (!f_and_big_interval[1]) {
        ibz_gcd(&gcd, &(*omega)[0], &(*omega)[1]);
        f_and_big_interval[0] = 0;
        ibz_set(&tmp, 5);
        ibz_set(&fivef, 1);
        ibz_mul(&cst2, &cst, &cst);
        // is there a faster way to do this pecisely enough?
        while (ibz_divides(&gcd, &tmp) && (ibz_cmp(&fivef, &cst2) < 0)) {
            ibz_div(&gcd, &r, &gcd, &tmp);
            assert(ibz_is_zero(&r));
            ibz_mul(&fivef, &fivef, &tmp);
            f_and_big_interval[0] = f_and_big_interval[0] + 1;
        }
        if (f_and_big_interval[0] != 0) {
            ibz_div(&(*omega)[0], &r, &(*omega)[0], &fivef);
            assert(ibz_is_zero(&r));
            ibz_div(&(*omega)[1], &r, &(*omega)[1], &fivef);
            assert(ibz_is_zero(&r));
            ibz_mul_2exp(&r, &(*omega)[0], 1);
            ibz_sub(&r, &r, &(*omega)[1]);
            // r needed for test in 2 lines
            quat_alg_elem_set(&fivevec, 1, 2, 1, 0, 0);
            // tmp is still 5
            if (ibz_divides(&r, &tmp)) {
                ibz_neg(&fivevec.coord[1], &fivevec.coord[1]);
            }
            // maybe implement fast exponentiation if needed
            for (uint32_t i = 1; i < f_and_big_interval[0]; i++) {
                ibz_vec_2_gaussian_mul((ibz_vec_2_t *)&fivevec.coord[0],
                                       (ibz_vec_2_t *)&fivevec.coord[0],
                                       (ibz_vec_2_t *)&fivevec.coord[0]);
            }
            ibz_vec_2_gaussian_mul(omega, (ibz_vec_2_t *)&fivevec.coord[0], omega);
            quat_alg_mul(gen, &fivevec, gen, alg);
        }
        // set ouput
        ibz_mul(&(*ab)[0], &fivef, &(*ab)[0]);
        ibz_mul(&(*ab)[1], &fivef, &(*ab)[1]);
        ibz_sub(&tmp, &(*ab)[1], &(*ab)[0]);
        f_and_big_interval[1] = (ibz_cmp(&tmp, &cst) > 0);
    } else {
        f_and_big_interval[0] = 0;
    }

    ok = 1;
fin:;
    ibz_vec_2_finalize(&B);
    quat_alg_elem_finalize(&fivevec);
    ibz_finalize(&li);
    ibz_finalize(&Ba);
    ibz_finalize(&gcd);
    ibz_finalize(&tmp);
    ibz_finalize(&r);
    ibz_finalize(&fivef);
    ibz_finalize(&cst);
    ibz_finalize(&cst2);
    return (ok);
}

// assume that a, b are squares of what the yare in paper
// big_interval is 1 if b-a>=ClogN^2 as in line 6 loop1
// major mistake to be eliminated: mul by roots of 5 on bounds, not div
void
quat_qlapoty_loop_one(quat_alg_elem_t *alpha0,
                      ibz_vec_2_t *short_solution,
                      const quat_alg_elem_t *gen,
                      const ibz_t *N,
                      const ibz_vec_2_t *absq,
                      const ibz_vec_2_t *omega,
                      const quat_alg_t *alg,
                      uint32_t f,
                      int big_interval)
{
    ibz_vec_2_t k, onebounds, twobounds, gamma, st;
    ibz_t tmp, q, r, sum, fivef, ksq, ngamma, gcdomega;
    quat_alg_elem_t work;
    ibz_init(&tmp);
    ibz_init(&sum);
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&fivef);
    ibz_init(&gcdomega);
    ibz_init(&ksq);
    ibz_init(&ngamma);
    ibz_vec_2_init(&k);
    ibz_vec_2_init(&onebounds);
    ibz_vec_2_init(&twobounds);
    ibz_vec_2_init(&gamma);
    ibz_vec_2_init(&st);
    quat_alg_elem_init(&work);

    // sampling bounds for k1
    ibz_set(&fivef, 5);
    ibz_pow(&fivef, &fivef, f);
    // gcd omega
    ibz_gcd(&gcdomega, &(*omega[0]), &(*omega)[1]);

    // sqrt{5^fa^2/2}
    ibz_div_2exp(&tmp, &(*absq)[0], 1);
    ibz_sqrt_floor(&onebounds[0], &tmp);
    // rounding up
    ibz_mul(&sum, &onebounds[0], &onebounds[0]);
    if (!(ibz_cmp(&sum, &tmp) == 0) && ibz_is_even(&(*absq)[0]))
        ibz_add(&onebounds[0], &onebounds[0], &ibz_const_one);

    // bsqrt{ 5^f} with sqrt rounded up this time
    ibz_sqrt_floor(&onebounds[1], &(*absq)[1]);
    assert(ibz_cmp(&onebounds[0], &onebounds[1]) < 0);

    while (1) {
        // sample k1
        ibz_rand_interval(&k[0], &onebounds[0], &onebounds[1]);
        // sqrt(max (0,(a^2-k1^2)) rounded up
        ibz_mul(&ksq, &k[0], &k[0]);
        if (ibz_cmp(&(*absq)[0], &ksq) > 0) {
            ibz_sub(&tmp, &(*absq)[0], &ksq);
            ibz_sqrt_floor(&twobounds[0], &tmp);
            ibz_mul(&sum, &twobounds[0], &twobounds[0]);
            if (ibz_cmp(&sum, &tmp) < 0)
                ibz_add(&twobounds[0], &twobounds[0], &ibz_const_one);
        } else {
            ibz_set(&twobounds[0], 0);
        }
        // sqrt(b^2-k1^2) rounded down
        ibz_sub(&tmp, &(*absq)[1], &ksq);
        ibz_sqrt_floor(&twobounds[1], &tmp);
        // sample k2
        ibz_rand_interval(&k[1], &twobounds[0], &twobounds[1]);
        // end samping
        if (ibz_is_even(&k[0]) == ibz_is_even(&k[1]))
            continue;
        ibz_gcd(&tmp, &k[1], &k[0]);
        if (!ibz_is_one(&tmp))
            continue;
        ibz_mul(&tmp, &k[1], &k[1]);
        ibz_add(&tmp, &ksq, &tmp);
        if (big_interval) {
            ibz_gcd(&tmp, &tmp, N);
            if (!ibz_is_one(&tmp))
                continue;
        } else {
            ibz_gcd(&tmp, &tmp, &gcdomega);
            if (!ibz_is_one(&tmp))
                continue;
        }
        if (ibz_is_odd(&k[0]))
            ibz_swap(&k[0], &k[1]);
        if (f & 1) {
            ibz_swap(&k[0], &k[1]);
            ibz_neg(&k[0], &k[0]);
        }
        // gaussian gcd (continue line 9)

        ibz_neg(&k[1], &k[1]);
        ibz_vec_2_gaussian_gcd(&gamma, &k, omega);
        ibz_neg(&k[1], &k[1]);
        // ngamma
        ibz_sum_two_squares(&ngamma, &gamma[0], &gamma[1]);
        while (!ibz_is_one(&ngamma)) {
            // update k
            ibz_vec_2_gaussian_mul(&gamma, &gamma, &gamma);
            ibz_vec_2_gaussian_mul(&k, &k, &gamma);
            ibz_div(&k[0], &r, &k[0], &ngamma);
            assert(ibz_is_zero(&r));
            ibz_div(&k[1], &r, &k[1], &ngamma);
            assert(ibz_is_zero(&r));
            // prepare for next test
            ibz_neg(&k[1], &k[1]);
            ibz_vec_2_gaussian_gcd(&gamma, &k, omega);
            ibz_neg(&k[1], &k[1]);
            // ngamma
            ibz_sum_two_squares(&ngamma, &gamma[0], &gamma[1]);
        }
        ibz_vec_2_gaussian_mul(&gamma, &gamma, &gamma);
        ibz_vec_2_gaussian_mul(&k, &k, &gamma);
        ibz_div(&k[0], &r, &k[0], &ngamma);
        assert(ibz_is_zero(&r));
        ibz_div(&k[1], &r, &k[1], &ngamma);
        assert(ibz_is_zero(&r));
        // alpha0
        quat_alg_elem_set(&work, 1, 0, 0, 0, 0);
        ibz_vec_2_copy((ibz_vec_2_t *)&work.coord[0], &k);
        quat_alg_mul(&work, &work, gen, alg);
        assert(ibz_cmp(&work.denom, &ibz_const_two) == 0);
        ibz_gcd(&tmp, &work.coord[0], &work.coord[1]);
        ibz_gcd(&tmp, N, &tmp);
        if (!ibz_is_one(&tmp))
            continue;
        break;
    }
    ibz_vec_2_gaussian_mul(&st, &k, omega);
    quat_alg_elem_copy(alpha0, &work);
    ibz_vec_2_copy(short_solution, &st);

    ibz_finalize(&fivef);
    ibz_finalize(&tmp);
    ibz_finalize(&sum);
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&ksq);
    ibz_finalize(&ngamma);
    ibz_finalize(&gcdomega);
    ibz_vec_2_finalize(&k);
    ibz_vec_2_finalize(&onebounds);
    ibz_vec_2_finalize(&twobounds);
    ibz_vec_2_finalize(&gamma);
    ibz_vec_2_finalize(&st);
    quat_alg_elem_finalize(&work);
}

// For loop 2
//  Rinv is the inverse of R from algo to get alpha0
//  int quat_qlapoty_sample_lambda(ibz_t *l, ibz_t *s, ibz_t *t, ibz_t *sc, ibz_t *tc, const ibz_mat_2x2_t *R, const
//  ibz_mat_2x2_t *Rinv, const ibz_t *Rdet, const ibz_t *aalpha0inv, const quat_alg_elem_t *alpha0, const ibz_t *r,
//  const ibz_t *N, const ibz_cornacchia_extended_params_t *cornacchia_params, const ibz_t *twoe) need alpha0,
//  aalpha0inv,cparams, N, 2e, R, Rinv, Rdet

// Rinv is the inverse of R from algo to get alpha0
int
quat_qlapoty_loop_two(ibz_t *l,
                      ibz_t *s,
                      ibz_t *t,
                      ibz_t *sc,
                      ibz_t *tc,
                      const quat_alg_elem_t *alpha0,
                      const ibz_t *r,
                      const ibz_vec_2_t *st0,
                      const ibz_t *N,
                      const ibz_cornacchia_extended_params_t *cornacchia_params,
                      const ibz_t *twoe)
{
    int found = 0;
    ibz_t lambda, lambdainv, k, tmp, rtwo, twoemdN, twoN, quadN, sum, z, bound, Rdet, deltaN;
    ibz_vec_2_t target, target0, targetl;
    ibz_mat_2x2_t R, Rinv;
    ibz_init(&k);
    ibz_init(&Rdet);
    ibz_init(&tmp);
    ibz_init(&sum);
    ibz_init(&rtwo);
    ibz_init(&twoN);
    ibz_init(&quadN);
    ibz_init(&twoemdN);
    ibz_init(&lambda);
    ibz_init(&lambdainv);
    ibz_init(&z);
    ibz_init(&deltaN);
    ibz_init(&bound);
    ibz_vec_2_init(&target);
    ibz_vec_2_init(&target0);
    ibz_vec_2_init(&targetl);
    ibz_mat_2x2_init(&R);
    ibz_mat_2x2_init(&Rinv);
    // Delta seems opposite to note, check what is right
    assert(ibz_cmp(N, &ibz_const_zero) > 0);
    ibz_set(&deltaN, 1 + 2 * (ibz_is_odd(r) == 0));
    ibz_mul(&deltaN, &deltaN, N);
    // setup lattices
    quat_qlapoty_get_short_basis(&R, &target0, st0, alpha0, N);
    ibz_mat_2x2_inv_with_det_as_denom(&Rinv, &Rdet, &R);
    ibz_neg(&target0[0], &target0[0]);
    ibz_neg(&target0[1], &target0[1]);

    ibz_set(&lambda, 1);
    ibz_neg(&lambda, &lambda);
    ibz_add(&rtwo, r, r);
    ibz_sub(&twoemdN, twoe, &deltaN); // 2^e - delta N
    ibz_add(&twoN, N, N);             // 2N
    ibz_add(&quadN, &twoN, &twoN);    // 4N
    ibz_div(&bound, &z, twoe, r);
    ibz_sqrt_floor(&bound, &bound);

    while (ibz_cmp(&lambda, &bound) <= 0) {
        // get lambda
        ibz_add(&lambda, &lambda, &ibz_const_two);
        ibz_gcd(&tmp, &lambda, &quadN);
        if (!ibz_is_one(&tmp))
            continue;
        // get s,t
        ibz_mul(&tmp, &rtwo, &lambda); // tmp = 2 r lambda
        ibz_mul(&tmp, &tmp, &lambda);  // tmp = 2 r lambda^2
        ibz_sub(&sum, &twoemdN, &tmp); // sum = 2^e-deltaN-2 r lambda^2

        ibz_invmod(&lambdainv, &lambda, &quadN); //  sum = (2^e -deltaN -2r lambda^2 )
        ibz_mul(&tmp, &sum, &lambdainv);         // tmp =  (2^e -deltaN -2r lambda^2 )(lambda)^{-1})
        ibz_mod(&tmp, &tmp, &quadN);             // tmp = (2^e -deltaN -2r lambda^2 )lambdainv mod 4N
        ibz_mul(&targetl[0], &target0[0], &tmp); // -targetl target0 tmp
        ibz_mul(&targetl[1], &target0[1], &tmp);
        ibz_neg(&targetl[0], &targetl[0]);
        ibz_neg(&targetl[1], &targetl[1]);
        ibz_mod(&targetl[0], &targetl[0], &quadN);
        ibz_mod(&targetl[1], &targetl[1], &quadN);
        ibz_mat_2x2_eval(&target, &Rinv, &targetl);
        // inv is actually Rinv * 2N (matrix determinant), since it is integer.
        // So divide and round
        ibz_rounded_div(&target[0], &target[0], &Rdet);
        ibz_rounded_div(&target[1], &target[1], &Rdet);
        ibz_mat_2x2_eval(&target, &R, &target);
        ibz_sub(&target[0], &targetl[0], &target[0]);
        ibz_sub(&target[1], &targetl[1], &target[1]);

#ifndef NDEBUG
        ibz_t x, y;
        ibz_mat_2x2_t test;
        ibz_mat_2x2_init(&test);
        ibz_init(&x);
        ibz_init(&y);
        // inverse test
        ibz_mat_2x2_inv_with_det_as_denom(&test, &x, &R);
        assert(ibz_cmp(&((Rinv))[0][0], &test[0][0]) == 0);
        assert(ibz_cmp(&((Rinv))[0][1], &test[0][1]) == 0);
        assert(ibz_cmp(&((Rinv))[1][0], &test[1][0]) == 0);
        assert(ibz_cmp(&((Rinv))[1][1], &test[1][1]) == 0);
        assert(ibz_cmp(&Rdet, &x) == 0);
        // inverse test lambda
        ibz_mul(&x, &lambda, &lambdainv);
        ibz_mod(&x, &x, &quadN);
        assert(ibz_is_one(&x));
        // homogeneous test
        ibz_vec_2_t test_v;
        ibz_vec_2_init(&test_v);
        ibz_vec_2_set(&test_v, 5, -7);
        ibz_mat_2x2_eval(&test_v, &R, &test_v);
        ibz_mul(&x, &alpha0->coord[0], &test_v[0]);
        ibz_mul(&y, &alpha0->coord[1], &test_v[1]);
        ibz_add(&x, &x, &y);
        ibz_mul(&x, &x, &lambda);
        ibz_div(&y, &x, &x, &quadN);
        assert(ibz_is_zero(&x));
        ibz_vec_2_finalize(&test_v);

        // equation test
        ibz_mul(&x, &alpha0->coord[0], &target[0]);
        ibz_mul(&y, &alpha0->coord[1], &target[1]);
        ibz_add(&x, &x, &y);
        ibz_mul(&x, &x, &lambda);
        ibz_sub(&x, &sum, &x);
        // ibz_add(&x,&sum,&x);
        ibz_div(&x, &y, &x, &quadN);
        assert(ibz_is_zero(&y));
        ibz_finalize(&x);
        ibz_finalize(&y);
        ibz_mat_2x2_finalize(&test);
#endif

        // Test solution is alternating sign
        if (ibz_is_even(&target[0]) == ibz_is_even(&target[1]))
            continue;
        // get k
        ibz_add(&k, &sum, &deltaN);                   // k = 2^e-deltaN -2 r lambda^2 + deltaN = 2^e-2r lambda^2
        ibz_mul(&tmp, &target[0], &alpha0->coord[0]); // tmp = 2as
        ibz_mul(&sum, &target[1], &alpha0->coord[1]); // sum=2bt
        ibz_add(&sum, &tmp, &sum);                    // sum = 2as+2bt
        ibz_mul(&tmp, &sum, &lambda);                 // tmp = 2as lambda +2bt lambda
        ibz_sub(&k, &k, &tmp);                        // k = 2^e-2 r lambda^2- 2as lambda - 2bt lambda
        ibz_div(&k, &tmp, &k, N);
        assert(ibz_is_zero(&tmp));
        assert(ibz_is_odd(&k));

        // z
        ibz_add(&z, &k, &k);
        ibz_mul(&sum, &target[0], &target[0]);
        ibz_mul(&tmp, &target[1], &target[1]);

        ibz_sub(&z, &z, &tmp);
        ibz_sub(&z, &z, &sum);

        // z size
        if (ibz_cmp(&z, &ibz_const_zero) < 0)
            continue;

        // z 1 mod 4
        assert((ibz_get(&z) & 3) == 1);

        // Odd solutions
        assert((ibz_mod8(r) + ((ibz_mod8(&k) - 1) % 4) / 2) % 2 == 1);

        // Cornacchia
        found = ibz_cornacchia_extended(sc, tc, &z, cornacchia_params);
        if (found)
            break;
    }
    if (found) {
        ibz_copy(s, &target[0]);
        ibz_copy(t, &target[1]);
        ibz_copy(l, &lambda);
#ifndef NDEBUG
        ibz_t tn1, tn2, tt;
        ibz_init(&tt);
        ibz_init(&tn1);
        ibz_init(&tn2);
        ibz_sum_two_squares(&tt, sc, tc);
        assert(ibz_cmp(&z, &tt) == 0);
        ibz_sum_two_squares(&tt, s, t);
        ibz_add(&tt, &z, &tt);
        assert(ibz_is_even(&tt));
        ibz_div_2exp(&tt, &tt, 1);
        assert(ibz_cmp(&k, &tt) == 0);
        ibz_mul(&tt, &k, N);
        // lattice part: tn2 = lambda(2as + 2bt)
        assert(ibz_cmp(&alpha0->denom, &ibz_const_two) == 0);
        ibz_mul(&tn1, &alpha0->coord[0], s);
        ibz_mul(&tn2, &alpha0->coord[1], t);
        ibz_add(&tn2, &tn2, &tn1);
        ibz_mul(&tn2, &tn2, l);
        // other part: tn1 = 2^e-2rlambda^2-tn2
        ibz_mul(&tn1, l, l);
        ibz_mul(&tn1, &tn1, r);
        ibz_add(&tn1, &tn1, &tn1);
        ibz_sub(&tn1, twoe, &tn1);
        ibz_sub(&tn1, &tn1, &tn2);
        assert(ibz_cmp(&tt, &tn1) == 0);

        ibz_finalize(&tt);
        ibz_finalize(&tn1);
        ibz_finalize(&tn2);
#endif
    }

    ibz_finalize(&Rdet);
    ibz_finalize(&k);
    ibz_finalize(&tmp);
    ibz_finalize(&sum);
    ibz_finalize(&lambda);
    ibz_finalize(&lambdainv);
    ibz_finalize(&z);
    ibz_finalize(&rtwo);
    ibz_finalize(&twoN);
    ibz_finalize(&quadN);
    ibz_finalize(&twoemdN);
    ibz_finalize(&bound);
    ibz_finalize(&deltaN);
    ibz_vec_2_finalize(&target);
    ibz_vec_2_finalize(&target0);
    ibz_vec_2_finalize(&targetl);
    ibz_mat_2x2_finalize(&R);
    ibz_mat_2x2_finalize(&Rinv);
    return (found);
}

int
quat_qlapoty(quat_alg_elem_t *mu1,
             quat_alg_elem_t *mu2,
             quat_alg_elem_t *theta,
             quat_alg_elem_t *smallest,
             const quat_left_ideal_t *lideal,
             const quat_alg_t *alg,
             int max_counter_alpha,
             int gen_sampling_bound_bits,
             int two_power,
             const ibz_cornacchia_extended_params_t *cornacchia_params)
{
    int found = 0;
    ibz_t s, t, r, sc, tc, lambda, twoe;
    quat_alg_elem_t gen, temp, alpha0;
    quat_left_ideal_t small;
    ibz_vec_2_t ab, omega, st0;
    ibz_mat_2x2_t Li;
    uint32_t f_and_big_interval[2] = { 0, 0 };
    ibz_init(&s);
    ibz_init(&t);
    ibz_init(&sc);
    ibz_init(&tc);
    ibz_init(&r);
    ibz_init(&lambda);
    ibz_init(&twoe);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&temp);
    quat_alg_elem_init(&alpha0);
    quat_left_ideal_init(&small);
    ibz_vec_2_init(&ab);
    ibz_vec_2_init(&omega);
    ibz_vec_2_init(&st0);
    ibz_mat_2x2_init(&Li);

    quat_lideal_shortest_equivalent(&small, smallest, lideal, alg);
    quat_lattice_hnf(&small.lattice);
    // compute gen
    ibz_copy(&gen.coord[0], &small.lattice.basis[0][2]);
    ibz_copy(&gen.coord[1], &small.lattice.basis[1][2]);
    ibz_copy(&gen.coord[2], &ibz_const_one);
    ibz_copy(&gen.coord[3], &ibz_const_zero);
    ibz_copy(&gen.denom, &small.lattice.denom);
    // compute 2d basis and minima
    quat_qlapoty_gen_to_dim2_lattice(&Li, NULL, &gen, &small.norm);
    quat_dim2_lattice_short_basis(&Li, &Li, &ibz_const_one);
#ifndef NDEBUG
    // test omega is in lattice, and lattice is correct
    ibz_t tmp1, sum1;
    ibz_init(&tmp1);
    ibz_init(&sum1);
    ibz_mul(&tmp1, &omega[0], &gen.coord[0]);
    ibz_mul(&sum1, &omega[1], &gen.coord[1]);
    ibz_add(&sum1, &tmp1, &sum1);
    ibz_mul_2exp(&tmp1, &small.norm, 2);
    ibz_mul(&sum1, &sum1, &tmp1);
    assert(ibz_is_zero(&sum1));
    ibz_finalize(&tmp1);
    ibz_finalize(&sum1);
#endif
    // bound adjustment
    ibz_mul_2exp(&twoe, &ibz_const_one, two_power);
    found = quat_qlapoty_initialize_bounds(&ab, &omega, &gen, f_and_big_interval, &small.norm, &Li, alg, &twoe);
    if (!found)
        goto fin;
#ifndef NDEBUG
    quat_qlapoty_gen_to_dim2_lattice(&Li, NULL, &gen, &small.norm);
    assert(quat_dim2_lattice_contains(&Li, &omega[0], &omega[1]));
#endif
    // find alpha0
    quat_qlapoty_loop_one(
        &alpha0, &st0, &gen, &small.norm, &ab, &omega, alg, f_and_big_interval[0], f_and_big_interval[1]);
    quat_alg_norm(&r, &t, &alpha0, alg);
    assert(ibz_is_one(&t));
    ibz_div(&r, &t, &r, &small.norm);
    assert(ibz_is_zero(&t));
    // find lambda
    found = quat_qlapoty_loop_two(&lambda, &s, &t, &sc, &tc, &alpha0, &r, &st0, &small.norm, cornacchia_params, &twoe);
    if (!found)
        goto fin;

    if (ibz_is_even(&s) != ibz_is_even(&sc))
        ibz_swap(&sc, &tc);
    assert(ibz_is_even(&s) == ibz_is_even(&sc));
    assert(ibz_is_even(&t) == ibz_is_even(&tc));

    // set output
    quat_alg_elem_set(mu1, 2, 0, 0, 0, 0);
    quat_alg_elem_set(mu2, 2, 0, 0, 0, 0);
    ibz_add(&mu1->coord[0], &s, &sc);
    ibz_add(&mu1->coord[1], &t, &tc);
    ibz_sub(&mu2->coord[0], &s, &sc);
    ibz_sub(&mu2->coord[1], &t, &tc);
    quat_alg_elem_mul_by_scalar(&temp, &lambda, &alpha0);
    quat_alg_elem_mul_by_scalar(mu1, &small.norm, mu1);
    quat_alg_elem_mul_by_scalar(mu2, &small.norm, mu2);
    quat_alg_add(mu1, mu1, &temp);
    quat_alg_add(mu2, mu2, &temp);
#ifndef NDEBUG
    ibz_t tn1, tn2, tt;
    ibz_init(&tt);
    ibz_init(&tn1);
    ibz_init(&tn2);
    quat_alg_norm(&tn1, &tt, mu1, alg);
    assert(ibz_is_one(&tt));
    quat_alg_norm(&tn2, &tt, mu2, alg);
    assert(ibz_is_one(&tt));
    ibz_add(&tn1, &tn1, &tn2);
    ibz_div(&tn2, &tt, &tn1, &small.norm);
    assert(ibz_is_zero(&tt));
    assert(ibz_cmp(&tn2, &twoe) == 0);
    ibz_finalize(&tt);
    ibz_finalize(&tn1);
    ibz_finalize(&tn2);
#endif
    quat_alg_conj(&temp, mu1);
    ibz_mul(&temp.denom, &temp.denom, &small.norm);
    quat_alg_mul(theta, mu2, &temp, alg);
fin:;

    ibz_finalize(&s);
    ibz_finalize(&t);
    ibz_finalize(&sc);
    ibz_finalize(&tc);
    ibz_finalize(&r);
    ibz_finalize(&lambda);
    ibz_finalize(&twoe);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&temp);
    quat_alg_elem_finalize(&alpha0);
    quat_left_ideal_finalize(&small);
    ibz_vec_2_finalize(&ab);
    ibz_vec_2_finalize(&omega);
    ibz_vec_2_finalize(&st0);
    ibz_mat_2x2_finalize(&Li);
    return (found);
}
