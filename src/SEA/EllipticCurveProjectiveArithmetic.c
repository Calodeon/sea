#include "EllipticCurveArithmetic.h"
#include "EllipticCurveProjectiveArithmetic.h"
#include "ModularPolynomial.h"
#include "DivisionPolynomial.h"
#include "toolbox.h"

#include <stdio.h>
#include <string.h>

/*************************
 *   private functions   *
 *************************/

static void
_ECProj_add_wiki(ECProjPoint * P3, const ECProjPoint * P1,
                 const ECProjPoint * P2, const EllipticCurve * E,
                 fmpz_t * lambda)
{
    /* lambda[0] = tmp */
    fmpz_t *tmp = lambda + 0;

    /* lambda[1] = Y_2 * Z_1 */
    fmpz_mul(lambda[1], P2->y, P1->z);
    fmpz_mod(lambda[1], lambda[1], E->p);

    /* lambda[2] = Y_1 * Z_2 */
    fmpz_mul(lambda[2], P1->y, P2->z);
    fmpz_mod(lambda[2], lambda[2], E->p);

    /* lambda[3] = X_2 * Z_1 */
    fmpz_mul(lambda[3], P2->x, P1->z);
    fmpz_mod(lambda[3], lambda[3], E->p);

    /* lambda[4] = X_1 * Z_2 */
    fmpz_mul(lambda[4], P1->x, P2->z);
    fmpz_mod(lambda[4], lambda[4], E->p);

    /* since P1 != P2, lambda[3] = lambda[4] implies P1 = -P2 */
    /* TODO check this sooner */
    if (fmpz_equal(lambda[3], lambda[4]))
    {
        P3->is_zero = 1;
        return;
    }

    /* lambda[5] = lambda[1] - lambda[2] */
    fmpz_sub(lambda[5], lambda[1], lambda[2]);
    fmpz_mod(lambda[5], lambda[5], E->p);

    /* lambda[6] = lambda[3] - lambda[4] */
    fmpz_sub(lambda[6], lambda[3], lambda[4]);
    fmpz_mod(lambda[6], lambda[6], E->p);

    /* lambda[7] = Z_1 * Z_2 */
    fmpz_mul(lambda[7], P1->z, P2->z);
    fmpz_mod(lambda[7], lambda[7], E->p);

    /* lambda[8] = lambda[5]^2 * lambda[7] */
    fmpz_powm_ui(lambda[8], lambda[5], 2, E->p);
    fmpz_mul(lambda[8], lambda[8], lambda[7]);
    fmpz_mod(lambda[8], lambda[8], E->p);

    /* lambda[8] -= lambda[6]^3 */
    fmpz_powm_ui(*tmp, lambda[6], 3, E->p);
    fmpz_sub(lambda[8], lambda[8], *tmp);
    fmpz_mod(lambda[8], lambda[8], E->p);

    /* lambda[8] -= 2*lambda[6]^2 * lambda[2] */
    fmpz_powm_ui(*tmp, lambda[6], 2, E->p);
    fmpz_mul_ui(*tmp, *tmp, 2);
    fmpz_mul(*tmp, *tmp, lambda[4]);
    fmpz_sub(lambda[8], lambda[8], *tmp);
    fmpz_mod(lambda[8], lambda[8], E->p);

    /* X3 = lambda[6] * lambda[8] */
    fmpz_mul(P3->x, lambda[6], lambda[8]);
    fmpz_mod(P3->x, P3->x, E->p);

    /* Y3 = lambda[5] * (lambda[6]^2 * lambda[4] - lambda[8]) */
    fmpz_powm_ui(P3->y, lambda[6], 2, E->p);
    fmpz_mul(P3->y, P3->y, lambda[4]);
    fmpz_sub(P3->y, P3->y, lambda[8]);
    fmpz_mul(P3->y, P3->y, lambda[5]);
    fmpz_mod(P3->y, P3->y, E->p);

    /* Y3 -= lambda[6]^3 * lambda[2] */
    fmpz_powm_ui(*tmp, lambda[6], 3, E->p);
    fmpz_mul(*tmp, *tmp, lambda[2]);
    fmpz_sub(P3->y, P3->y, *tmp);
    fmpz_mod(P3->y, P3->y, E->p);

    /* Z3 = lambda[6]^3 * lambda[7] */
    fmpz_powm_ui(P3->z, lambda[6], 3, E->p);
    fmpz_mul(P3->z, P3->z, lambda[7]);
    fmpz_mod(P3->z, P3->z, E->p);
}

#ifdef ECPROJ_BLAKE_IMPL
static void
_ECProj_add_Blake(ECProjPoint * P3, const ECProjPoint * P1,
                  const ECProjPoint * P2, const EllipticCurve * E,
                  fmpz_t * lambda)
{
    /* lambda[0] = tmp */
    fmpz_t *tmp = lambda + 0;

    /* lambda[1] = X_1 * Z_2^2 */
    fmpz_powm_ui(lambda[1], P2->z, 2, E->p);
    fmpz_mul(lambda[1], lambda[1], P1->x);
    fmpz_mod(lambda[1], lambda[1], E->p);

    /* lambda[2] = X_2 * Z_1^2 */
    fmpz_powm_ui(lambda[2], P1->z, 2, E->p);
    fmpz_mul(lambda[2], lambda[2], P2->x);
    fmpz_mod(lambda[2], lambda[2], E->p);

    /* lambda[3] = lambda[1] - lambda[2] */
    fmpz_sub(lambda[3], lambda[1], lambda[2]);
    fmpz_mod(lambda[3], lambda[3], E->p);

    /* since P1 != P2, lambda[3] = 0 implies P1 = -P2 */
    if (fmpz_equal_ui(lambda[3], 0))
    {
        P3->is_zero = 1;
        return;
    }

    /* lambda[4] = Y_1 * Z_2^3 */
    fmpz_powm_ui(lambda[4], P2->z, 3, E->p);
    fmpz_mul(lambda[4], lambda[4], P1->y);
    fmpz_mod(lambda[4], lambda[4], E->p);

    /* lambda[5] = Y_2 * Z_1^3 */
    fmpz_powm_ui(lambda[5], P1->z, 3, E->p);
    fmpz_mul(lambda[5], lambda[5], P2->y);
    fmpz_mod(lambda[5], lambda[5], E->p);

    /* lambda[6] = lambda[4] - lambda[5] */
    fmpz_sub(lambda[6], lambda[4], lambda[5]);
    fmpz_mod(lambda[6], lambda[6], E->p);

    /* lambda[7] = lambda[1] + lambda[2] */
    fmpz_add(lambda[7], lambda[1], lambda[2]);
    fmpz_mod(lambda[7], lambda[7], E->p);

    /* lambda[8] = lambda[4] + lambda[5] */
    fmpz_add(lambda[8], lambda[4], lambda[5]);
    fmpz_mod(lambda[8], lambda[8], E->p);

    /* Z_3 = Z_1 * Z_2 * lambda[3] */
    fmpz_mul(P3->z, P1->z, P2->z);
    fmpz_mul(P3->z, P3->z, lambda[3]);
    fmpz_mod(P3->z, P3->z, E->p);

    /* X_3 = lambda[6]^2 - lambda[7] * lambda[3]^2 */
    fmpz_powm_ui(P3->x, lambda[6], 2, E->p);
    fmpz_powm_ui(*tmp, lambda[3], 2, E->p);
    fmpz_mul(*tmp, *tmp, lambda[7]);
    fmpz_sub(P3->x, P3->x, *tmp);
    fmpz_mod(P3->x, P3->x, E->p);

    /* lambda[9] = lambda[7] * lambda[3]^2 - 2*X_3 */
    fmpz_powm_ui(lambda[9], lambda[3], 2, E->p);
    fmpz_mul(lambda[9], lambda[9], lambda[7]);
    fmpz_mul_ui(*tmp, P3->x, 2);
    fmpz_sub(lambda[9], lambda[9], *tmp);
    fmpz_mod(lambda[9], lambda[9], E->p);

    /* Y_3 = (lambda[9] * lambda[6] - lambda[8] * lambda[3]^3)/2 */
    fmpz_mul(P3->y, lambda[9], lambda[6]);
    fmpz_powm_ui(*tmp, lambda[3], 3, E->p);
    fmpz_mul(*tmp, *tmp, lambda[8]);
    fmpz_sub(P3->y, P3->y, *tmp);
    fmpz_mod(P3->y, P3->y, E->p);
    fmpz_divexact_ui(P3->y, P3->y, 2);
}
#endif /* ECPROJ_BLAKE_IMPL */

static void
_ECProj_double_wiki(ECProjPoint * P2,
                    const ECProjPoint * P1, const EllipticCurve * E,
                    fmpz_t * lambda)
{
    /* tmp = lambda[0] */
    fmpz_t *tmp = lambda + 0;

    /* lambda[1] = 3 * X_1^2 */
    fmpz_powm_ui(lambda[1], P1->x, 2, E->p);
    fmpz_mul_ui(lambda[1], lambda[1], 3);
    fmpz_mod(lambda[1], lambda[1], E->p);

    /* lambda[1] += A * Z_1^2 */
    fmpz_powm_ui(*tmp, P1->z, 2, E->p);
    fmpz_mul(*tmp, *tmp, E->A);
    fmpz_add(lambda[1], lambda[1], *tmp);
    fmpz_mod(lambda[1], lambda[1], E->p);

    /* lambda[2] = Y_1 * Z_1 */
    fmpz_mul(lambda[2], P1->y, P1->z);
    fmpz_mod(lambda[2], lambda[2], E->p);

    /* lambda[3] = X_1 * Y_1 * lambda[2] */
    fmpz_mul(lambda[3], P1->x, P1->y);
    fmpz_mul(lambda[3], lambda[3], lambda[2]);
    fmpz_mod(lambda[3], lambda[3], E->p);

    /* lambda[4] = lambda[1]^2 - 8 * lambda[3] */
    fmpz_powm_ui(lambda[4], lambda[1], 2, E->p);
    fmpz_mul_ui(*tmp, lambda[3], 8);
    fmpz_sub(lambda[4], lambda[4], *tmp);
    fmpz_mod(lambda[4], lambda[4], E->p);

    /* X_2 = 2 * lambda[4] * lambda[2] */
    fmpz_mul(P2->x, lambda[4], lambda[2]);
    fmpz_mul_ui(P2->x, P2->x, 2);
    fmpz_mod(P2->x, P2->x, E->p);

    /* Y_2 = 8 * Y_1^2 * lambda[2]^2 */
    fmpz_powm_ui(P2->y, P1->y, 2, E->p);
    fmpz_powm_ui(*tmp, lambda[2], 2, E->p);
    fmpz_mul(P2->y, P2->y, *tmp);
    fmpz_mul_ui(P2->y, P2->y, 8);
    fmpz_mod(P2->y, P2->y, E->p);

    /* Y_2 = lambda[1] * (4 * lambda[3] - lambda[4]) - Y_2 */
    fmpz_mul_ui(*tmp, lambda[3], 4);
    fmpz_sub(*tmp, *tmp, lambda[4]);
    fmpz_mul(*tmp, *tmp, lambda[1]);
    fmpz_sub(P2->y, *tmp, P2->y);
    fmpz_mod(P2->y, P2->y, E->p);

    /* Z_2 = 8 * lambda[2]^3 */
    fmpz_powm_ui(P2->z, lambda[2], 3, E->p);
    fmpz_mul_ui(P2->z, P2->z, 8);
    fmpz_mod(P2->z, P2->z, E->p);

}

#ifdef ECPROJ_BLAKE_IMPL
static void
_ECProj_double_Blake(ECProjPoint * P2, const ECProjPoint * P1,
                     const EllipticCurve * E, fmpz_t * lambda)
{
    /* tmp = lambda[0] */
    fmpz_t *tmp = lambda + 0;

    /* lambda[1] = 3 * X_1^2 */
    fmpz_powm_ui(lambda[1], P1->x, 2, E->p);
    fmpz_mul_ui(lambda[3], lambda[3], 3);
    fmpz_mod(lambda[3], lambda[3], E->p);

    /* lambda[1] += A * Z_1^4 */
    fmpz_powm_ui(*tmp, P1->x, 4, E->p);
    fmpz_mul(*tmp, *tmp, E->A);
    fmpz_add(lambda[1], lambda[1], *tmp);
    fmpz_mod(lambda[1], lambda[1], E->p);

    /* Z_2 = 2 * Y_1 * Z_1 */
    fmpz_mul(P2->z, P1->y, P1->z);
    fmpz_mul_ui(P2->z, P2->z, 2);
    fmpz_mod(P2->z, P2->z, E->p);

    /* lambda[2] = 4 * X_1 * Y_1^2 */
    fmpz_powm_ui(lambda[2], P1->y, 2, E->p);
    fmpz_mul(lambda[2], lambda[2], P1->x);
    fmpz_mul_ui(lambda[2], lambda[2], 4);
    fmpz_mod(lambda[2], lambda[2], E->p);

    /* X_2 = lambda[1]^2 - 2 * lambda[2] */
    fmpz_powm_ui(P2->x, lambda[1], 2, E->p);
    fmpz_sub(P2->x, P2->x, lambda[2]);
    fmpz_sub(P2->x, P2->x, lambda[2]);
    fmpz_mod(P2->x, P2->x, E->p);

    /* lambda[3] = 8 * Y_1^4 */
    fmpz_powm_ui(lambda[3], P1->y, 4, E->p);
    fmpz_mul_ui(lambda[3], lambda[3], 8);
    fmpz_mod(lambda[3], lambda[3], E->p);

    /* Y_2 = lambda[1] * (lambda[2] - X_3) - lambda[3] */
    fmpz_sub(P2->y, lambda[2], P2->x);
    fmpz_mul(P2->y, P2->y, lambda[1]);
    fmpz_sub(P2->y, P2->y, lambda[3]);
    fmpz_mod(P2->y, P2->y, E->p);

}
#endif /* ECPROJ_BLAKE_IMPL */

/*************************
 *   public functions    *
 *************************/

/* Create and alloc a new point */
ECProjPoint *
new_ECProjPoint()
{
    ECProjPoint *P = malloc(sizeof(ECProjPoint));
    P->is_zero = 1;
    fmpz_init(P->x);
    fmpz_init(P->y);
    fmpz_init(P->z);
    return P;
}

/* Create and alloc a copy of a point */
ECProjPoint *
new_ECProjPoint_copy(const ECProjPoint * P1)
{
    ECProjPoint *P = malloc(sizeof(ECProjPoint));
    P->is_zero = P1->is_zero;
    if (!P->is_zero)
    {
        fmpz_init_set(P->x, P1->x);
        fmpz_init_set(P->y, P1->y);
        fmpz_init_set(P->z, P1->z);
    }
    return P;
}

/* free the memory */
void
free_ECProjPoint(ECProjPoint * P)
{
    fmpz_clear(P->x);
    fmpz_clear(P->y);
    fmpz_clear(P->z);
    free(P);
}

/* print a point "(x : y : z)" */
void
print_ECProjPoint(ECProjPoint * P)
{
    if (P->is_zero)
        printf("(O : 1 : 0)");

    else
    {
        flint_printf("(");
        fmpz_print(P->x);
        flint_printf(" : ");
        fmpz_print(P->y);
        flint_printf(" : ");
        fmpz_print(P->z);
        flint_printf(")");
    }
}

/* returns 1 (true) if P == Q, 0 (false) otherwise */
int
EllipticCurve_projective_equals(const ECProjPoint * P, const ECProjPoint * Q,
                                const EllipticCurve * E)
{
    if (P->is_zero && Q->is_zero)
        return 1;               /* true */

    else if (P->is_zero || Q->is_zero)
        return 0;               /* false */

    fmpz_t P_x;
    fmpz_t P_y;
    fmpz_t Q_x;
    fmpz_t Q_y;

    fmpz_init(P_x);
    fmpz_init(P_y);
    fmpz_init(Q_x);
    fmpz_init(Q_y);

    fmpz_mul(P_x, P->x, Q->z);
    fmpz_mod(P_x, P_x, E->p);
    fmpz_mul(P_y, P->y, Q->z);
    fmpz_mod(P_y, P_y, E->p);

    fmpz_mul(Q_x, Q->x, P->z);
    fmpz_mod(Q_x, Q_x, E->p);
    fmpz_mul(Q_y, Q->y, P->z);
    fmpz_mod(Q_y, Q_y, E->p);

    int equals = (fmpz_cmp(P_x, Q_x) == 0) && (fmpz_cmp(P_y, Q_y) == 0);

    fmpz_clear(P_x);
    fmpz_clear(P_y);
    fmpz_clear(Q_x);
    fmpz_clear(Q_y);

    return equals;
}

/*
 * Returns 1 (true) if projective point P is equivalent to affine point a, or 
 * 0 (false) otherwise.
 */
int
EllipticCurve_projective_isEquivalentTo(const ECProjPoint * P,
                                        const ECPoint * a,
                                        const EllipticCurve * E)
{
    if (P->is_zero || fmpz_is_zero(P->z))
        return a->is_zero;

    fmpz_t P_x;
    fmpz_t P_y;

    fmpz_init(P_x);
    fmpz_init(P_y);

    fmpz_mul(P_x, a->x, P->z);
    fmpz_mod(P_x, P_x, E->p);
    fmpz_mul(P_y, a->y, P->z);
    fmpz_mod(P_y, P_y, E->p);

    int equals = (fmpz_cmp(P_x, P->x) == 0) && (fmpz_cmp(P_y, P->y) == 0);

    fmpz_clear(P_x);
    fmpz_clear(P_y);

    return equals;
}

/* checks Y^2 * Z = X^3 + A * X * Z^2 + B * Z^3 */
int
EllipticCurve_projective_isOnCurve(const ECProjPoint * P,
                                   const EllipticCurve * E)
{
    fmpz_t RHS;                 /* Right-Hand Side */
    fmpz_t LHS;                 /* Left-Hand Side */
    fmpz_t acc;                 /* accumulator */

    fmpz_init(LHS);
    fmpz_init(RHS);
    fmpz_init(acc);

    /* LHS = y^2 * z mod p */
    fmpz_powm_ui(LHS, P->y, 2, E->p);
    fmpz_mul(LHS, LHS, P->z);
    fmpz_mod(LHS, LHS, E->p);

    /* RHS = x^3 mod p */
    fmpz_powm_ui(RHS, P->x, 3, E->p);

    /* RHS += z^2 * x * A mod p */
    fmpz_powm_ui(acc, P->z, 2, E->p);
    fmpz_mul(acc, acc, P->x);
    fmpz_mul(acc, acc, E->A);
    fmpz_add(RHS, RHS, acc);
    fmpz_mod(RHS, RHS, E->p);

    /* RHS += z^3 * B mod p */
    fmpz_powm_ui(acc, P->z, 3, E->p);
    fmpz_mul(acc, acc, E->B);
    fmpz_add(RHS, RHS, acc);
    fmpz_mod(RHS, RHS, E->p);

    int equals = fmpz_equal(LHS, RHS);

    fmpz_clear(LHS);
    fmpz_clear(RHS);
    fmpz_clear(acc);

    return equals;
}

/* Y^2 * Z = X^3 + A * X * Z^2 + B * Z^3
 * finds (X,Y) on Y^2 = X^3 + AX + B and picks random Z to multiply X and Y with
 */
void
EllipticCurve_projective_randomPoint(ECProjPoint * P, flint_rand_t state, 
                                     const EllipticCurve * E)
{
    /* finds regular point on curve */
    ECPoint *Q = new_ECPoint();
    EllipticCurve_randomPoint(Q, state, E);

    /* multiplies regular point with random Z != 0 */
    P->is_zero = 0;
    do
    {
        fmpz_randm(P->z, state, E->p);
    } while (fmpz_equal_ui(P->z, 0));

    fmpz_mul(P->x, Q->x, P->z);
    fmpz_mod(P->x, P->x, E->p);
    fmpz_mul(P->y, Q->y, P->z);
    fmpz_mod(P->y, P->y, E->p);

    free_ECPoint(Q);
}

void
EllipticCurve_projective_add(ECProjPoint * P3, const ECProjPoint * P1,
                             const ECProjPoint * P2, const EllipticCurve * E)
{
    slong i;
    ECProjPoint *P4 = new_ECProjPoint();
    fmpz_t *lambda = malloc(ECPROJ_LAMBDA_SIZE * sizeof(fmpz_t));
    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_init(lambda[i]);

    EllipticCurve_projective_add_memory(P3, P1, P2, E, P4, lambda);

    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_clear(lambda[i]);

    free(lambda);
    free_ECProjPoint(P4);
}

void
EllipticCurve_projective_double(ECProjPoint * P2, const ECProjPoint * P1,
                                const EllipticCurve * E)
{
    slong i;
    ECProjPoint *P4 = new_ECProjPoint();
    fmpz_t *lambda = malloc(ECPROJ_LAMBDA_SIZE * sizeof(fmpz_t));
    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_init(lambda[i]);

    EllipticCurve_projective_double_memory(P2, P1, E, P4, lambda);

    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_clear(lambda[i]);

    free(lambda);
    free_ECProjPoint(P4);
}

void
EllipticCurve_projective_add_memory(ECProjPoint * P3, const ECProjPoint * P1,
                                    const ECProjPoint * P2,
                                    const EllipticCurve * E, ECProjPoint * P4,
                                    fmpz_t * lambda)
{
    if (P1->is_zero || fmpz_is_zero(P1->z))
    {
        P3->is_zero = P2->is_zero;
        fmpz_set(P3->x, P2->x);
        fmpz_set(P3->y, P2->y);
        fmpz_set(P3->z, P2->z);
        return;
    }
    else if (P2->is_zero || fmpz_is_zero(P2->z))
    {
        P3->is_zero = P1->is_zero;
        fmpz_set(P3->x, P1->x);
        fmpz_set(P3->y, P1->y);
        fmpz_set(P3->z, P1->z);
        return;
    }

    /* if P1 = P2, then P3 = [2] P1 and doubling is used */
    if (EllipticCurve_projective_equals(P1, P2, E))
    {
        EllipticCurve_projective_double_memory(P3, P1, E, P4, lambda);
        return;
    }

    P4->is_zero = 0;

    /*
     *  The following is a working implementation of
     *  https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Standard_Projective_Coordinates#Point_Addition_.2812M_.2B_2S.29
     */

#ifndef ECPROJ_BLAKE_IMPL

    _ECProj_add_wiki(P4, P1, P2, E, lambda);

#endif /* ECPROJ_BLAKE_IMPL */

    /*
     *  The following is a non-working implementation of
     *  I. Blake, G. Seroussi, N. Smart, Elliptic Curves in Cryptography,
     *      Cambridge University Press, p.59
     */

#ifdef ECPROJ_BLAKE_IMPL

    _ECProj_add_Blake(P4, P1, P2, E, lambda);

#endif /* ECPROJ_BLAKE_IMPL */

    /* Make addition in place */
    fmpz_set(P3->x, P4->x);
    fmpz_set(P3->y, P4->y);
    fmpz_set(P3->z, P4->z);
    P3->is_zero = P4->is_zero;
}

void
EllipticCurve_projective_double_memory(ECProjPoint * P2,
                                       const ECProjPoint * P1,
                                       const EllipticCurve * E,
                                       ECProjPoint * P4, fmpz_t * lambda)
{
    if (P1->is_zero || fmpz_is_zero(P1->z))
    {
        P2->is_zero = 1;
        return;
    }

    P4->is_zero = 0;

    /*
     *  The following is a working implementation of
     *  https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Standard_Projective_Coordinates#Point_Doubling_.287M_.2B_5S_or_7M_.2B_3S.29
     */

#ifndef ECPROJ_BLAKE_IMPL

    _ECProj_double_wiki(P4, P1, E, lambda);

#endif /* ECPROJ_BLAKE_IMPL */

    /*
     *  The following is a non-working implementation of
     *  I. Blake, G. Seroussi, N. Smart, Elliptic Curves in Cryptography,
     *      Cambridge University Press, p.60
     */

#ifdef ECPROJ_BLAKE_IMPL

    _ECProj_double_Blake(P4, P1, E, lambda);

#endif /* ECPROJ_BLAKE_IMPL */

    /* Make addition in place */
    fmpz_set(P2->x, P4->x);
    fmpz_set(P2->y, P4->y);
    fmpz_set(P2->z, P4->z);
    P2->is_zero = P4->is_zero;
}

void
EllipticCurve_projective_negation(ECProjPoint * P2, const ECProjPoint * P1,
                                  const EllipticCurve * E)
{
    if (P1->is_zero)
    {
        P2->is_zero = 1;
        return;
    }

    fmpz_set(P2->x, P1->x);
    fmpz_neg(P2->y, P1->y);
    fmpz_set(P2->z, P1->z);
}

void
EllipticCurve_projective_to_affine(ECPoint * R, const ECProjPoint * P,
                                   const EllipticCurve * E)
{
    R->is_zero = P->is_zero;

    if (fmpz_equal_ui(P->z, 1))
    {
        fmpz_set(R->x, P->x);
        fmpz_set(R->y, P->y);
        return;
    }

    fmpz_t invZ;
    fmpz_init(invZ);
    fmpz_invmod(invZ, P->z, E->p);
    fmpz_mul(R->x, P->x, invZ);
    fmpz_mod(R->x, R->x, E->p);
    fmpz_mul(R->y, P->y, invZ);
    fmpz_mod(R->y, R->y, E->p);
    fmpz_clear(invZ);
}

void
EllipticCurve_affine_to_projective(ECProjPoint * R, const ECPoint * P,
                                   const EllipticCurve * E)
{
    R->is_zero = P->is_zero;
    fmpz_set(R->x, P->x);
    fmpz_set(R->y, P->y);
    fmpz_set_ui(R->z, 1);
}

void
EllipticCurve_projective_scalarMul(ECProjPoint * P2, const ECProjPoint * P1,
                                   const fmpz_t m, const EllipticCurve * E)
{
    slong i;
    if (fmpz_is_zero(m) || P1->is_zero)
    {
        P2->is_zero = 1;
        return;
    }

    fmpz_t n;
    fmpz_init(n);
    fmpz_abs(n, m);
    ECProjPoint *P3 = new_ECProjPoint_copy(P1);
    ECProjPoint *P4 = new_ECProjPoint();
    fmpz_t *lambda = malloc(ECPROJ_LAMBDA_SIZE * sizeof(fmpz_t));

    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_init(lambda[i]);

    for (i = fmpz_sizeinbase(n, 2) - 2; i >= 0; --i)
    {
        EllipticCurve_projective_double_memory(P3, P3, E, P4, lambda);

        if (fmpz_tstbit(n, i))
            EllipticCurve_projective_add_memory(P3, P3, P1, E, P4, lambda);
    }

    /* Make scalar multiplication in place */
    P2->is_zero = P3->is_zero;
    fmpz_set(P2->x, P3->x);
    if (fmpz_sgn(m) == -1)
        fmpz_sub(P2->y, E->p, P3->y);

    else
        fmpz_set(P2->y, P3->y);

    fmpz_set(P2->z, P3->z);

    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
    {
        fmpz_clear(lambda[i]);
    }
    free(lambda);
    fmpz_clear(n);
    free_ECProjPoint(P4);
    free_ECProjPoint(P3);

}

/*
 * Generates precomputed data for efficient scalar multiplication, when the
 * same point P will be multiplied by a lot of different scalars.
 * After the multiplications, free_precomp_multipleScalarMultiplications
 * should be called.
 */

ECProjPoint ***
precomp_projective_multipleScalarMultiplications(const ECProjPoint * P,
                                                 int word_size, int exp_size,
                                                 const EllipticCurve * curve)
{
    slong i, j;

    /* nbr of chunks */
    int n_words = (exp_size - 1) / word_size + 1;

    /* number of possible non-zero */
    int n_patterns = (1 << word_size) - 1;

    /* words: 0001, 0010, 0011, 0100... */
    ECProjPoint ***Ps = malloc(sizeof(ECProjPoint **) * n_patterns);

    for (j = 0; j < n_patterns; ++j)
        Ps[j] = malloc(sizeof(ECProjPoint *) * (n_words));

    Ps[0][0] = new_ECProjPoint_copy(P);
    if (n_patterns > 1)
    {
        Ps[1][0] = new_ECProjPoint();
        EllipticCurve_projective_double(Ps[1][0], Ps[0][0], curve);
    }
    for (j = 2; j < n_patterns; ++j)
    {
        Ps[j][0] = new_ECProjPoint();
        EllipticCurve_projective_add(Ps[j][0], Ps[j - 1][0], Ps[0][0], curve);
    }

    for (i = 1; i < n_words; ++i)
    {
        Ps[0][i] = new_ECProjPoint();
        EllipticCurve_projective_double(Ps[0][i],
                                        Ps[(1 << (word_size - 1)) - 1][i - 1],
                                        curve);

        if (n_patterns > 1)
        {
            Ps[1][i] = new_ECProjPoint();
            EllipticCurve_projective_double(Ps[1][i], Ps[0][i], curve);
        }

        for (j = 2; j < n_patterns; ++j)
        {
            Ps[j][i] = new_ECProjPoint();
            EllipticCurve_projective_add(Ps[j][i], Ps[j - 1][i], Ps[0][i],
                                         curve);
        }
    }

    return Ps;
}

void
free_precomp_projective_multipleScalarMultiplications(ECProjPoint *** Ps,
                                                      int word_size,
                                                      int exp_size)
{
    slong i, j;

    /* nbr of chunks */
    int n_words = (exp_size - 1) / word_size + 1;

    /* number of possible non-zero */
    int n_patterns = (1 << word_size) - 1;

    /* words: 0001, 0010, 0011, 0100... */
    for (j = 0; j < n_patterns; ++j)
    {
        for (i = 0; i < n_words; ++i)
            free_ECProjPoint(Ps[j][i]);

        free(Ps[j]);
    }
    free(Ps);
}

/*
 * Optimized for multiplication by positive integer, using pre-allocated memory
 */
 
void
projective_scalarMultiplications_withPrecomp_optimized(ECProjPoint * P,
                                                       const fmpz_t m,
                                                       int word_size,
                                                       ECProjPoint *** Ps,
                                                       const EllipticCurve *
                                                       curve, ECProjPoint * P4,
                                                       fmpz_t * lambda)
{
    slong i, j;
    P->is_zero = 1;

    if (fmpz_is_zero(m))
        return;

    int n_chunks = fmpz_sizeinbase(m, 1 << word_size);
    int chunk;

    for (i = 0; i < n_chunks; ++i)
    {
        chunk = 0;

        for (j = 0; j < word_size; ++j)
            chunk |= (fmpz_tstbit(m, (i * word_size) + j) << j);

        if (chunk)
            EllipticCurve_projective_add_memory(P, P, Ps[chunk - 1][i], curve,
                                                P4, lambda);
    }
}

void
projective_scalarMultiplications_withPrecomp(ECProjPoint * P,
                                             const fmpz_t m, int word_size,
                                             ECProjPoint *** Ps,
                                             const EllipticCurve * curve)
{
    slong i;
    fmpz_t n;
    fmpz_init(n);
    fmpz_abs(n, m);

    ECProjPoint *P4 = new_ECProjPoint();
    fmpz_t *lambda = malloc(ECPROJ_LAMBDA_SIZE * sizeof(fmpz_t));
    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_init(lambda[i]);

    projective_scalarMultiplications_withPrecomp_optimized(P, n, word_size, Ps,
                                                           curve, P4, lambda);

    if (fmpz_sgn(m) == -1)
        fmpz_sub(P->y, curve->p, P->y);

    for (i = 0; i < ECPROJ_LAMBDA_SIZE; ++i)
        fmpz_clear(lambda[i]);

    free(lambda);
    free_ECProjPoint(P4);
}
