
#include "EllipticCurveArithmetic.h"
#include "ModularPolynomial.h"
#include "DivisionPolynomial.h"
#include "toolbox.h"
#include <stdio.h>
#include <string.h>


/* Creates a new point */

ECPoint *
new_ECPoint()
{
    ECPoint *P = malloc(sizeof(ECPoint));
    P->is_zero = 1;
    fmpz_init(P->x);
    fmpz_init(P->y);
    return P;
}

/* Creates a copy of a point */

ECPoint *
new_ECPoint_copy(const ECPoint * P1)
{
    ECPoint *P = malloc(sizeof(ECPoint));
    P->is_zero = P1->is_zero;
    if (!P->is_zero)
    {
        fmpz_init_set(P->x, P1->x);
        fmpz_init_set(P->y, P1->y);
    }
    return P;
}

/* Frees the memory */

void
free_ECPoint(ECPoint * P)
{
    fmpz_clear(P->x);
    fmpz_clear(P->y);
    free(P);
}

/* Prints a point "(x, y)" */

void
print_ECPoint(ECPoint * P)
{
    if (P->is_zero)
    {
        printf("O");
    }
    else
    {
        printf("(");
        fmpz_print(P->x);
        printf(", ");
        fmpz_print(P->y);
        printf(")");
    }
}

void
EllipticCurve_randomPoint(ECPoint * P, 
                          flint_rand_t state, const EllipticCurve * E)
{
    fmpz_t y_sqr;
    fmpz_init(y_sqr);

    P->is_zero = 0;

    while (1)
    {
        fmpz_randm(P->x, state, E->p);

        fmpz_mul(y_sqr, P->x, P->x);
        fmpz_mod(y_sqr, y_sqr, E->p);
        fmpz_add(y_sqr, y_sqr, E->A);
        fmpz_mul(y_sqr, y_sqr, P->x);
        fmpz_mod(y_sqr, y_sqr, E->p);
        fmpz_add(y_sqr, y_sqr, E->B);
        fmpz_mod(y_sqr, y_sqr, E->p);

        if (fmpz_jacobi(y_sqr, E->p) == 1)
        {
            fmpz_sqrtmod(P->y, y_sqr, E->p);
            fmpz_clear(y_sqr);
            return;
        }
    }
}

/*
 * Generates precomputed data for efficient scalar multiplication, when the
 * same point P will be multiplied by a lot of different scalars.
 * After the multiplications, free_precomp_multipleScalarMultiplications
 * should be called.
 */

ECPoint ***
precomp_multipleScalarMultiplications(const ECPoint * P,
                                      int word_size, int exp_size,
                                      const EllipticCurve * curve)
{
    slong i, j;

    /* nbr of chunks */
    int n_words = (exp_size - 1) / word_size + 1;

    /* number of possible non-zero */
    int n_patterns = (1 << word_size) - 1;

    /* words: 0001, 0010, 0011, 0100... */
    ECPoint ***Ps = malloc(sizeof(ECPoint **) * n_patterns);

    for (j = 0; j < n_patterns; ++j)
        Ps[j] = malloc(sizeof(ECPoint *) * (n_words));

    Ps[0][0] = new_ECPoint_copy(P);
    if (n_patterns > 1)
    {
        Ps[1][0] = new_ECPoint();
        EllipticCurve_double(Ps[1][0], Ps[0][0], curve);
    }
    for (j = 2; j < n_patterns; ++j)
    {
        Ps[j][0] = new_ECPoint();
        EllipticCurve_add(Ps[j][0], Ps[j - 1][0], Ps[0][0], curve);
    }

    for (i = 1; i < n_words; ++i)
    {
        Ps[0][i] = new_ECPoint();
        EllipticCurve_double(Ps[0][i], Ps[(1 << (word_size - 1)) - 1][i - 1],
                             curve);

        if (n_patterns > 1)
        {
            Ps[1][i] = new_ECPoint();
            EllipticCurve_double(Ps[1][i], Ps[0][i], curve);
        }

        for (j = 2; j < n_patterns; ++j)
        {
            Ps[j][i] = new_ECPoint();
            EllipticCurve_add(Ps[j][i], Ps[j - 1][i], Ps[0][i], curve);
        }
    }

    return Ps;
}

void
free_precomp_multipleScalarMultiplications(ECPoint *** Ps, int word_size,
                                           int exp_size)
{
    slong j, i;

    /* nbr of chunks */
    int n_words = (exp_size - 1) / word_size + 1;

    /* number of possible non-zero */
    int n_patterns = (1 << word_size) - 1;

    /* words: 0001, 0010, 0011, 0100... */
    for (j = 0; j < n_patterns; ++j)
    {
        for (i = 0; i < n_words; ++i)
        {
            free_ECPoint(Ps[j][i]);
        }
        free(Ps[j]);
    }
    free(Ps);
}

/*
 * Optimized for multiplication by positive integer, using pre-allocated memory
 */

void
scalarMultiplications_withPrecomp_optimized(ECPoint * P, const fmpz_t m,
                                            int word_size, ECPoint *** Ps,
                                            const EllipticCurve * curve,
                                            ECPoint * P4, fmpz_t tmp, fmpz_t s)
{
    slong i, j;
    P->is_zero = 1;

    if (fmpz_is_zero(m))
        return;
    else {

        int n_chunks = fmpz_sizeinbase(m, 1 << word_size), chunk;

        for (i = 0; i < n_chunks; ++i)
        {
            chunk = 0;

            for (j = 0; j < word_size; ++j)
                chunk |= (fmpz_tstbit(m, (i * word_size) + j) << j);

            if (chunk)
                EllipticCurve_add_memory(P, P, Ps[chunk - 1][i], curve, P4,
                                         tmp, s);
        }
    }
}

void
scalarMultiplications_withPrecomp(ECPoint * P, const fmpz_t m,
                                  int word_size, ECPoint *** Ps,
                                  const EllipticCurve * curve)
{
    fmpz_t n;

    fmpz_t tmp, s;
    ECPoint *P4 = new_ECPoint();

    fmpz_init(n);
    fmpz_abs(n, m);
    fmpz_init(tmp);
    fmpz_init(s);

    scalarMultiplications_withPrecomp_optimized(P, n, word_size, Ps, curve, P4,
                                                tmp, s);

    if (fmpz_sgn(m) == -1)
        fmpz_sub(P->y, curve->p, P->y);

    fmpz_clear(n);
    fmpz_clear(tmp);
    fmpz_clear(s);
    free_ECPoint(P4);
}

/* m*(x1,y1) */

void
EllipticCurve_scalarMul(ECPoint * P2, const ECPoint * P1, const fmpz_t m,
                        const EllipticCurve * E)
{

    if (fmpz_is_zero(m) || P1->is_zero)
    {
        P2->is_zero = 1;
        return;
    }
    else {

        slong i;
        fmpz_t tmp, s, n;
        ECPoint *P4 = new_ECPoint();
        ECPoint *P3 = new_ECPoint_copy(P1);

        fmpz_init(tmp);
        fmpz_init(s);
        fmpz_init(n);
        fmpz_abs(n, m);

        for (i = fmpz_sizeinbase(n, 2) - 2; i >= 0; i--)
        {
            EllipticCurve_double_memory(P3, P3, E, P4, tmp, s);

            if (fmpz_tstbit(n, i))
                EllipticCurve_add_memory(P3, P3, P1, E, P4, tmp, s);
        }

        P2->is_zero = P3->is_zero;
        fmpz_set(P2->x, P3->x);

        if (fmpz_sgn(m) == -1)
            fmpz_sub(P2->y, E->p, P3->y);
        else
            fmpz_set(P2->y, P3->y);

        fmpz_clear(n);
        fmpz_clear(tmp);
        fmpz_clear(s);
        free_ECPoint(P4);
        free_ECPoint(P3);
    }
}

/* 2*(x1,y1) */

void
EllipticCurve_double_memory(ECPoint * P2, const ECPoint * P1,
                            const EllipticCurve * E, ECPoint * P4, fmpz_t tmp,
                            fmpz_t s)
{
    if (P1->is_zero || fmpz_is_zero(P1->y))
        P2->is_zero = 1;
    else
    {
        P2->is_zero = 0;

        fmpz_mul(s, P1->x, P1->x);
        fmpz_mod(s, s, E->p);
        fmpz_mul_ui(s, s, 3);
        fmpz_add(s, s, E->A);
        fmpz_mul_ui(tmp, P1->y, 2);
        fmpz_invmod(tmp, tmp, E->p);
        fmpz_mul(s, s, tmp);
        fmpz_mod(s, s, E->p);

        fmpz_mul(P4->x, s, s);
        fmpz_mul_ui(tmp, P1->x, 2);
        fmpz_sub(P4->x, P4->x, tmp);
        fmpz_mod(P4->x, P4->x, E->p);

        fmpz_sub(P4->y, P1->x, P4->x);
        fmpz_mul(P4->y, P4->y, s);
        fmpz_sub(P4->y, P4->y, P1->y);
        fmpz_mod(P4->y, P4->y, E->p);

        fmpz_set(P2->x, P4->x);
        fmpz_set(P2->y, P4->y);
    }
}

/* P3 = P1 + P2 on E */
void
EllipticCurve_add_memory(ECPoint * P3, const ECPoint * P1, const ECPoint * P2,
                         const EllipticCurve * E, ECPoint * P4, fmpz_t tmp,
                         fmpz_t s)
{
    if (P1->is_zero)
    {
        P3->is_zero = P2->is_zero;
        fmpz_set(P3->x, P2->x);
        fmpz_set(P3->y, P2->y);
        return;
    }
    else if (P2->is_zero)
    {
        P3->is_zero = P1->is_zero;
        fmpz_set(P3->x, P1->x);
        fmpz_set(P3->y, P1->y);
        return;
    }

    fmpz_sub(tmp, P1->x, P2->x);

    if (fmpz_invmod(tmp, tmp, E->p))
    {
        P3->is_zero = 0;

        fmpz_sub(s, P1->y, P2->y);
        fmpz_mul(s, s, tmp);
        fmpz_mod(s, s, E->p);

        fmpz_mul(P4->x, s, s);
        fmpz_sub(P4->x, P4->x, P1->x);
        fmpz_sub(P4->x, P4->x, P2->x);
        fmpz_mod(P4->x, P4->x, E->p);

        fmpz_sub(P4->y, P2->x, P4->x);
        fmpz_mul(P4->y, P4->y, s);
        fmpz_sub(P4->y, P4->y, P2->y);

        fmpz_mod(P3->y, P4->y, E->p);
        fmpz_set(P3->x, P4->x);
    }
    else if (fmpz_equal(P1->y, P2->y))
        EllipticCurve_double(P3, P1, E);
    else
        P3->is_zero = 1;
}

/* 2*(x1,y1) */
void
EllipticCurve_double(ECPoint * P2, const ECPoint * P1, const EllipticCurve * E)
{
    if (P1->is_zero || fmpz_is_zero(P1->y))
        P2->is_zero = 1;
    else
    {
        fmpz_t tmp, s;
        ECPoint *P4 = new_ECPoint();

        fmpz_init(tmp);
        fmpz_init(s);

        fmpz_mul(s, P1->x, P1->x);
        fmpz_mod(s, s, E->p);
        fmpz_mul_ui(s, s, 3);
        fmpz_add(s, s, E->A);
        fmpz_mul_ui(tmp, P1->y, 2);
        fmpz_invmod(tmp, tmp, E->p);
        fmpz_mul(s, s, tmp);
        fmpz_mod(s, s, E->p);

        fmpz_mul(P4->x, s, s);
        fmpz_mul_ui(tmp, P1->x, 2);
        fmpz_sub(P4->x, P4->x, tmp);
        fmpz_mod(P4->x, P4->x, E->p);

        fmpz_sub(P4->y, P1->x, P4->x);
        fmpz_mul(P4->y, P4->y, s);
        fmpz_sub(P4->y, P4->y, P1->y);
        fmpz_mod(P4->y, P4->y, E->p);

        fmpz_set(P2->x, P4->x);
        fmpz_set(P2->y, P4->y);

        P2->is_zero = 0;

        fmpz_clear(tmp);
        fmpz_clear(s);
        free_ECPoint(P4);
    }
}

/* P3 = P1 + P2 on E */
void
EllipticCurve_add(ECPoint * P3, const ECPoint * P1, const ECPoint * P2,
                  const EllipticCurve * E)
{
    if (P1->is_zero)
    {
        P3->is_zero = P2->is_zero;
        fmpz_set(P3->x, P2->x);
        fmpz_set(P3->y, P2->y);
        return;
    }
    else if (P2->is_zero)
    {
        P3->is_zero = P1->is_zero;
        fmpz_set(P3->x, P1->x);
        fmpz_set(P3->y, P1->y);
        return;
    }
    else {

        fmpz_t tmp;

        fmpz_init(tmp);
        fmpz_sub(tmp, P1->x, P2->x);

        if (fmpz_invmod(tmp, tmp, E->p))
        {
            fmpz_t s;
            ECPoint *P4 = new_ECPoint();

            fmpz_init(s);

            fmpz_sub(s, P1->y, P2->y);
            fmpz_mul(s, s, tmp);
            fmpz_mod(s, s, E->p);

            fmpz_mul(P4->x, s, s);
            fmpz_mod(P4->x, P4->x, E->p);
            fmpz_sub(P4->x, P4->x, P1->x);
            fmpz_sub(P4->x, P4->x, P2->x);
            fmpz_mod(P4->x, P4->x, E->p);

            fmpz_sub(P4->y, P2->x, P4->x);
            fmpz_mul(P4->y, P4->y, s);
            fmpz_sub(P4->y, P4->y, P2->y);
            fmpz_mod(P4->y, P4->y, E->p);

            fmpz_set(P3->x, P4->x);
            fmpz_set(P3->y, P4->y);

            P3->is_zero = 0;

            fmpz_clear(tmp);
            fmpz_clear(s);
            free_ECPoint(P4);
        }
        else if (fmpz_equal(P1->y, P2->y))
        {
            fmpz_clear(tmp);
            EllipticCurve_double(P3, P1, E);
        }
        else
        {
            P3->is_zero = 1;
            fmpz_clear(tmp);
        }
    }
}
