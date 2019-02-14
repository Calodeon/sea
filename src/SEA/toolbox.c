
#include <stdio.h>
#include <time.h>

#include "toolbox.h"

static clock_t global_timer[100];

void
tic(int i)
{
    global_timer[i] = clock();
}

void
tac(int i)
{
    /* printf("[%.2f ms]\t", (1000.*(float)(clock() - global_timer[i])/CLOCKS_PER_SEC)); */
    printf("[%d ms]\t",
           (int) (1000. * (float) (clock() - global_timer[i]) / CLOCKS_PER_SEC));
}

void
powerXMod(fmpz_mod_poly_t result, const fmpz_t e, const fmpz_mod_poly_t f)
{
    fmpz_mod_poly_t finv;

    fmpz_mod_poly_init(finv, fmpz_mod_poly_modulus(result));

    fmpz_mod_poly_reverse(finv, f, f->length);
    fmpz_mod_poly_inv_series_newton(finv, finv, f->length);

    fmpz_mod_poly_zero(result);
    fmpz_mod_poly_powmod_x_fmpz_preinv(result, e, f, finv);

    fmpz_mod_poly_clear(finv);
}

void
nextPrime(fmpz_t p, const fmpz_t n)
{
    if (fmpz_is_even(n))
        fmpz_add_ui(p, n, 1);
    else
        fmpz_add_ui(p, n, 2);

    while (!fmpz_is_probabprime(p))
        fmpz_add_ui(p, p, 2);
}

/* n >= 2 */
int
nextPrime_ui(int n)
{
    int p;

    if ((n & 1) == 0)
        p = n + 1;
    else
        p = n + 2;

    while (!n_is_oddprime_small(p))
        p += 2;

    return p;
}

int
rootsQuadraticPolynomial(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p)
{
    int d = 0;
    fmpz_t delta,tmp, inv_2a;
    fmpz_init(delta);
    fmpz_init(tmp);
    fmpz_init(inv_2a);

    fmpz_mod_poly_get_coeff_fmpz(delta, poly, 2);
    fmpz_mod_poly_get_coeff_fmpz(tmp, poly, 0);
    fmpz_mul(delta, delta, tmp);
    fmpz_mod(delta, delta, p);
    fmpz_mul_ui(delta, delta, 4);

    fmpz_mod_poly_get_coeff_fmpz(tmp, poly, 1);
    fmpz_mul(tmp, tmp, tmp);
    fmpz_mod(tmp, tmp, p);

    fmpz_sub(delta, tmp, delta);
    fmpz_mod(delta, delta, p);

    if(!fmpz_sqrtmod(delta,delta,p))
    {
        d = 0;
        *roots = NULL;
    }
    else
    {
        d = 2;
        fmpz *find_roots = malloc(sizeof(fmpz) * 2);
        fmpz_init(find_roots);
        fmpz_init(find_roots+1);

        fmpz_mod_poly_get_coeff_fmpz(inv_2a, poly, 2);
        fmpz_mul_ui(inv_2a, inv_2a, 2);
        fmpz_invmod(inv_2a, inv_2a, p);

        fmpz_mod_poly_get_coeff_fmpz(find_roots, poly, 1);
        fmpz_add(find_roots, find_roots, delta);
        fmpz_neg(find_roots,find_roots);
        fmpz_mod(find_roots, find_roots, p);
        fmpz_mul(find_roots, find_roots, inv_2a);
        fmpz_mod(find_roots, find_roots, p);

        fmpz_mod_poly_get_coeff_fmpz(find_roots+1, poly, 1);
        fmpz_sub(find_roots+1, delta, find_roots+1);
        fmpz_mod(find_roots+1, find_roots+1, p);
        fmpz_mul(find_roots+1, find_roots+1, inv_2a);
        fmpz_mod(find_roots+1, find_roots+1, p);

        *roots = find_roots;
    }


    fmpz_clear(delta);
    fmpz_clear(tmp);
    fmpz_clear(inv_2a);

    return d;
}

int
rootLinearPolynomial(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p)
{
    fmpz_t inv_a;

    fmpz_init(inv_a);
    fmpz_mod_poly_get_coeff_fmpz(inv_a, poly, 1);
    fmpz_invmod(inv_a, inv_a, p);

    fmpz *find_roots = malloc(sizeof(fmpz));
    fmpz_init(find_roots);

    fmpz_mod_poly_get_coeff_fmpz(find_roots, poly, 0);
    fmpz_neg(find_roots,find_roots);
    fmpz_mul(find_roots, find_roots, inv_a);
    fmpz_mod(find_roots, find_roots, p);

    *roots = find_roots;
    fmpz_clear(inv_a);
    return 1;
}

int
rootsInFp(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p)
{
    if (fmpz_mod_poly_degree(poly) > 2)
    {
        slong i;
        fmpz_mod_poly_factor_t fac;
        fmpz_mod_poly_t polynomial, x_p;

        fmpz_mod_poly_init(polynomial, p);
        fmpz_mod_poly_init(x_p, p);

        fmpz_mod_poly_set_coeff_ui(polynomial, 1, 1);   /* polynomial = x */

        powerXMod(x_p, p, poly);

        fmpz_mod_poly_sub(polynomial, x_p, polynomial); /* polynomial = x^p - x */
        fmpz_mod_poly_gcd(polynomial, polynomial, poly);


        int d = fmpz_mod_poly_degree(polynomial);

        if (d > 2)
        {
            fmpz *find_roots = malloc(sizeof(fmpz) * d);

            fmpz_mod_poly_factor_init(fac);
            fmpz_mod_poly_factor(fac, polynomial);
            for (i = 0; i < d; ++i)
            {
                fmpz_init(find_roots + i);
                fmpz_mod_poly_get_coeff_fmpz(find_roots + i, fac->poly + i, 0);
                fmpz_sub(find_roots + i, p, find_roots + i);
            }
            *roots = find_roots;
            fmpz_mod_poly_factor_clear(fac);

            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);

            return d;
        }
        else if (d == 2)
        {
            d = rootsQuadraticPolynomial(roots,polynomial,p);
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return d;
        }
        else if (d == 1)
        {
            rootLinearPolynomial(roots,polynomial,p);
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return 1;
        }
        else
        {
            *roots = NULL;
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return 0;
        }
    }
    else if (fmpz_mod_poly_degree(poly) == 2)
    {
        return rootsQuadraticPolynomial(roots,poly,p);
    }
    else if (fmpz_mod_poly_degree(poly) == 1) {
        return rootLinearPolynomial(roots,poly,p);
    }
    else
    {
        *roots = NULL;
        return 0;
    }
}


int
linearRoots(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p)
{
    slong i;
    fmpz_mod_poly_factor_t fac;

    fmpz_mod_poly_factor_init(fac);

    fmpz_mod_poly_factor(fac, poly);

    int d = fmpz_mod_poly_degree(poly);

    if (d != 0)
    {
        fmpz *find_roots = malloc(sizeof(fmpz) * d);

        for (i = 0; i < d; ++i)
        {
            fmpz_init(find_roots + i);
            fmpz_mod_poly_get_coeff_fmpz(find_roots + i, fac->poly + i, 0);
            fmpz_sub(find_roots + i, p, find_roots + i);
        }
        *roots = find_roots;
    }

    fmpz_mod_poly_factor_clear(fac);

    return d;
}

/*
 * Computes numerator of A(B/C), i.e.
 * C^n*A(B/C) = a_n*B^n + a_{n-1}B^{n-1}C^1 + ... + a_0C^n, n = deg(A)
 */
void
numeratorRationalComposition(fmpz_mod_poly_t result, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_poly_t C)
{
    slong i;
    fmpz_t coeff;
    int n = fmpz_mod_poly_degree(A);

    fmpz_mod_poly_t *Ci = malloc(sizeof(fmpz_mod_poly_t) * (n + 1));
    fmpz_mod_poly_init(Ci[0], fmpz_mod_poly_modulus(A));
    fmpz_mod_poly_set_coeff_ui(Ci[0], 0, 1);

    fmpz_init(coeff);

    for (i = 1; i <= n; ++i)
    {
        fmpz_mod_poly_init(Ci[i], fmpz_mod_poly_modulus(A));
        fmpz_mod_poly_mul(Ci[i], Ci[i - 1], C);
    }

    fmpz_mod_poly_t result_tmp, Bi;
    fmpz_mod_poly_init(result_tmp, fmpz_mod_poly_modulus(A));
    fmpz_mod_poly_get_coeff_fmpz(coeff, A, 0);
    fmpz_mod_poly_scalar_mul_fmpz(Ci[n], Ci[n], coeff);
    fmpz_mod_poly_set(result_tmp, Ci[n]);
    fmpz_mod_poly_init(Bi, fmpz_mod_poly_modulus(A));
    fmpz_mod_poly_set(Bi, B);

    for (i = 1; i <= n; ++i)
    {
        fmpz_mod_poly_get_coeff_fmpz(coeff, A, i);
        fmpz_mod_poly_mul(Ci[n - i], Ci[n - i], Bi);
        fmpz_mod_poly_scalar_mul_fmpz(Ci[n - i], Ci[n - i], coeff);
        fmpz_mod_poly_mul(Bi, Bi, B);
        fmpz_mod_poly_add(result_tmp, result_tmp, Ci[n - i]);
    }

    fmpz_mod_poly_set(result, result_tmp);

    for (i = 0; i <= n; ++i)
        fmpz_mod_poly_clear(Ci[i]);

    free(Ci);
    fmpz_mod_poly_clear(result_tmp);
    fmpz_mod_poly_clear(Bi);
    fmpz_clear(coeff);
}

int
findNonSquareInFiniteField(int l)
{
    slong i;
    for (i = 1; i < l; ++i)
        if (n_jacobi(i, l) == -1)
            return i;

    return 0;
}
