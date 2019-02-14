#include <stdio.h>
#include <string.h>
#include "EllipticCurvePointCounting.h"
#include "ModularPolynomial.h"
#include "DivisionPolynomial.h"
#include "toolbox.h"

/*
 * Compute isogenous curve and kernel with Atkin modular polynomials
 * Algo 27 from "Mathematics of public key cryptography" (Galbraith) pp. 527
 */

int
computeIsogeny_AtkinPolynomials(fmpz_mod_poly_t D, fmpz_t A2_, fmpz_t B2_,
                                const EllipticCurve * E1, const fmpz_t j1,
                                const fmpz_t j2, const fmpz_t f1,
                                const ulong l, const fmpz_mod_poly_t psi_j1,
                                const fmpz_mod_poly_t psi_f1,
                                const ModularPolynomial * psi,
                                int check_result)
{
    fmpz_t m, j1_prime, j2_prime, s, p1, k, m2, k2, tmp, delta, rg, f1_prime;
    fmpz_t A2, B2;
    fmpz_t psiJ1, psiF1, psiJ2, psiF2, psiFF1, psiFF2, psiJJ1, psiJJ2,
        psiFJ1, psiFJ2;
    ulong i, n, j, d = (l - 1) / 2;

    fmpz_mod_poly_t derivative;

    fmpz_init(A2);
    fmpz_init(B2);
    fmpz_init(m);
    fmpz_init(tmp);
    fmpz_init(rg);
    fmpz_init(m2);
    fmpz_init_set_ui(k2, 1728);
    fmpz_init_set_ui(k, 1728);
    fmpz_init(j1_prime);
    fmpz_init(j2_prime);
    fmpz_init(f1_prime);
    fmpz_init(p1);
    fmpz_init(s);
    fmpz_init(delta);
    fmpz_init(psiJ1);
    fmpz_init(psiJ2);
    fmpz_init(psiF1);
    fmpz_init(psiF2);
    fmpz_init(psiFF1);
    fmpz_init(psiFF2);
    fmpz_init(psiJJ1);
    fmpz_init(psiJJ2);
    fmpz_init(psiFJ1);
    fmpz_init(psiFJ2);
    fmpz_mod_poly_init(derivative, E1->p);

    fmpz_mod_poly_derivative(derivative, psi_j1);
    fmpz_mod_poly_evaluate_fmpz(psiF1, derivative, f1);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiFF1, derivative, f1);

    fmpz_mod_poly_derivative(derivative, psi_f1);
    fmpz_mod_poly_evaluate_fmpz(psiJ1, derivative, j1);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiJJ1, derivative, j1);

    fmpz_mod_poly_derivative(derivative, psi_f1);
    fmpz_mod_poly_evaluate_fmpz(psiJ2, derivative, j2);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiJJ2, derivative, j2);

    evaluate_ModularPolynomialAtY(derivative, psi, j2);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiF2, derivative, f1);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiFF2, derivative, f1);

    evaluate_ModularPolynomialDXAtX(derivative, psi, f1);
    fmpz_mod_poly_derivative(derivative, derivative);
    fmpz_mod_poly_evaluate_fmpz(psiFJ1, derivative, j1);
    fmpz_mod_poly_evaluate_fmpz(psiFJ2, derivative, j2);

    /* fmpz_mod_poly_evaluate_fmpz(phiX12, phiX1, f2); */

    /* m = 18 B1/A1 mod p */
    fmpz_invmod(m, E1->A, E1->p);
    fmpz_mul_ui(m, m, 18);
    fmpz_mul(m, m, E1->B);
    fmpz_mod(m, m, E1->p);

    /* j' = mj */
    fmpz_mul(j1_prime, m, j1);
    fmpz_mod(j1_prime, j1_prime, E1->p);

    /* k = j' / (1728 - j) */
    fmpz_sub(k, k, j1);
    fmpz_invmod(k, k, E1->p);
    fmpz_mul(k, k, j1_prime);
    fmpz_mod(k, k, E1->p);

    /* f1' = - j' * psiJ1 / psiF1 */
    fmpz_invmod(f1_prime, psiF1, E1->p);
    fmpz_mul(f1_prime, f1_prime, psiJ1);
    fmpz_neg(f1_prime, f1_prime);
    fmpz_mod(f1_prime, f1_prime, E1->p);
    fmpz_mul(f1_prime, f1_prime, j1_prime);
    fmpz_mod(f1_prime, f1_prime, E1->p);

    /* j2' = - f1' * psiF2 / psiJ2 */
    fmpz_mul_ui(j2_prime, psiJ2, l);
    fmpz_invmod(j2_prime, j2_prime, E1->p);
    fmpz_mul(j2_prime, j2_prime, psiF2);
    fmpz_neg(j2_prime, j2_prime);
    fmpz_mod(j2_prime, j2_prime, E1->p);
    fmpz_mul(j2_prime, j2_prime, f1_prime);
    fmpz_mod(j2_prime, j2_prime, E1->p);

    /* m2 = j2' / j2 */
    fmpz_invmod(m2, j2, E1->p);
    fmpz_mul(m2, m2, j2_prime);
    fmpz_mod(m2, m2, E1->p);

    /* k2 = j2' / (1728 - j2) */
    fmpz_sub(k2, k2, j2);
    fmpz_invmod(k2, k2, E1->p);
    fmpz_mul(k2, k2, j2_prime);
    fmpz_mod(k2, k2, E1->p);

    /* A2 = l^4 * m2 * k2 / 48 */
    fmpz_set_ui(A2, l * l * l * l);
    fmpz_mul(A2, A2, m2);
    fmpz_mul(A2, A2, k2);
    fmpz_mod(A2, A2, E1->p);
    fmpz_set_ui(tmp, 48);
    fmpz_invmod(tmp, tmp, E1->p);
    fmpz_mul(A2, A2, tmp);
    fmpz_mod(A2, A2, E1->p);

    /* B2 = l^6 * m2^2 * k2 / 864 */
    fmpz_set_ui(B2, l);
    fmpz_pow_ui(B2, B2, 6);
    fmpz_mul(B2, B2, m2);
    fmpz_mul(B2, B2, m2);
    fmpz_mod(A2, A2, E1->p);
    fmpz_mul(B2, B2, k2);
    fmpz_mod(A2, A2, E1->p);
    fmpz_set_ui(tmp, 864);
    fmpz_invmod(tmp, tmp, E1->p);
    fmpz_mul(B2, B2, tmp);
    fmpz_mod(B2, B2, E1->p);

    /* compute rg */
    fmpz_mul(rg, psiFF2, f1_prime);
    fmpz_mod(rg, rg, E1->p);
    fmpz_mul(rg, rg, f1_prime);
    fmpz_mod(rg, rg, E1->p);

    fmpz_mul(tmp, psiJJ2, j2_prime);
    fmpz_mod(tmp, tmp, E1->p);
    fmpz_mul(tmp, tmp, j2_prime);
    fmpz_mod(tmp, tmp, E1->p);
    fmpz_mul_ui(tmp, tmp, l * l);

    fmpz_add(rg, rg, tmp);

    fmpz_mul(tmp, psiFJ2, f1_prime);
    fmpz_mod(tmp, tmp, E1->p);
    fmpz_mul(tmp, tmp, j2_prime);
    fmpz_mod(tmp, tmp, E1->p);
    fmpz_mul_ui(tmp, tmp, 2 * l);

    fmpz_add(rg, rg, tmp);

    fmpz_mul(tmp, f1_prime, psiF2);
    fmpz_invmod(tmp, tmp, E1->p);

    fmpz_mul(rg, rg, tmp);
    fmpz_neg(rg, rg);
    fmpz_mod(rg, rg, E1->p);

    /* done computing rg, now compute r = rg - delta */
    fmpz_invmod(delta, j1_prime, E1->p);
    fmpz_mul(delta, delta, psiFF1);
    fmpz_mod(delta, delta, E1->p);
    fmpz_mul(delta, delta, f1_prime);
    fmpz_mod(delta, delta, E1->p);
    fmpz_mul(delta, delta, f1_prime);
    fmpz_mod(delta, delta, E1->p);

    fmpz_mul(tmp, j1_prime, psiJJ1);
    fmpz_mod(tmp, tmp, E1->p);

    fmpz_add(delta, delta, tmp);

    fmpz_mul(tmp, f1_prime, psiFJ1);
    fmpz_mul_ui(tmp, tmp, 2);
    fmpz_mod(tmp, tmp, E1->p);

    fmpz_add(delta, delta, tmp);

    fmpz_invmod(tmp, psiJ1, E1->p);
    fmpz_mul(delta, delta, tmp);
    fmpz_mod(delta, delta, E1->p);

    fmpz_sub(rg, rg, delta);
    fmpz_mod(rg, rg, E1->p);

    /* done computing r, now compute p1 */

    fmpz_mul_ui(p1, rg, 6);

    fmpz_mul_ui(tmp, k2, l);
    fmpz_sub(tmp, tmp, k);
    fmpz_mul_ui(tmp, tmp, 3);

    fmpz_sub(p1, p1, tmp);

    fmpz_mul_ui(tmp, m2, l);
    fmpz_sub(tmp, tmp, m);
    fmpz_mul_ui(tmp, tmp, 4);

    fmpz_add(p1, p1, tmp);

    fmpz_set_ui(tmp, 12);
    fmpz_invmod(tmp, tmp, E1->p);
    fmpz_mul_ui(p1, p1, l);
    fmpz_mul(p1, p1, tmp);

    fmpz_mod(p1, p1, E1->p);

    int extraSpaceForCheck = check_result ? 3 : 0;

    fmpz t[d + 1 + extraSpaceForCheck], c[d + 1 + extraSpaceForCheck];

    /* t(0) */
    fmpz_init_set_ui(t, d);
    /* c(0) */
    fmpz_init(c);

    for (i = 1; i <= d + extraSpaceForCheck; ++i)
    {
        fmpz_init(t + i);
        fmpz_init(c + i);
    }

    /* t(1) */
    fmpz_set_ui(tmp, 2);
    fmpz_invmod(tmp, tmp, E1->p);
    fmpz_mul(t + 1, p1, tmp);
    fmpz_mod(t + 1, t + 1, E1->p);

    if (d > 1)
    {
        /* t(2) */
        fmpz_set_ui(tmp, 30);
        fmpz_invmod(tmp, tmp, E1->p);
        fmpz_mul_ui(t + 2, E1->A, 10 * d - 1);
        fmpz_add(t + 2, t + 2, A2);
        fmpz_neg(t + 2, t + 2);
        fmpz_mul(t + 2, t + 2, tmp);
        fmpz_mod(t + 2, t + 2, E1->p);

        /* c(1) */
        fmpz_mul(c + 1, t, E1->A);
        fmpz_mul_ui(c + 1, c + 1, 2);
        fmpz_mul_ui(tmp, t + 2, 6);
        fmpz_add(c + 1, c + 1, tmp);
        fmpz_mod(c + 1, c + 1, E1->p);

        if (d > 2)
        {
            /* t(3) */
            fmpz_mul(t + 3, t + 1, E1->A);
            fmpz_mod(t + 3, t + 3, E1->p);
            fmpz_mul_ui(t + 3, t + 3, 42);

            fmpz_mul_ui(tmp, E1->B, 28 * d - 1);
            fmpz_add(t + 3, t + 3, tmp);

            fmpz_add(t + 3, B2, t + 3);
            fmpz_neg(t + 3, t + 3);

            fmpz_set_ui(tmp, 70);
            fmpz_invmod(tmp, tmp, E1->p);
            fmpz_mul(t + 3, t + 3, tmp);
            fmpz_mod(t + 3, t + 3, E1->p);

            /* c(2) */
            fmpz_mul(c + 2, E1->A, t + 1);
            fmpz_mul_ui(c + 2, c + 2, 6);
            fmpz_mul(tmp, E1->B, t);
            fmpz_mul_ui(tmp, tmp, 4);
            fmpz_add(c + 2, c + 2, tmp);
            fmpz_mul_ui(tmp, t + 3, 10);
            fmpz_add(c + 2, c + 2, tmp);
            fmpz_mod(c + 2, c + 2, E1->p);
        }

    }


    for (n = 2; n < d + extraSpaceForCheck; ++n)
    {
        fmpz_zero(s);
        for (i = 1; i < n; ++i)
        {
            fmpz_mul(tmp, c + i, c + n - i);
            fmpz_add(s, s, tmp);
            fmpz_mod(s, s, E1->p);
        }
        fmpz_mul_ui(s, s, 3);

        fmpz_mul(c + n + 1, E1->A, c + n - 1);
        fmpz_mod(c + n + 1, c + n + 1, E1->p);
        fmpz_mul_ui(c + n + 1, c + n + 1, (2 * n - 1) * (n - 1));
        fmpz_mul(tmp, E1->B, c + n - 2);
        fmpz_mul_ui(tmp, tmp, (2 * n - 2) * (n - 2));
        fmpz_mod(tmp, tmp, E1->p);
        fmpz_add(c + n + 1, c + n + 1, tmp);
        fmpz_sub(c + n + 1, s, c + n + 1);
        fmpz_set_ui(tmp, (n - 1) * (2 * n + 5));
        fmpz_invmod(tmp, tmp, E1->p);
        fmpz_mul(c + n + 1, c + n + 1, tmp);
        fmpz_mod(c + n + 1, c + n + 1, E1->p);
    }

    for (n = 3; n < d + extraSpaceForCheck; ++n)
    {
        fmpz_mul(t + n + 1, E1->A, t + n - 1);
        fmpz_mod(t + n + 1, t + n + 1, E1->p);
        fmpz_mul_ui(t + n + 1, t + n + 1, (4 * n - 2));
        fmpz_mul(tmp, E1->B, t + n - 2);
        fmpz_mul_ui(tmp, tmp, (4 * n - 4));
        fmpz_mod(tmp, tmp, E1->p);
        fmpz_add(t + n + 1, t + n + 1, tmp);
        fmpz_sub(t + n + 1, c + n, t + n + 1);
        fmpz_set_ui(tmp, 4 * n + 2);
        fmpz_invmod(tmp, tmp, E1->p);
        fmpz_mul(t + n + 1, t + n + 1, tmp);
        fmpz_mod(t + n + 1, t + n + 1, E1->p);
    }

    fmpz_mod_poly_zero(D);
    fmpz_one(c);

    for (n = 1; n <= d; ++n)
    {
        fmpz_zero(s);
        for (i = 1; i <= n; ++i)
        {
            fmpz_mul(tmp, t + i, c + n - i);

            if (i % 2 == 0)
                fmpz_neg(tmp, tmp);

            fmpz_add(s, s, tmp);
            fmpz_mod(s, s, E1->p);
        }
        fmpz_set_ui(tmp, n);
        fmpz_invmod(tmp, tmp, E1->p);
        fmpz_mul(c + n, s, tmp);
        fmpz_mod(c + n, c + n, E1->p);
        /* mod ??? */
    }

    for (i = 0; i <= d; ++i)
    {
        if (i % 2 == 1)
            fmpz_sub(tmp, E1->p, c + i);
        else
            fmpz_set(tmp, c + i);

        fmpz_mod_poly_set_coeff_fmpz(D, d - i, tmp);
    }

    int feed_back = 1;          /* no problem to report */

    if (check_result)
    {
        /* TODO : check whether i = 0 gives any info or never detects any error... */
        for (i = 0; i < extraSpaceForCheck + 1; ++i)
        {
            fmpz_set(s, t + d + i);
            for (j = 0; j < d; ++j)
            {
                fmpz_mod_poly_get_coeff_fmpz(tmp, D, j);
                fmpz_mul(tmp, tmp, t + j + i);
                fmpz_mod(tmp, tmp, E1->p);
                fmpz_add(s, s, tmp);
                fmpz_mod(s, s, E1->p);
            }
            feed_back = fmpz_is_zero(s);
        }
    }

    fmpz_set(A2_, A2);
    fmpz_set(B2_, B2);

    for (i = 0; i <= d + extraSpaceForCheck; ++i)
    {
        fmpz_clear(t + i);
        fmpz_clear(c + i);
    }

    fmpz_mod_poly_clear(derivative);

    fmpz_clear(A2);
    fmpz_clear(B2);
    fmpz_clear(tmp);
    fmpz_clear(rg);
    fmpz_clear(m);
    fmpz_clear(m2);
    fmpz_clear(k2);
    fmpz_clear(k);
    fmpz_clear(j1_prime);
    fmpz_clear(j2_prime);
    fmpz_clear(f1_prime);
    fmpz_clear(p1);
    fmpz_clear(s);
    fmpz_clear(delta);
    fmpz_clear(psiJ1);
    fmpz_clear(psiJ2);
    fmpz_clear(psiF1);
    fmpz_clear(psiF2);
    fmpz_clear(psiFF1);
    fmpz_clear(psiFF2);
    fmpz_clear(psiJJ1);
    fmpz_clear(psiJJ2);
    fmpz_clear(psiFJ1);
    fmpz_clear(psiFJ2);

    return feed_back;           /* no problem ! */
}
