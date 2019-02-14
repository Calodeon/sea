
#include "EllipticCurve.h"
#include "EllipticCurveArithmetic.h"
#include "toolbox.h"
#include <stdio.h>

EllipticCurve *
new_EllipticCurve(const fmpz_t p, const fmpz_t A, const fmpz_t B)
{
    EllipticCurve *curve = malloc(sizeof(EllipticCurve));
    fmpz_init_set(curve->p, p);
    fmpz_init_set(curve->A, A);
    fmpz_init_set(curve->B, B);
    return curve;
}

void
free_EllipticCurve(EllipticCurve * curve)
{
    fmpz_clear(curve->p);
    fmpz_clear(curve->A);
    fmpz_clear(curve->B);
    free(curve);
}

/* Prints the short Weierstrass equation of the curve */
void
print_EllipticCurve(const EllipticCurve * curve)
{
    printf("Y^2 = X^3 + ");
    fmpz_print(curve->A);
    printf("X + ");
    fmpz_print(curve->B);
}

int
jInvariant(fmpz_t j, const EllipticCurve * curve)
{
    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_powm_ui(j, curve->A, 3, curve->p);     /* j = A^3 */
    fmpz_mul_ui(j, j, 4);                       /* j = 4A^3 */
    fmpz_powm_ui(tmp, curve->B, 2, curve->p);   /* tmp = B^2 */
    fmpz_mul_ui(tmp, tmp, 27);                  /* tmp = 27B^2 */
    fmpz_add(tmp, tmp, j);                      /* tmp = 4A^3 + 27B^2 */

    if (fmpz_invmod(tmp, tmp, curve->p) == 0)
    {                           /* tmp = (4A^3 + 27B^2)^-1 mod p */
        fmpz_clear(tmp);
        return 0;
    }

    fmpz_mul(j, j, tmp);        /* j = 4A^3/(4A^3 + 27B^2) */
    fmpz_mul_ui(j, j, 1728);    /* j = 1728 * 4A^3/(4A^3 + 27B^2) */
    fmpz_mod(j, j, curve->p);   /* reduce mod p */

    fmpz_clear(tmp);
    return 1;
}

int
EllipticCurve_isSupersingular_proba(const EllipticCurve * curve, flint_rand_t state)
{
    ECPoint * Q, * R;
    fmpz_t m;
    slong i;

    Q = new_ECPoint();
    R = new_ECPoint();
    fmpz_init(m);

    fmpz_add_ui(m, curve->p, 1UL);

    for (i = 0; i < EC_RANDOM_POINTS_TO_TEST_SUPERSINGULARITY; ++i)
    {
        EllipticCurve_randomPoint(Q, state, curve);
        EllipticCurve_scalarMul(R, Q, m, curve);
        if (!(R->is_zero))
        {
            free_ECPoint(Q);
            free_ECPoint(R);
            fmpz_clear(m);

            return 0;
        }
    }

    free_ECPoint(Q);
    free_ECPoint(R);
    fmpz_clear(m);

    return 1;
}
