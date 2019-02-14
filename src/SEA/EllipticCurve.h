#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include "fmpz.h"

#define EC_RANDOM_POINTS_TO_TEST_SUPERSINGULARITY 10

/* Elliptic curve */
typedef struct
{
    fmpz_t p;                   /* the order of the prime field */
    /* y^2 = x^3 + Ax + B */
    fmpz_t A;
    fmpz_t B;

} EllipticCurve;

EllipticCurve *new_EllipticCurve(const fmpz_t p, const fmpz_t A,
                                 const fmpz_t B);

void free_EllipticCurve(EllipticCurve * curve);

void print_EllipticCurve(const EllipticCurve * curve);

/* returns 0 if y^2 = x^3 + Ax + B is singular, non-zero otherwise */
int jInvariant(fmpz_t j, const EllipticCurve * curve);

int EllipticCurve_isSupersingular_proba(const EllipticCurve * curve,
                                  flint_rand_t state);

#endif /* ELLIPTIC_CURVE_H */
