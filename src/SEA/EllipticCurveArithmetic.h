
#ifndef ELLIPTIC_CURVE_ARITHMETIC_H
#define ELLIPTIC_CURVE_ARITHMETIC_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ModularPolynomial.h"
#include "toolbox.h"

/* Affine point on a curve */
typedef struct
{
    fmpz_t x;
    fmpz_t y;
    int is_zero;

} ECPoint;

#include "EllipticCurve.h"

ECPoint *new_ECPoint_copy(const ECPoint * P1);

ECPoint *new_ECPoint();

void free_ECPoint(ECPoint * P);

void print_ECPoint(ECPoint * P);

void EllipticCurve_randomPoint(ECPoint * P, flint_rand_t state,
                               const EllipticCurve * E);

void EllipticCurve_scalarMul(ECPoint * P2, const ECPoint * P1, const fmpz_t m,
                             const EllipticCurve * E);

void EllipticCurve_double(ECPoint * P2, const ECPoint * P1,
                          const EllipticCurve * E);

void EllipticCurve_double_memory(ECPoint * P2, const ECPoint * P1,
                                 const EllipticCurve * E, ECPoint * P4,
                                 fmpz_t tmp, fmpz_t s);

void EllipticCurve_add(ECPoint * P3, const ECPoint * P1, const ECPoint * P2,
                       const EllipticCurve * E);

void EllipticCurve_add_memory(ECPoint * P3, const ECPoint * P1,
                              const ECPoint * P2, const EllipticCurve * E,
                              ECPoint * P4, fmpz_t tmp, fmpz_t s);

ECPoint ***precomp_multipleScalarMultiplications(const ECPoint * P,
                                                 int word_size, int exp_size,
                                                 const EllipticCurve * curve);

void free_precomp_multipleScalarMultiplications(ECPoint *** Ps, int word_size,
                                                int exp_size);

void scalarMultiplications_withPrecomp_optimized(ECPoint * P, const fmpz_t m,
                                                 int word_size, ECPoint *** Ps,
                                                 const EllipticCurve * curve,
                                                 ECPoint * P4, fmpz_t tmp,
                                                 fmpz_t s);

void scalarMultiplications_withPrecomp(ECPoint * P, const fmpz_t m,
                                       int word_size, ECPoint *** Ps,
                                       const EllipticCurve * curve);

#endif /* ELLIPTIC_CURVE_ARITHMETIC_H */
