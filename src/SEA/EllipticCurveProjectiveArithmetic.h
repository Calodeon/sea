#ifndef ELLIPTIC_CURVE_PROJECTIVE_ARITHMETIC_H
#define ELLIPTIC_CURVE_PROJECTIVE_ARITHMETIC_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

#include "EllipticCurveArithmetic.h"
#include "ModularPolynomial.h"
#include "toolbox.h"

/* Projective point on a curve */
typedef struct
{
    fmpz_t x;
    fmpz_t y;
    fmpz_t z;
    int is_zero;

} ECProjPoint;

#ifdef ECPROJ_BLAKE_IMPL
#define ECPROJ_LAMBDA_SIZE 10
#endif
#ifndef ECPROJ_BLAKE_IMPL
#define ECPROJ_LAMBDA_SIZE 9
#endif


#include "EllipticCurve.h"

ECProjPoint *new_ECProjPoint_copy(const ECProjPoint * P1);
ECProjPoint *new_ECProjPoint();

void free_ECProjPoint(ECProjPoint * P);
void print_ECProjPoint(ECProjPoint * P);

int EllipticCurve_projective_equals(const ECProjPoint * P,
                                    const ECProjPoint * Q,
                                    const EllipticCurve * E);

int EllipticCurve_projective_isEquivalentTo(const ECProjPoint * P,
                                            const ECPoint * a,
                                            const EllipticCurve * E);

int EllipticCurve_projective_isOnCurve(const ECProjPoint * P,
                                       const EllipticCurve * E);

void EllipticCurve_projective_randomPoint(ECProjPoint * P,
                                          flint_rand_t state,
                                          const EllipticCurve * E);

void EllipticCurve_projective_scalarMul(ECProjPoint * P2,
                                        const ECProjPoint * P1, const fmpz_t m,
                                        const EllipticCurve * E);

void EllipticCurve_projective_double(ECProjPoint * P2, const ECProjPoint * P1,
                                     const EllipticCurve * E);

void EllipticCurve_projective_double_memory(ECProjPoint * P2,
                                            const ECProjPoint * P1,
                                            const EllipticCurve * E,
                                            ECProjPoint * P4, fmpz_t * lambda);

void EllipticCurve_projective_add(ECProjPoint * P3, const ECProjPoint * P1,
                                  const ECProjPoint * P2,
                                  const EllipticCurve * E);

void EllipticCurve_projective_add_memory(ECProjPoint * P3,
                                         const ECProjPoint * P1,
                                         const ECProjPoint * P2,
                                         const EllipticCurve * E,
                                         ECProjPoint * P4, fmpz_t * lambda);

ECProjPoint ***precomp_projective_multipleScalarMultiplications(const
                                                                ECProjPoint *
                                                                P,
                                                                int word_size,
                                                                int exp_size,
                                                                const
                                                                EllipticCurve *
                                                                curve);

void free_precomp_projective_multipleScalarMultiplications(ECProjPoint *** Ps,
                                                           int word_size,
                                                           int exp_size);

void projective_scalarMultiplications_withPrecomp_optimized(ECProjPoint * P,
                                                            const fmpz_t m,
                                                            int word_size,
                                                            ECProjPoint *** Ps,
                                                            const EllipticCurve
                                                            * curve,
                                                            ECProjPoint * P4,
                                                            fmpz_t * lambda);

void projective_scalarMultiplications_withPrecomp(ECProjPoint * P,
                                                  const fmpz_t m,
                                                  int word_size,
                                                  ECProjPoint *** Ps,
                                                  const EllipticCurve * curve);

void EllipticCurve_projective_negation(ECProjPoint * P2,
                                       const ECProjPoint * P1,
                                       const EllipticCurve * E);

void EllipticCurve_projective_to_affine(ECPoint * R, const ECProjPoint * P,
                                        const EllipticCurve * E);

void EllipticCurve_affine_to_projective(ECProjPoint * R, const ECPoint * P,
                                        const EllipticCurve * E);

#endif /* ELLIPTIC_CURVE_PROJECTIVE_ARITHMETIC_H */
