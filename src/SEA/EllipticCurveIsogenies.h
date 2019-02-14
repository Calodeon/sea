
#ifndef ELLIPTIC_CURVE_ISOGENIES_H
#define ELLIPTIC_CURVE_ISOGENIES_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ModularPolynomial.h"
#include "EllipticCurve.h"

/*
 * Compute isogenous curve and kernel with Atkin modular polynomials
 * Algo 27 from "Mathematics of public key cryptography" (Galbraith) pp. 527
 *
 * psi_j1 = psi(F,j1), psi_f1 = psi(f1,J), psi is the Atkin modular polynomial
 */

int computeIsogeny_AtkinPolynomials(fmpz_mod_poly_t D, fmpz_t A2_, fmpz_t B2_,
                                    const EllipticCurve * E1, const fmpz_t j1,
                                    const fmpz_t j2, const fmpz_t f1,
                                    const ulong l,
                                    const fmpz_mod_poly_t psi_j1,
                                    const fmpz_mod_poly_t psi_f1,
                                    const ModularPolynomial * psi,
                                    int check_result);

#endif /* ELLIPTIC_CURVE_ISOGENIES_H */
