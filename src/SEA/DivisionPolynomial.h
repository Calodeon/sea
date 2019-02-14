/* division polynomials */

#ifndef DIVISION_POLYNOMIAL_H
#define DIVISION_POLYNOMIAL_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

#include "EllipticCurve.h"

typedef struct
{
    fmpz_mod_poly_t *f;         /* the polynomials, f_n = psi_n if n odd,
                                   and psi_n/2y if n even */
    const EllipticCurve *curve;
    fmpz_mod_poly_t square_y;
    fmpz_mod_poly_t square_2y;
    fmpz_mod_poly_t quartic_2y;
} DivisionPolynomials;

DivisionPolynomials *new_DivisionPolynomials(const EllipticCurve * curve);

void print_DivisionPolynomials(const DivisionPolynomials * poly);

void free_DivisionPolynomials(DivisionPolynomials * poly);

/*
 * check==true => will check if h is really a factor of l-th division polynomial
 */

int findSignOfTrace(int *eigenvalue, const DivisionPolynomials * divPolys,
                const int l,
                const fmpz_mod_poly_t h,
                const int lambda_guess,
                const int tl_guess);

int findTrace(int *eigenvalue, const DivisionPolynomials * divPolys,
              const int l, const int n, const int lambda_prev,
              const fmpz_mod_poly_t h, int twist_secure);

void isogenyPhiFromPsi(fmpz_mod_poly_t phi, const fmpz_mod_poly_t psi,
                       const EllipticCurve * curve);

#endif /* DIVISION_POLYNOMIAL_H */
