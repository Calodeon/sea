/* reading modular polynomials */

#ifndef MODULAR_POLYNOMIAL_H
#define MODULAR_POLYNOMIAL_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "params.h"

/*
 * (X,Y) corresponds to (F1,F2) (symmetric, like Weber) or (F,j) (asymmetric,
 * like Atkin)
 */

typedef struct
{
    int l;                  /* the prime */
    fmpz_poly_t *coef;      /* coefficient (in terms of a polynomial in X) of
                               Y^n for n = 0, ..., maxDeg */
    unsigned int maxDeg;
} ModularPolynomial;

ModularPolynomial *read_ModularPolynomialSymmetric(const int l,
                                                   char path_model[]);

ModularPolynomial *read_ModularPolynomialAsymmetric(const int l,
                                                    char path_model[]);

/* returns Phi(f,Y) */
void evaluate_ModularPolynomialAtX(fmpz_mod_poly_t result,
                                   const ModularPolynomial * poly,
                                   const fmpz_t f);

/* returns Phi(X,f) */
void evaluate_ModularPolynomialAtY(fmpz_mod_poly_t result,
                                   const ModularPolynomial * poly,
                                   const fmpz_t f);

/* returns d/dX[Phi(X,Y)](f,Y) */
void evaluate_ModularPolynomialDXAtX(fmpz_mod_poly_t result,
                                     const ModularPolynomial * poly,
                                     const fmpz_t f);

void weberToJDerivative(fmpz_t result, const fmpz_t f, const fmpz_t p);

void weberToJDerivativeDouble(fmpz_t FP, fmpz_t FPP, const fmpz_t f,
                              const fmpz_t p);

void print_ModularPolynomial(const ModularPolynomial * poly);

void free_ModularPolynomial(ModularPolynomial * poly);

#endif /* MODULAR_POLYNOMIAL_H */
