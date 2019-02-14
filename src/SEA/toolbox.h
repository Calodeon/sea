
#ifndef TOOLBOX_H
#define TOOLBOX_H

#include "fmpz.h"
#include "params.h"
#include "ModularPolynomial.h"

void tic(int i);

void tac(int i);                /* prints time in ms since last tic */

void powerXMod(fmpz_mod_poly_t result, const fmpz_t e,
               const fmpz_mod_poly_t f);

void nextPrime(fmpz_t p, const fmpz_t n);

int nextPrime_ui(int n);

int rootsInFp(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p);

int linearRoots(fmpz ** roots, const fmpz_mod_poly_t poly, const fmpz_t p);


void numeratorRationalComposition(fmpz_mod_poly_t result,
                                  const fmpz_mod_poly_t A,
                                  const fmpz_mod_poly_t B,
                                  const fmpz_mod_poly_t C);

int findNonSquareInFiniteField(int l);

#endif /* TOOLBOX_H */
