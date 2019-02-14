#ifndef ELLIPTIC_CURVE_POINT_COUNTING_H
#define ELLIPTIC_CURVE_POINT_COUNTING_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

#include "ModularPolynomial.h"

#include "EllipticCurve.h"
#include "EllipticCurveArithmetic.h"

/***********************************
 *  SEA: Counting rational points  *
 ***********************************/

/*
 * Classification of abort strategies
 *
 *   NONE: compute the number of points, no early abort
 *   PRIME: abort if the curve has an order which is not prime
 *   MONTGOMERY_4: abort if the curve is not Montgomery, or if its order is 
 *                 not 4*q for a prime q
 *   MONTGOMERY_8: abort if the curve is not Montgomery, or if its order is 
 *                 not 4*q or 8*q for a prime q
 */

typedef enum
{ NONE, PRIME, MONTGOMERY_4, MONTGOMERY_8 } AbortStrategy;

/*
 * A note on montgomery + small cofactor + twist-security
 *
 * We have:
 *   E is Montgomery => 4 | #E(F_p)
 *   E is Montgomery => twist of E is Montgomery
 *
 * So for Montgomery curves, both the curve and the twist will have a cofactor
 * at least 4.
 * Furthermore, if p = 1 mod 4, 8 will be a factor of either one of them,
 * and if p = 3 mod 4, 8 is a cofactor of one if and only if it is a cofactor of
 * the other
 */

/*
 * SEA_init should be called once before the first time a point counting is performed (loads
 * the modular polynomials)
 */
void SEA_init();


/*
 * Frees the memory (the modular polynomials)
 */
void SEA_clear();


typedef struct
{
    int l;
    int nonSquare;
    int generator1;
    int generator2;
    ModularPolynomial *modular;
} PrimeData;

PrimeData *SEA_init_smallPrimes(int nbr_small_primes);

int SEA_read_primeData(PrimeData * dest, FILE *file);

void SEA_free_smallPrimes(PrimeData * data, int nbr_small_primes);

/*
 * the SEA_* functions : return -2 in case of error, -1 in case of early abort,
 * 0 if the cardinality has been computed
 */

int SEA_gmp(mpz_t n, const mpz_t p, const mpz_t A, const mpz_t B);

int SEA_simple(fmpz_t n, const fmpz_t p, const fmpz_t A, const fmpz_t B);

int SEA(fmpz_t n, const fmpz_t p, const fmpz_t A, const fmpz_t B,
        AbortStrategy strategy, int twist_secure);

int SEA_AtkinPolynomials(fmpz_t cardinality, const EllipticCurve * curve,
                         PrimeData primes[], flint_rand_t state,
                         AbortStrategy strategy, int twist_secure);

#endif /* ELLIPTIC_CURVE_POINT_COUNTING_H */
