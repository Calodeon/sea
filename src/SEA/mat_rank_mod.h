#ifndef MAT_RANK_MOD_H
#define MAT_RANK_MOD_H

#include "fmpz.h"
#include "params.h"
#include "ModularPolynomial.h"
#include "fmpz_mat.h"

#define E(j,k) fmpz_mat_entry(A,j,k)

slong mat_rank_mod(fmpz_mat_t A, const fmpz_t p);

#endif /* MAT_RANK_MOD_H */
