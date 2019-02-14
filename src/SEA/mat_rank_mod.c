#include <stdio.h>
#include "mat_rank_mod.h"

/* similar to fmpz_mat_rref_mod, except the matrix is not reduced */
slong
mat_rank_mod(fmpz_mat_t A, const fmpz_t p)
{
    fmpz_t t, inv;
    slong m, n, j, k, rank, r, pivot_row, pivot_col;

    if (fmpz_mat_is_empty(A))
        return 0;

    m = A->r;
    n = A->c;
    rank = pivot_row = pivot_col = 0;

    fmpz_init(t);
    fmpz_init(inv);

    while (pivot_row < m && pivot_col < n)
    {
        r = fmpz_mat_find_pivot_any(A, pivot_row, m, pivot_col);

        if (r == -1)
        {
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
            fmpz_mat_swap_rows(A, NULL, pivot_row, r);

        rank++;

        fmpz_invmod(inv, E(pivot_row, pivot_col), p);

        /* pivot row */
        for (k = pivot_col + 1; k < n; k++)
        {
            fmpz_mul(E(pivot_row, k), E(pivot_row, k), inv);
            fmpz_mod(E(pivot_row, k), E(pivot_row, k), p);
        }
        fmpz_one(E(pivot_row, pivot_col));

        /* other rows */
        for (j = pivot_row + 1; j < m; j++)
        {
            for (k = pivot_col + 1; k < n; k++)
            {
                fmpz_mul(t, E(j, pivot_col), E(pivot_row, k));
                fmpz_sub(E(j, k), E(j, k), t);
                fmpz_mod(E(j, k), E(j, k), p);
            }
            fmpz_zero(E(j, pivot_col));
        }

        pivot_row++;
        pivot_col++;
    }

    fmpz_clear(inv);
    fmpz_clear(t);

    return rank;
}
