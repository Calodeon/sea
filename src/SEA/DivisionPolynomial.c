/* division polynomials */

#include "DivisionPolynomial.h"
#include "toolbox.h"
#include <stdio.h>
#include <string.h>


typedef struct
{
    fmpz_mod_poly_t *f;         /* The polynomials, f_n = psi_n if n odd,
                                   and psi_n/2y if n even */
    fmpz_mod_poly_t *f_square;  /* The polynomials, f_n = psi_n if n odd,
                                   and psi_n/2y if n even */
    fmpz_mod_poly_t *f_cube;    /* The polynomials, f_n = psi_n if n odd,
                                   and psi_n/2y if n even */
    int *init;                  /* Polynomials computed from 0 to N-1 */
    int Nalloc;                 /* Size of the array psi */
    fmpz_mod_poly_t h;
    fmpz_mod_poly_t hinv;
    const fmpz *p;
    const fmpz_mod_poly_struct *quartic_2y;
} Fs;


DivisionPolynomials *
new_DivisionPolynomials(const EllipticCurve * curve)
{
    fmpz_t tmp;
    DivisionPolynomials *poly = malloc(sizeof(DivisionPolynomials));

    fmpz_init(tmp);
    poly->curve = curve;
    poly->f = malloc(sizeof(fmpz_mod_poly_t) * 5);

    fmpz_mod_poly_init(poly->f[0], curve->p);   /* 0 */

    fmpz_mod_poly_init2(poly->f[1], curve->p, 1);
    fmpz_mod_poly_set_ui(poly->f[1], 1);    /* 1 */

    fmpz_mod_poly_init2(poly->f[2], curve->p, 1);
    fmpz_mod_poly_set_ui(poly->f[2], 1);    /* 1 */

    fmpz_mod_poly_init2(poly->f[3], curve->p, 4);
    fmpz_mul(tmp, curve->A, curve->A);
    fmpz_neg(tmp, tmp);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[3], 0, tmp);
    fmpz_mul_ui(tmp, curve->B, 12);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[3], 1, tmp);
    fmpz_mul_ui(tmp, curve->A, 6);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[3], 2, tmp);
    fmpz_mod_poly_set_coeff_ui(poly->f[3], 4, 3);

    fmpz_mod_poly_init2(poly->f[4], curve->p, 6);
    fmpz_mul(tmp, curve->A, curve->A);
    fmpz_mul(tmp, tmp, curve->A);
    fmpz_mul_ui(tmp, tmp, 2);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[4], 0, tmp);
    fmpz_mul(tmp, curve->B, curve->B);
    fmpz_mul_ui(tmp, tmp, 16);
    fmpz_add(fmpz_poly_lead(poly->f[4]), fmpz_poly_lead(poly->f[4]), tmp);
    fmpz_neg(fmpz_poly_lead(poly->f[4]), fmpz_poly_lead(poly->f[4]));

    fmpz_mul(tmp, curve->A, curve->B);
    fmpz_mul_ui(tmp, tmp, 8);
    fmpz_neg(tmp, tmp);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[4], 1, tmp);

    fmpz_mul(tmp, curve->A, curve->A);
    fmpz_mul_ui(tmp, tmp, 10);
    fmpz_neg(tmp, tmp);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[4], 2, tmp);

    fmpz_mul_ui(tmp, curve->B, 40);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[4], 3, tmp);

    fmpz_mul_ui(tmp, curve->A, 10);
    fmpz_mod_poly_set_coeff_fmpz(poly->f[4], 4, tmp);

    fmpz_mod_poly_set_coeff_ui(poly->f[4], 6, 2);

    fmpz_mod_poly_init2(poly->square_y, curve->p, 3);
    fmpz_mod_poly_init2(poly->square_2y, curve->p, 3);
    fmpz_mod_poly_init2(poly->quartic_2y, curve->p, 6);

    fmpz_mod_poly_set_coeff_fmpz(poly->square_y, 0, curve->B);
    fmpz_mod_poly_set_coeff_fmpz(poly->square_y, 1, curve->A);
    fmpz_mod_poly_set_coeff_ui(poly->square_y, 3, 1);
    fmpz_set_ui(tmp, 4);

    fmpz_mod_poly_scalar_mul_fmpz(poly->square_2y, poly->square_y, tmp);

    fmpz_mod_poly_mul(poly->quartic_2y, poly->square_2y, poly->square_2y);

    fmpz_clear(tmp);

    return poly;
}

void fRecursive(Fs * fs, int i);

void fRecursive_square(Fs * fs, int i);

void fRecursive_cube(Fs * fs, int i);

void
fRecursive_external(fmpz_mod_poly_t f, Fs * fs, int i)
{

    fmpz_mod_poly_t tmp;
    int k = (i >> 1);

    fmpz_mod_poly_init(tmp, fs->p);

    if ((i & 1) == 0)
    {                           /* i even */
        fRecursive_square(fs, k - 1);
        fRecursive(fs, k);
        fRecursive_square(fs, k + 1);
        fRecursive(fs, k + 2);
        fRecursive(fs, k - 2);
        fmpz_mod_poly_mulmod_preinv(f, fs->f[k + 2], fs->f_square[k - 1],
                                        fs->h, fs->hinv);
        fmpz_mod_poly_mulmod_preinv(tmp, fs->f[k - 2], fs->f_square[k + 1],
                                        fs->h, fs->hinv);
        fmpz_mod_poly_sub(f, f, tmp);
        fmpz_mod_poly_mulmod_preinv(f, f, fs->f[k], fs->h, fs->hinv);
    }
    else
    {                           /* i odd */
        fRecursive(fs, k - 1);
        fRecursive_cube(fs, k);
        fRecursive_cube(fs, k + 1);
        fRecursive(fs, k + 2);
        if ((k & 1) == 0)
        {                       /* k even */
            fmpz_mod_poly_mulmod_preinv(f, fs->f[k + 2], fs->f_cube[k],
                                            fs->h, fs->hinv);
            fmpz_mod_poly_mulmod_preinv(f, f, fs->quartic_2y, fs->h,
                                            fs->hinv);
            fmpz_mod_poly_mulmod_preinv(tmp, fs->f[k - 1],
                                            fs->f_cube[k + 1], fs->h,
                                            fs->hinv);
            fmpz_mod_poly_sub(f, f, tmp);
        }
        else
        {                       /* k odd */
            fmpz_mod_poly_mulmod_preinv(f, fs->f[k + 2], fs->f_cube[k],
                                            fs->h, fs->hinv);
            fmpz_mod_poly_mulmod_preinv(tmp, fs->f[k - 1],
                                            fs->f_cube[k + 1], fs->h,
                                            fs->hinv);
            fmpz_mod_poly_mulmod_preinv(tmp, tmp, fs->quartic_2y, fs->h,
                                            fs->hinv);
            fmpz_mod_poly_sub(f, f, tmp);
        }
    }
    fmpz_mod_poly_clear(tmp);
}

void
fRecursive(Fs * fs, int i)
{
    if (fs->init[i])
        return;                 /* already computed */

    fs->init[i] = 1;

    fmpz_mod_poly_init(fs->f[i], fs->p);
    fRecursive_external(fs->f[i], fs, i);
}

void
fRecursive_square(Fs * fs, int i)
{

    if (fs->init[i] > 1)
        return;                 /* already computed */
    if (fs->init[i] == 0)
        fRecursive(fs, i);      /* we'll need this */

    fs->init[i] = 2;

    fmpz_mod_poly_init(fs->f_square[i], fs->p);
    //printf("TAG #1: %ld %ld\n", fmpz_mod_poly_degree(fs->f[i]), fmpz_mod_poly_degree(fs->h));
    fmpz_mod_poly_mulmod_preinv(fs->f_square[i], fs->f[i], fs->f[i], fs->h,
                                    fs->hinv);
    //printf("TAG #2\n");
}

void
fRecursive_cube(Fs * fs, int i)
{
    if (fs->init[i] > 2)
        return;                     /* already computed */
    if (fs->init[i] < 2)
        fRecursive_square(fs, i);   /* we'll need this */

    fs->init[i] = 3;

    fmpz_mod_poly_init(fs->f_cube[i], fs->p);
    fmpz_mod_poly_mulmod_preinv(fs->f_cube[i], fs->f_square[i], fs->f[i],
                                    fs->h, fs->hinv);
}

/*
 * fmpz_mod_poly *get_f(fmpz_mod_poly_t *f, int i, const fmpz_mod_poly_t h,
 *                      const fmpz_mod_poly_t hinv, const DivisionPolynomials
 *                      *divPolys)
 * {
 *     if (f[i]) 
 *         return f[i];
 *     else
 *         return f[i];
 * }
 */

/*
 * Computes the polynomial phi (numerator of x coordinate) from psi, the
 * kernel polynomial
 */
void
isogenyPhiFromPsi(fmpz_mod_poly_t phi, const fmpz_mod_poly_t psi,
                  const EllipticCurve * curve)
{
    fmpz_mod_poly_t phi_tmp, psi_prime, poly, poly2;
    fmpz_t tmp;

    fmpz_mod_poly_init(phi_tmp, curve->p);
    fmpz_mod_poly_init(psi_prime, curve->p);
    fmpz_mod_poly_init(poly, curve->p);
    fmpz_mod_poly_init(poly2, curve->p);

    fmpz_init(tmp);

    fmpz_mod_poly_derivative(psi_prime, psi);

    fmpz_mod_poly_get_coeff_fmpz(tmp, psi, fmpz_mod_poly_degree(psi) - 1);
    fmpz_mul_ui(tmp, tmp, 2);
    fmpz_mod_poly_set_coeff_fmpz(phi_tmp, 0, tmp);
    fmpz_mod_poly_set_coeff_ui(phi_tmp, 1, fmpz_mod_poly_degree(psi) * 2 + 1);
    fmpz_mod_poly_mul(phi_tmp, phi_tmp, psi);
    fmpz_mod_poly_mul(phi_tmp, phi_tmp, psi);

    fmpz_mod_poly_zero(poly);
    fmpz_mod_poly_set_coeff_ui(poly, 2, 6);
    fmpz_mul_ui(tmp, curve->A, 2);
    fmpz_mod_poly_set_coeff_fmpz(poly, 0, tmp);
    fmpz_mod_poly_mul(poly, poly, psi);
    fmpz_mod_poly_mul(poly, poly, psi_prime);

    fmpz_mod_poly_sub(phi_tmp, phi_tmp, poly);

    fmpz_mod_poly_zero(poly);
    fmpz_mod_poly_set_coeff_ui(poly, 3, 1);
    fmpz_mod_poly_set_coeff_fmpz(poly, 1, curve->A);
    fmpz_mod_poly_set_coeff_fmpz(poly, 0, curve->B);
    fmpz_set_ui(tmp, 4);
    fmpz_mod_poly_scalar_mul_fmpz(poly, poly, tmp);

    fmpz_mod_poly_mul(poly2, psi_prime, psi_prime);
    fmpz_mod_poly_derivative(psi_prime, psi_prime);
    fmpz_mod_poly_mul(psi_prime, psi_prime, psi);
    fmpz_mod_poly_sub(poly2, poly2, psi_prime);
    fmpz_mod_poly_mul(poly, poly, poly2);

    fmpz_mod_poly_add(phi_tmp, phi_tmp, poly);
    fmpz_mod_poly_set(phi, phi_tmp);

    fmpz_clear(tmp);
    fmpz_mod_poly_clear(phi_tmp);
    fmpz_mod_poly_clear(psi_prime);
    fmpz_mod_poly_clear(poly);
    fmpz_mod_poly_clear(poly2);
}

Fs *
new_Fs(int alloc, const fmpz_mod_poly_t h, const fmpz_t p,
       const fmpz_mod_poly_t quartic_2y, const DivisionPolynomials * divPolys)
{
    slong i;
    Fs *fs = malloc(sizeof(Fs));

    fs->Nalloc = alloc;

    fs->f = malloc(sizeof(fmpz_mod_poly_t) * alloc);
    fs->f_square = malloc(sizeof(fmpz_mod_poly_t) * alloc);
    fs->f_cube = malloc(sizeof(fmpz_mod_poly_t) * alloc);
    fs->init = malloc(sizeof(int) * alloc);

    fs->p = p;
    fs->quartic_2y = quartic_2y;

    fmpz_mod_poly_init(fs->h, fmpz_mod_poly_modulus(h));
    fmpz_mod_poly_set(fs->h, h);

    for (i = 0; i < 5; ++i)
        fs->init[i] = 3;

    for (i = 5; i < alloc; ++i)
        fs->init[i] = 0;

    fmpz_mod_poly_init(fs->hinv, p);

    fmpz_mod_poly_reverse(fs->hinv, h, h->length);
    fmpz_mod_poly_inv_series_newton(fs->hinv, fs->hinv, h->length);

    fmpz_mod_poly_init(fs->f[0], p);            /* 0 */
    fmpz_mod_poly_init(fs->f_square[0], p);     /* 0 */
    fmpz_mod_poly_init(fs->f_cube[0], p);       /* 0 */

    fmpz_mod_poly_init2(fs->f[1], p, 1);
    fmpz_mod_poly_set_ui(fs->f[1], 1);          /* 1 */
    fmpz_mod_poly_init(fs->f_square[1], p);
    fmpz_mod_poly_set_ui(fs->f_square[1], 1);   /* 1 */
    fmpz_mod_poly_init(fs->f_cube[1], p);
    fmpz_mod_poly_set_ui(fs->f_cube[1], 1);     /* 1 */

    fmpz_mod_poly_init2(fs->f[2], p, 1);
    fmpz_mod_poly_set_ui(fs->f[2], 1);          /* 1 */
    fmpz_mod_poly_init(fs->f_square[2], p);
    fmpz_mod_poly_set_ui(fs->f_square[2], 1);   /* 1 */
    fmpz_mod_poly_init(fs->f_cube[2], p);
    fmpz_mod_poly_set_ui(fs->f_cube[2], 1);     /* 1 */

    fmpz_mod_poly_init(fs->f[3], p);
    fmpz_mod_poly_init(fs->f[4], p);

    fmpz_mod_poly_init(fs->f_square[3], p);
    fmpz_mod_poly_init(fs->f_square[4], p);

    fmpz_mod_poly_init(fs->f_cube[3], p);
    fmpz_mod_poly_init(fs->f_cube[4], p);

    if (fmpz_mod_poly_degree(h) > 4)
    {
        fmpz_mod_poly_set(fs->f[3], divPolys->f[3]);
        if (fmpz_mod_poly_degree(h) > 6)
            fmpz_mod_poly_set(fs->f[4], divPolys->f[4]);
        else
            fmpz_mod_poly_rem(fs->f[4], divPolys->f[4], h);
    }
    else
    {
        fmpz_mod_poly_rem(fs->f[3], divPolys->f[3], h);
        fmpz_mod_poly_rem(fs->f[4], divPolys->f[4], h);
    }

    fmpz_mod_poly_mulmod_preinv(fs->f_square[3], fs->f[3], fs->f[3], h,
                                    fs->hinv);
    fmpz_mod_poly_mulmod_preinv(fs->f_square[4], fs->f[4], fs->f[4], h,
                                    fs->hinv);

    fmpz_mod_poly_mulmod_preinv(fs->f_cube[3], fs->f_square[3], fs->f[3],
                                    h, fs->hinv);
    fmpz_mod_poly_mulmod_preinv(fs->f_cube[4], fs->f_square[4], fs->f[4],
                                    h, fs->hinv);

    return fs;
}

void
free_Fs(Fs * fs)
{
    slong i;

    for (i = 0; i < fs->Nalloc; ++i)
    {
        if (fs->init[i])
            fmpz_mod_poly_clear(fs->f[i]);
        if (fs->init[i] > 1)
            fmpz_mod_poly_clear(fs->f_square[i]);
        if (fs->init[i] > 2)
            fmpz_mod_poly_clear(fs->f_cube[i]);
    }

    fmpz_mod_poly_clear(fs->h);
    fmpz_mod_poly_clear(fs->hinv);

    free(fs->init);
    free(fs->f_square);
    free(fs->f_cube);
    free(fs->f);
    free(fs);
}

/*
 * Checks if t mod l is tl_guess or l-tl_guess
 *
 * Fast method: ASSUMES THAT l % 4 == 3
 */
int
findSignOfTrace_3mod4(int *eigenvalue, const DivisionPolynomials * divPolys,
                const int l,
                const fmpz_mod_poly_t h,
                const fmpz_mod_poly_t square_y,
                const int lambda_guess,
                const int tl_guess)
{
    fmpz_t a, b;
    fmpz_init(a);
    fmpz_init(b);
    int lambda = lambda_guess;
    int tl = tl_guess;
    /* Fast method: no need to compute Y^p */
    int res_pow, jacobi;
    fmpz_mod_poly_resultant_euclidean(b, h, square_y);
    fmpz_sub_ui(a, divPolys->curve->p, 1);
    fmpz_fdiv_q_2exp(a, a, 1);  /* a <- (p-1)/2 */
    fmpz_powm(b, b, a, divPolys->curve->p);

    res_pow = fmpz_is_one(b) ? 1 : -1;

    jacobi = n_powmod(lambda, fmpz_mod_poly_degree(h), l);

    if (jacobi != 1)
        jacobi = -1;

    if (res_pow == jacobi)
        *eigenvalue = lambda;
    else
    {
        tl = l - tl;
        *eigenvalue = l - lambda;
    }

    fmpz_clear(a);
    fmpz_clear(b);
    return tl;
}

/*
 * Checks if t mod l is tl_guess or l-tl_guess
 * To do this, look at Y coordinate
 *
 * Slow method: to be used when l % 4 != 3
 */
int
findSignOfTrace_slow(int *eigenvalue, const DivisionPolynomials * divPolys,
                const int l,
                const fmpz_mod_poly_t h,
                const fmpz_mod_poly_t square_y,
                const fmpz_mod_poly_t square_2y,
                const int lambda_guess,
                const int tl_guess,
                Fs *fs)
{
    fmpz_t a, b;
    fmpz_mod_poly_t tmp;
    fmpz_mod_poly_t tmp2;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_mod_poly_init(tmp, divPolys->curve->p);
    fmpz_mod_poly_init(tmp2, divPolys->curve->p);
    int lambda = lambda_guess;
    int tl = tl_guess;


    fmpz_mod_poly_struct *hinv = fs->hinv;
    fmpz_mod_poly_t *f = fs->f;
    fmpz_mod_poly_t *f_square = fs->f_square;

    fmpz_add_ui(a, divPolys->curve->p, 1);
    fmpz_fdiv_q_2exp(a, a, 1);  /* a <- (p+1)/2 */

    fmpz_mod_poly_set(tmp, square_y);

    /* BIG EXPONENTIATION */
    if (TALKATIVE)
        tic(7);

    /* tmp <- y^{p+1} */
    fmpz_mod_poly_powmod_fmpz_binexp_preinv(tmp, tmp, a, h, hinv);

    if (TALKATIVE >= 3)
    {
        printf("y^p");
        tac(7);
    }

    /* we need f[lambda + 2] */
    fRecursive(fs, lambda - 1);
    fRecursive_square(fs, lambda);
    fRecursive(fs, lambda + 1);
    fRecursive(fs, lambda + 2);

    fmpz_set_ui(b, 4);
    fmpz_mod_poly_scalar_mul_fmpz(tmp, tmp, b);
    fmpz_mod_poly_mulmod_preinv(tmp, tmp, f_square[lambda], h, hinv);
    fmpz_mod_poly_mulmod_preinv(tmp, tmp, f[lambda], h, hinv);

    if ((lambda & 1) == 0)
        fmpz_mod_poly_mulmod_preinv(tmp, tmp, square_2y, h, hinv);

    fmpz_mod_poly_mulmod_preinv(tmp2, f[lambda + 2], f[lambda - 1], h,
                                    hinv);
    fmpz_mod_poly_mulmod_preinv(tmp2, tmp2, f[lambda - 1], h, hinv);
    if ((lambda & 1) == 1)
        fmpz_mod_poly_mulmod_preinv(tmp2, tmp2, square_2y, h, hinv);

    fmpz_mod_poly_sub(tmp, tmp, tmp2);

    /* if f[lambda-2] not defined */
    if (lambda == 1)
        fmpz_mod_poly_neg(tmp2, f[lambda + 1]);
    else
    {
        fRecursive(fs, lambda - 2);
        fmpz_mod_poly_mulmod_preinv(tmp2, f[lambda - 2], f[lambda + 1],
                                        h, hinv);
    }

    fmpz_mod_poly_mulmod_preinv(tmp2, tmp2, f[lambda + 1], h, hinv);

    if ((lambda & 1) == 1)
        fmpz_mod_poly_mulmod_preinv(tmp2, tmp2, square_2y, h, hinv);

    fmpz_mod_poly_add(tmp, tmp, tmp2);
    fmpz_mod_poly_gcd(tmp, tmp, h);

    if (fmpz_mod_poly_degree(tmp) < 1)
    {
        tl = l - tl;
        *eigenvalue = l - lambda;
    }
    else
        *eigenvalue = lambda;

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_mod_poly_clear(tmp);
    fmpz_mod_poly_clear(tmp2);
    return tl;
}


/*
 * Checks if t mod l is tl_guess or l-tl_guess
 *
 */
int
findSignOfTrace(int *eigenvalue, const DivisionPolynomials * divPolys,
                const int l,
                const fmpz_mod_poly_t h,
                const int lambda_guess,
                const int tl_guess)
{
    int lambda = lambda_guess;
    int tl = tl_guess;

    if (l % 4 == 3)
    {
        fmpz_mod_poly_t square_y;
        fmpz_mod_poly_init(square_y, divPolys->curve->p);

        if (fmpz_mod_poly_degree(h) > 3) 
            fmpz_mod_poly_set(square_y, divPolys->square_y);
        else
            fmpz_mod_poly_rem(square_y, divPolys->square_y, h);

        tl = findSignOfTrace_3mod4(eigenvalue, divPolys,l,h,square_y,lambda,tl);
        fmpz_mod_poly_clear(square_y);
    }
    else
    {
        int upper_bound = (l - 1) >> 1;
        int alloc = (l == 3) ? 5 : (upper_bound + 3);

        if (lambda > upper_bound)
        {
            lambda = l-lambda;
            tl = l-tl;
        }

        fmpz_mod_poly_t square_y, square_2y, quartic_2y;
        Fs *fs;

        fmpz_mod_poly_init(square_y, divPolys->curve->p);
        fmpz_mod_poly_init(square_2y, divPolys->curve->p);
        fmpz_mod_poly_init(quartic_2y, divPolys->curve->p);

        if (fmpz_mod_poly_degree(h) > 6)
        {
            fmpz_mod_poly_set(square_y, divPolys->square_y);
            fmpz_mod_poly_set(square_2y, divPolys->square_2y);
            fmpz_mod_poly_set(quartic_2y, divPolys->quartic_2y);
        }
        else
        {
            fmpz_mod_poly_rem(quartic_2y, divPolys->quartic_2y, h);
            if (fmpz_mod_poly_degree(h) > 3)
            {
                fmpz_mod_poly_set(square_y, divPolys->square_y);
                fmpz_mod_poly_set(square_2y, divPolys->square_2y);
            }
            else
            {
                fmpz_mod_poly_rem(square_y, divPolys->square_y, h);
                fmpz_mod_poly_rem(square_2y, divPolys->square_2y, h);
            }
        }
        fs = new_Fs(alloc, h, divPolys->curve->p, quartic_2y, divPolys);

        tl = findSignOfTrace_slow(eigenvalue, divPolys,l,h,square_y,square_2y,lambda,tl,fs);

        fmpz_mod_poly_clear(square_y);
        fmpz_mod_poly_clear(square_2y);
        fmpz_mod_poly_clear(quartic_2y);

        free_Fs(fs);
    }

    return tl;
}

int
findTrace(int *eigenvalue, const DivisionPolynomials * divPolys, const int l,
          const int n, const int lambda_prev, const fmpz_mod_poly_t h,
          int twist_secure)
{
    slong i;
    fmpz_t a, b;

    /* initialize database of division polynomials */

    int ln1 = n_pow(l, n - 1), ln = ln1 * l;
    int lambda = 0;
    int upper_bound = (n == 1) ? ((l - 1) >> 1) : ln - 1;
    int alloc = (ln == 3) ? 5 : (upper_bound + 3);
    int tl, pl;

    Fs *fs;

    fmpz_mod_poly_struct *hinv;


    fmpz_mod_poly_t *f, *f_square;

    fmpz_mod_poly_t xpx, tmp, tmp2, square_y, square_2y, quartic_2y;
    fmpz_mod_poly_init(square_y, divPolys->curve->p);
    fmpz_mod_poly_init(square_2y, divPolys->curve->p);
    fmpz_mod_poly_init(quartic_2y, divPolys->curve->p);

    if (TALKATIVE)
        tic(4);

    if (fmpz_mod_poly_degree(h) > 6)
    {
        fmpz_mod_poly_set(square_y, divPolys->square_y);
        fmpz_mod_poly_set(square_2y, divPolys->square_2y);
        fmpz_mod_poly_set(quartic_2y, divPolys->quartic_2y);
    }
    else
    {
        fmpz_mod_poly_rem(quartic_2y, divPolys->quartic_2y, h);
        if (fmpz_mod_poly_degree(h) > 3)
        {
            fmpz_mod_poly_set(square_y, divPolys->square_y);
            fmpz_mod_poly_set(square_2y, divPolys->square_2y);
        }
        else
        {
            fmpz_mod_poly_rem(square_y, divPolys->square_y, h);
            fmpz_mod_poly_rem(square_2y, divPolys->square_2y, h);
        }

    }

    fs = new_Fs(alloc, h, divPolys->curve->p, quartic_2y, divPolys);

    hinv = fs->hinv;

    f = fs->f;
    f_square = fs->f_square;

    fmpz_mod_poly_init(tmp, divPolys->curve->p);
    fmpz_mod_poly_init(tmp2, divPolys->curve->p);

    fmpz_mod_poly_zero(tmp);
    fmpz_mod_poly_set_coeff_ui(tmp, 1, 1);

    fmpz_mod_poly_init(xpx, divPolys->curve->p);

    /* BIG EXPONENTIATION */

    if (TALKATIVE)
        tic(7);

    powerXMod(xpx, divPolys->curve->p, h);

    if (TALKATIVE >= 3)
    {
        printf("x^p");
        tac(7);
    }

    fmpz_mod_poly_sub(xpx, xpx, tmp);

    /* if l = 3 */
    if (fmpz_mod_poly_degree(h) == 1)
        fmpz_mod_poly_rem(xpx, xpx, h);

    tic(8);
    for (i = (n == 1) ? 1 : lambda_prev; i <= upper_bound; i += ln1)
    {
        fRecursive(fs, i - 1);
        fRecursive_square(fs, i);
        fRecursive(fs, i + 1);

        fmpz_mod_poly_mulmod_preinv(tmp, f_square[i], xpx, h, hinv);
        fmpz_mod_poly_mulmod_preinv(tmp2, f[i - 1], f[i + 1], h, hinv);

        if ((i & 1) == 0)        /* i is even */
            fmpz_mod_poly_mulmod_preinv(tmp, tmp, square_2y, h, hinv);
        else                     /* i is odd */
            fmpz_mod_poly_mulmod_preinv(tmp2, tmp2, square_2y, h, hinv);

        fmpz_mod_poly_add(tmp, tmp, tmp2);
        fmpz_mod_poly_gcd(tmp, tmp, h);

        /* if (fmpz_mod_poly_degree(tmp) > 0) { */
        if (fmpz_mod_poly_equal(tmp, h))
        {
            lambda = i;
            break;
        }
    }

    if (TALKATIVE >= 3)
    {
        printf("eigen");
        tac(8);
    }

    if (lambda == 0)
    {
        free_Fs(fs);

        fmpz_mod_poly_clear(xpx);
        fmpz_mod_poly_clear(tmp);
        fmpz_mod_poly_clear(tmp2);
        fmpz_mod_poly_clear(square_y);
        fmpz_mod_poly_clear(square_2y);
        fmpz_mod_poly_clear(quartic_2y);

        if (TALKATIVE >= 3)
        {
            printf("waist");
            tac(4);
            fflush(stdout);
        }

        return -2;              /* error, should not be reached, unless did not check that */
                                /* h was a factor of the division polynomial */
    }

    fmpz_init_set_ui(a, ln);
    fmpz_init_set_ui(b, lambda);
    fmpz_invmod(a, b, a);
    fmpz_mul(a, a, divPolys->curve->p);
    fmpz_mod_ui(a, a, ln);
    fmpz_add_ui(a, a, lambda);
    tl = fmpz_mod_ui(a, a, ln);
    pl = fmpz_mod_ui(a, divPolys->curve->p, ln);

    /*
     * if l divides pl + 1 - tl, then either the curve or its twist has order
     * divisible by l
     */

    if (twist_secure && (((pl + 1 - tl) % l == 0) || ((pl + 1 + tl) % l == 0)))
    {

        fmpz_clear(a);
        fmpz_clear(b);

        free_Fs(fs);

        fmpz_mod_poly_clear(xpx);
        fmpz_mod_poly_clear(tmp);
        fmpz_mod_poly_clear(tmp2);
        fmpz_mod_poly_clear(square_y);
        fmpz_mod_poly_clear(square_2y);
        fmpz_mod_poly_clear(quartic_2y);

        return -1;              /* early abort */
    }

    //printf("\nGuess1: %d %d\n",lambda,tl);
    if (n != 1)
    {
        *eigenvalue = lambda;
    }

    /*
     * Now, need to check if t mod l is tl of l-tl !
     * Fast method if l % 4 == 3
     */
    else if (l % 4 == 3)
        tl = findSignOfTrace_3mod4(eigenvalue, divPolys,l,h,square_y,lambda,tl);
    else
        tl = findSignOfTrace_slow(eigenvalue, divPolys,l,h,square_y,square_2y,lambda,tl,fs);

    fmpz_clear(a);
    fmpz_clear(b);

    free_Fs(fs);

    fmpz_mod_poly_clear(xpx);
    fmpz_mod_poly_clear(tmp);
    fmpz_mod_poly_clear(tmp2);
    fmpz_mod_poly_clear(square_y);
    fmpz_mod_poly_clear(square_2y);
    fmpz_mod_poly_clear(quartic_2y);

    return tl;
}

void
print_DivisionPolynomials(const DivisionPolynomials * poly)
{
    /* fmpz_mod_poly_print_pretty(poly->f, "x"); */
}

void
free_DivisionPolynomials(DivisionPolynomials * poly)
{
    slong i;
    if (poly)
    {
        for (i = 0; i < 5; ++i)
            fmpz_mod_poly_clear(poly->f[i]);

        fmpz_mod_poly_clear(poly->square_y);
        fmpz_mod_poly_clear(poly->square_2y);
        fmpz_mod_poly_clear(poly->quartic_2y);
        free(poly->f);
        free(poly);
    }
}
