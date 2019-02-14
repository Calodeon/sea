
#include <stdio.h>
#include <string.h>

#include "DivisionPolynomial.h"
#include "EllipticCurvePointCounting.h"
#include "EllipticCurveIsogenies.h"
#include "fmpz_mat.h"
#include "fq.h"
#include "mat_rank_mod.h"
#include "perm.h"
#include "toolbox.h"


// TODO: get rid of these magic numbers
/*
 * To minimize collisions, should be about sqrt(the bound for BSGS)
 */

#define HASH_TABLE_SIZE ((1<<14) + 1)

/* 
 * Here fmpz_get_ui is expected to return the least significant limb of the
 * input fmpz...
 */

#define HASH(P) (((P)->is_zero)? 1<<14: (fmpz_get_ui(((P)->x)) & ((1<<14) - 1)))

/*
 * Atkin primes larger than this bound will be ignored
 */

#define LARGEST_ATKIN 150

static PrimeData *SEA_smallPrimesData = NULL;

/*
 * Information from Theorem 1 of "Elliptic Curves with the Montgomery-Form and
 * Their Cryptographic Applications"; 0 means order mod 8 == 0
 *
 * The binary representation of the indices are the boolean values of
 *   0001: p = 3 mod 4
 *   0010: A+2 is a quad. residue
 *   0100: A-2 is a quad. residue
 *   1000: B is a quad. residue
 */

static const int montgomery_mod8[16] =
    { 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0 };

/*
 * To be called once before the first time a point counting is performed (loads
 * the modular polynomials)
 */

void
SEA_init()
{
    slong i;
    if (SEA_smallPrimesData == NULL)
    {
        SEA_smallPrimesData = SEA_init_smallPrimes(MAX_NBR_SMALL_PRIMES);
        if (SEA_smallPrimesData != NULL)
        {
            for (i = 0; i < INIT_NBR_MODULAR_POLY; ++i)
            {
                SEA_smallPrimesData[i].modular = read_ModularPolynomialAsymmetric(SEA_smallPrimesData[i].l, PATH_ATKIN);
            }
        }
    }
}

/*
 * Frees the memory (the modular polynomials)
 */

void
SEA_clear()
{
    if (SEA_smallPrimesData)
        SEA_free_smallPrimes(SEA_smallPrimesData, MAX_NBR_SMALL_PRIMES);
}



PrimeData *
SEA_init_smallPrimes(int nbr_small_primes)
{
    slong i;
    FILE *file;

    PrimeData *primeData = malloc(sizeof(PrimeData) * nbr_small_primes);

    file = fopen(PATH_SMALL_PRIMES, "r");
    if (!file)
    {
        printf("error opening file \"small_primes.txt\"\n");
        free(primeData);
        return NULL;
    }

    for (i = 0; i < nbr_small_primes; ++i)
    {
        int success = SEA_read_primeData(primeData + i, file);
        if (!success) {
            free(primeData);
            return NULL;
        }
    }

    fclose(file);

    return primeData;
}

int
SEA_read_primeData(PrimeData * dest, FILE *file)
{
    int status =
        fscanf(file, "%d %d (%d, %d)", &(dest->l), &(dest->nonSquare),
               &(dest->generator1), &(dest->generator2));

    int success = fscanf(file, "\n");

    if (success < 0 || status < 4)
    {
        printf("Not enough entries in \"small_primes.txt\"\n");
        return 0;
    }

    dest->modular = NULL;

    return 1;
}

void
SEA_free_smallPrimes(PrimeData * data, int nbr_small_primes)
{
    slong i;
    for (i = 0; i < nbr_small_primes; ++i) {
        if(data[i].modular){
            free_ModularPolynomial(data[i].modular);
        }
    }

    free(data);
}


/*
 * Simple interface for GMP
 */

int
SEA_gmp(mpz_t n, const mpz_t p, const mpz_t A, const mpz_t B)
{
    fmpz_t flint_n, flint_p, flint_A, flint_B;
    int feedback;

    fmpz_init(flint_n);
    fmpz_init(flint_p);
    fmpz_init(flint_A);
    fmpz_init(flint_B);

    fmpz_set_mpz(flint_p, p);
    fmpz_set_mpz(flint_A, A);
    fmpz_set_mpz(flint_B, B);

    feedback = SEA_simple(flint_n, flint_p, flint_A, flint_B);

    fmpz_get_mpz(n, flint_n);

    fmpz_clear(flint_n);
    fmpz_clear(flint_p);
    fmpz_clear(flint_A);
    fmpz_clear(flint_B);

    return feedback;
}

/*
 * Simple interface
 */

int
SEA_simple(fmpz_t n, const fmpz_t p, const fmpz_t A, const fmpz_t B)
{
    EllipticCurve *curve = new_EllipticCurve(p, A, B);
    int feedback;
    flint_rand_t state;

    flint_randinit(state);

    feedback = SEA_AtkinPolynomials(n, curve, SEA_smallPrimesData, state,
                                        PRIME, 1);

    free_EllipticCurve(curve);
    flint_randclear(state);

    return feedback;
}

/*
 * Full interface
 */

int
SEA(fmpz_t n, const fmpz_t p, const fmpz_t A, const fmpz_t B,
    AbortStrategy strategy, int twist_secure)
{
    EllipticCurve *curve = new_EllipticCurve(p, A, B);
    int feedback;
    flint_rand_t state;

    flint_randinit(state);

    feedback = SEA_AtkinPolynomials(n, curve, SEA_smallPrimesData, state,
                                        strategy, twist_secure);

    free_EllipticCurve(curve);
    flint_randclear(state);

    return feedback;
}

/*
 * Classification of splitting type of modular polynomial
 *
 *   ELKIES: (1,1,r,r,...,r)
 *   ATKIN: (r,r,...,r)
 *   COMPLETELY_SPLIT: (1,1,...,1)
 *   SINGLE_ROOT: (1,l)
 *   DROP: not interesting (large atkin)
 */

typedef enum
{ ELKIES, ATKIN, COMPLETELY_SPLIT, SINGLE_ROOT, DROP } prime_type;

typedef struct AtkinPrime AtkinPrime;
struct AtkinPrime
{
    int l;                      /* The prime */
    int N;                      /* Number of candidates for t mod l */
    int *tl;                    /* The candidates */
    AtkinPrime *next;           /* Next in the list */
};

AtkinPrime *
new_AtkinPrime(int l, AtkinPrime * next)
{
    AtkinPrime *prime = malloc(sizeof(AtkinPrime));
    prime->l = l;
    prime->N = 0;
    prime->tl = malloc(sizeof(int) * l);
    prime->next = next;
    return prime;
}

void
free_AtkinPrime(AtkinPrime * prime)
{
    free(prime->tl);
    free(prime);
}

void
free_AtkinPrime_deep(AtkinPrime * prime)
{
    AtkinPrime *ptr = prime;
    while (ptr)
    {
        prime = ptr;
        ptr = prime->next;
        free_AtkinPrime(prime);
    }
}

void
print_AtkinPrime_deep(AtkinPrime * prime)
{
    while (prime)
    {
        printf("Atkin[%d] -> %d\n", prime->l, prime->N);
        prime = prime->next;
    }
}

void
atkin_addCandidate(AtkinPrime * prime, int n)
{
    prime->tl[prime->N] = n;
    prime->N += 1;
}

typedef struct BSGSList BSGSList;
struct BSGSList
{
    fmpz_t ta;                  /* The candidate trace mod some atkin primes */
    fmpz_t ma;                  /* The product of those atkin primes */
    fmpz_t r;                   /* The corresponding r = (ta - t_elkies) / (m_elkies * ma) */
    ECPoint *P;
    BSGSList *next;
};

/* Order the lists by increasing order of "r" */
BSGSList *
orderBSGSListByRs(BSGSList * L)
{
    BSGSList *A, *B, *elem, *tail;
    int inSize, nbrMerges, aSize, bSize, i;

    if (!L)
        return NULL;

    inSize = 1;                 /* Size of the merging lists */

    while (1)
    {
        A = L;
        L = NULL;
        tail = NULL;

        nbrMerges = 0;

        while (A)
        {

            nbrMerges++;

            B = A;
            aSize = 0;
            for (i = 0; i < inSize; i++)
            {
                aSize++;
                B = B->next;
                if (!B)
                    break;
            }

            bSize = inSize;

            /* Merging the two lists */
            while (aSize > 0 || (bSize > 0 && B))
            {

                /* Find the next element */
                if (aSize == 0)
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }
                else if (bSize == 0 || !B)
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else if (fmpz_cmp(A->r, B->r) < 0)
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }

                /* Add it */
                if (tail)
                {
                    tail->next = elem;
                }
                else
                {
                    L = elem;
                }
                tail = elem;
            }

            A = B;
        }

        tail->next = NULL;

        /* One merge => we are done */
        if (nbrMerges <= 1)
            return L;

        /* If not done yet, merge lists twice the size */
        inSize *= 2;
    }

}

/* Insert L in L->next in O(n) */
AtkinPrime *
insertAtkinPrimes(AtkinPrime * L)
{
    if (!L)
        return NULL;            /* Empty list */
    else {

        AtkinPrime *tail = L->next;
        AtkinPrime *current = L->next;
        AtkinPrime *prev = NULL;

        /* If single element, return it */
        if (!current)
            return L;

        /* Browse list */
        while (current)
        {

            /* If score is bigger than the currently browsed element */
            if ((L->l) * (current->N) > (current->l) * (L->N))
            {
                if (current == tail)
                    tail = L;       /* Inserted element has the biggest score */

                L->next = current;
                if (prev)
                    prev->next = L;

                return tail;
            }

            prev = current;
            current = current->next;
        }

        /* Inserted element has the lowest score */
        prev->next = L;
        L->next = NULL;

        return tail;
    }
}

/* Order the list according to there score l/candidates */
AtkinPrime *
orderAtkinPrimesByScore(AtkinPrime * L)
{
    AtkinPrime *A, *B, *elem, *tail;
    int inSize, nbrMerges, aSize, bSize, i;

    if (!L)
        return NULL;

    inSize = 1;

    while (1)
    {
        A = L;
        L = NULL;
        tail = NULL;

        nbrMerges = 0;

        while (A)
        {

            nbrMerges++;

            B = A;
            aSize = 0;
            for (i = 0; i < inSize; i++)
            {
                aSize++;
                B = B->next;
                if (!B)
                    break;
            }

            bSize = inSize;

            /* Merging the two lists */
            while (aSize > 0 || (bSize > 0 && B))
            {

                /* Find the next element */
                if (aSize == 0)
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }
                else if (bSize == 0 || !B)
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else if ((A->l) * B->N >= A->N * (B->l))
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }

                /* Add it */
                if (tail)
                {
                    tail->next = elem;
                }
                else
                {
                    L = elem;
                }
                tail = elem;
            }

            A = B;
        }

        tail->next = NULL;

        /* One merge => we are done */
        if (nbrMerges <= 1)
            return L;

        /* If not done yet, merge lists twice the size */
        inSize *= 2;
    }

}

/*
 * Order the list from the prime with largest number of candidates 
 * to the smallest
 */

AtkinPrime *
orderAtkinPrimesByCandidates(AtkinPrime * L)
{
    AtkinPrime *A, *B, *elem, *tail;
    int inSize, nbrMerges, aSize, bSize, i;

    if (!L)
        return NULL;

    inSize = 1;

    while (1)
    {
        A = L;
        L = NULL;
        tail = NULL;

        nbrMerges = 0;

        while (A)
        {

            nbrMerges++;

            B = A;
            aSize = 0;
            for (i = 0; i < inSize; i++)
            {
                aSize++;
                B = B->next;
                if (!B)
                    break;
            }

            bSize = inSize;

            /* Merging the two lists */
            while (aSize > 0 || (bSize > 0 && B))
            {

                /* Find the next element */
                if (aSize == 0)
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }
                else if (bSize == 0 || !B)
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else if (A->N >= B->N)
                {
                    elem = A;
                    A = A->next;
                    aSize--;
                }
                else
                {
                    elem = B;
                    B = B->next;
                    bSize--;
                }

                /* Add it */
                if (tail)
                {
                    tail->next = elem;
                }
                else
                {
                    L = elem;
                }
                tail = elem;
            }

            A = B;
        }

        tail->next = NULL;

        /* One merge => we are done */
        if (nbrMerges <= 1)
            return L;

        /* If not done yet, merge lists twice the size */
        inSize *= 2;
    }

}

/*
 * Split L0 in two lists L1 and L2 such that the number of possible traces
 * mod the primes of L1 and L2 are close.
 * base2transfering will be set to the number of bits which can be removed
 * from K, by adding residues modulo powers of two on the side of ptrL1.
 * For example if K to big and the weights of the lists are unbalanced, e.g.
 * (15, 2^16), we can rebalance it as (15*2^6, 2^10), here base2transfering = 6
 */

void
evenlySplitAtkinPrimes(int *base2transfering, AtkinPrime ** ptrL1,
                       AtkinPrime ** ptrL2, AtkinPrime * L0, fmpz_t K)
{
    AtkinPrime *L1 = NULL, *L2 = NULL, *tmp;

    fmpz_t nL1, nL2;

    fmpz_init_set_ui(nL1, 1);
    fmpz_init_set(nL2, K);

    L0 = orderAtkinPrimesByCandidates(L0);

    while (L0)
    {
        if (fmpz_cmp(nL1, nL2) < 0)
        {
            tmp = L0;
            L0 = L0->next;
            tmp->next = L1;
            L1 = tmp;
            fmpz_mul_ui(nL1, nL1, L1->N);
        }
        else
        {
            tmp = L0;
            L0 = L0->next;
            tmp->next = L2;
            L2 = tmp;
            fmpz_mul_ui(nL2, nL2, L2->N);
        }
    }

    *base2transfering = 0;

    if (fmpz_cmp(nL1, nL2) < 0)
    {
        fmpz_t ratio;

        fmpz_init(ratio);

        fmpz_fdiv_q(ratio, nL2, nL1);
        fmpz_sqrt(ratio, ratio);

        *base2transfering = fmpz_sizeinbase(ratio, 2);

        fmpz_clear(ratio);
    }

    if (TALKATIVE >= 2)
    {
        fmpz_t prod;

        fmpz_init(prod);
        fmpz_mul(prod, nL1, nL2);

        printf("[");
        fmpz_print(prod);
        printf(" candidates]\n");

        fmpz_clear(prod);
    }
    if (TALKATIVE >= 2)
    {
        printf("(");
        fmpz_print(nL1);
        printf(" : ");
        fmpz_print(nL2);
        printf(") ");
    }

    *ptrL1 = L1;
    *ptrL2 = L2;

    fmpz_clear(nL1);
    fmpz_clear(nL2);
}

BSGSList *
new_BSGSList(BSGSList * trace, int tl, int l, BSGSList * next)
{
    BSGSList *bsgs = malloc(sizeof(BSGSList));

    fmpz_init(bsgs->r);
    fmpz_init(bsgs->ta);
    fmpz_init(bsgs->ma);
    bsgs->P = new_ECPoint();

    fmpz_CRT_ui(bsgs->ta, trace->ta, trace->ma, tl, l, 1);
    fmpz_mul_ui(bsgs->ma, trace->ma, l);

    bsgs->next = next;
    return bsgs;
}

BSGSList *
new_BSGSList_fromOneAtkinPrime(const AtkinPrime * atkin)
{
    BSGSList *bsgs = NULL, *ptr = NULL;
    slong i;

    for (i = 0; i < atkin->N; ++i)
    {
        ptr = bsgs;

        bsgs = malloc(sizeof(BSGSList));

        fmpz_init(bsgs->r);
        bsgs->P = new_ECPoint();

        fmpz_init_set_ui(bsgs->ta, (long) (atkin->tl[i]));
        fmpz_init_set_ui(bsgs->ma, (long) (atkin->l));
        bsgs->next = ptr;
    }

    return bsgs;
}

void
free_BSGSList(BSGSList * bsgs)
{

    fmpz_clear(bsgs->r);
    fmpz_clear(bsgs->ta);
    fmpz_clear(bsgs->ma);
    free_ECPoint(bsgs->P);
    free(bsgs);
}

void
free_BSGSList_deep(BSGSList * bsgs)
{
    BSGSList *ptr = bsgs;
    while (ptr)
    {
        bsgs = ptr;
        ptr = bsgs->next;
        free_BSGSList(bsgs);
    }
}

void
computeRInBSGSListCentered(BSGSList * bsgs, const fmpz_t t_elkies,
                           const fmpz_t m_elkies, const fmpz_t ma2)
{
    fmpz_t tmp, maDiv2;

    fmpz_init(tmp);
    fmpz_init(maDiv2);

    fmpz_mul(tmp, ma2, m_elkies);
    fmpz_invmod(tmp, tmp, bsgs->ma);
    fmpz_tdiv_q_2exp(maDiv2, bsgs->ma, 1);

    while (bsgs)
    {
        fmpz_sub(bsgs->r, bsgs->ta, t_elkies);
        fmpz_mul(bsgs->r, bsgs->r, tmp);
        fmpz_mod(bsgs->r, bsgs->r, bsgs->ma);

        if (fmpz_cmp(bsgs->r, maDiv2) > 0)
            fmpz_sub(bsgs->r, bsgs->r, bsgs->ma);

        bsgs = bsgs->next;
    }

    fmpz_clear(tmp);
    fmpz_clear(maDiv2);
}

void
computeRInBSGSListPositive(BSGSList * bsgs, const fmpz_t t_elkies,
                           const fmpz_t m_elkies, const fmpz_t ma1)
{
    fmpz_t tmp;
    fmpz_init(tmp);

    fmpz_mul(tmp, ma1, m_elkies);
    fmpz_invmod(tmp, tmp, bsgs->ma);

    while (bsgs)
    {
        fmpz_sub(bsgs->r, bsgs->ta, t_elkies);
        fmpz_mul(bsgs->r, bsgs->r, tmp);
        fmpz_mod(bsgs->r, bsgs->r, bsgs->ma);
        bsgs = bsgs->next;
    }

    fmpz_clear(tmp);
}

/*
 * Build a list of all the candidate traces for the list L of Atkin primes
 */

BSGSList *
buildBSGSListOfCandidateTraces_rec(const AtkinPrime * L)
{

    if (!L)
        return NULL;

    else if (!(L->next))
        return new_BSGSList_fromOneAtkinPrime(L);

    else {

        slong i;
        BSGSList *smaller_bsgs = buildBSGSListOfCandidateTraces_rec(L->next);
        BSGSList *extended = NULL, *smaller_iterator;
        fmpz_t new_ma, ma_ma_inv, l_l_inv, tmp, maDiv2;
        int next_exists;

        for (i = 0; i < L->N - 1; ++i)
        {
            smaller_iterator = smaller_bsgs;
            while (smaller_iterator)
            {
                extended =
                    new_BSGSList(smaller_iterator, L->tl[i], L->l, extended);
                smaller_iterator = smaller_iterator->next;
            }
        }

        smaller_iterator = smaller_bsgs;
        fmpz_init(new_ma);
        fmpz_init(tmp);
        fmpz_init(ma_ma_inv);
        fmpz_init(l_l_inv);
        fmpz_init(maDiv2);

        fmpz_mul_ui(new_ma, smaller_bsgs->ma, L->l);

        fmpz_set_ui(tmp, L->l);
        fmpz_invmod(ma_ma_inv, smaller_bsgs->ma, tmp);
        fmpz_mul(ma_ma_inv, ma_ma_inv, smaller_bsgs->ma);
        fmpz_invmod(l_l_inv, tmp, smaller_bsgs->ma);
        fmpz_mul_ui(l_l_inv, l_l_inv, L->l);

        fmpz_tdiv_q_2exp(maDiv2, smaller_bsgs->ma, 1);

        next_exists = 1;
        while (next_exists)
        {
            fmpz_mul(tmp, l_l_inv, smaller_iterator->ta);
            fmpz_mul_ui(smaller_iterator->ta, ma_ma_inv, L->tl[L->N - 1]);
            fmpz_add(smaller_iterator->ta, smaller_iterator->ta, tmp);
            fmpz_mod(smaller_iterator->ta, smaller_iterator->ta, new_ma);

            if (fmpz_cmp(smaller_iterator->ta, maDiv2) > 0)
                fmpz_sub(smaller_iterator->ta, smaller_iterator->ta, new_ma);

            fmpz_set(smaller_iterator->ma, new_ma);

            if (!(smaller_iterator->next))
            {
                smaller_iterator->next = extended;
                next_exists = 0;
            }
            else
            {
                smaller_iterator = smaller_iterator->next;
            }
        }

        fmpz_clear(new_ma);
        fmpz_clear(tmp);
        fmpz_clear(ma_ma_inv);
        fmpz_clear(l_l_inv);
        fmpz_clear(maDiv2);

        return smaller_bsgs;
    }
}

/*
 * Build a list of all the candidate traces for the list L of Atkin primes
 */

BSGSList *
buildBSGSListOfCandidateTraces(const AtkinPrime * L)
{

    BSGSList *bsgs = buildBSGSListOfCandidateTraces_rec(L);

    if (L == NULL)
    {
        bsgs = malloc(sizeof(BSGSList));

        fmpz_init(bsgs->r);
        bsgs->P = new_ECPoint();

        fmpz_init_set_ui(bsgs->ta, 0);
        fmpz_init_set_ui(bsgs->ma, 1);
        bsgs->next = NULL;
    }

    return bsgs;
}

static void
fmpz_mod_poly_to_fmpz_mat_col(fmpz_mat_t mat, slong col, fmpz_mod_poly_t poly)
{
    slong i;

    for (i = 0; i < poly->length; i++)
        fmpz_set(fmpz_mat_entry(mat, i, col), poly->coeffs + i);

    for (; i < mat->r; i++)
        fmpz_zero(fmpz_mat_entry(mat, i, col));

}


prime_type
splittingModularPolynomial(fmpz_mod_poly_t splittingPart, int *r,
                           const fmpz_mod_poly_t phi, const fmpz_t p)
{

    slong k;
    fmpz_mod_poly_t polynomial, x_p;
    const slong l = fmpz_mod_poly_degree(phi) - 1, n = l+1;

    fmpz_mod_poly_init(polynomial, p);
    fmpz_mod_poly_zero(polynomial);
    fmpz_mod_poly_set_coeff_ui(polynomial, 1, 1);   /* polynomial = x */

    fmpz_mod_poly_init(x_p, p);

    /* BIG EXPONENTIATION */

    if (TALKATIVE)
        tic(1);

    powerXMod(x_p, p, phi);     /* x_p = x^p */

    if (TALKATIVE >= 3){
        printf("x^p_mod"); tac(1);
    }

    fmpz_mod_poly_sub(polynomial, x_p, polynomial); /* polynomial = x^p - x */
    fmpz_mod_poly_gcd(polynomial, polynomial, phi);

    switch (fmpz_mod_poly_degree(polynomial))
    {
        case 0:
        {
            /*
             * Find the number of irreducible factors
             * Compute the matrix for the Berlekamp Map
             */

            if (TALKATIVE >= 2)
            {
                printf("Atkin\t");
                fflush(stdout);
            }

            if (l > LARGEST_ATKIN)
            {                   /* l is too large */
                fmpz_mod_poly_set(splittingPart, polynomial);
                fmpz_mod_poly_clear(x_p);
                fmpz_mod_poly_clear(polynomial);
                return DROP;
            }
            else {

                fmpz_mod_poly_t x_pi, x_pi2;
                fmpz_mat_t matrix;
                fmpz_t coeff, p_1;
                fmpz_mod_poly_t finv;
                int rank;


                tic(5);
                fmpz_init(coeff);
                fmpz_init(p_1);

                fmpz_mat_init(matrix, n, n);
                fmpz_mod_poly_init(x_pi, p);
                fmpz_mod_poly_init(x_pi2, p);
                fmpz_mod_poly_set_coeff_ui(x_pi, 0, 1);

                /* p_1 = p-1 */
                fmpz_init_set(p_1, p);
                fmpz_sub_ui(p_1, p_1, 1);
                fmpz_mod(p_1, p_1, p);

                fmpz_mod_poly_init(finv, fmpz_mod_poly_modulus(phi));
                fmpz_mod_poly_reverse(finv, phi, phi->length);
                fmpz_mod_poly_inv_series_newton(finv, finv, phi->length);

                for (k = 0; k < n; k++)
                {
                    /* Q - I */
                    fmpz_mod_poly_set(x_pi2, x_pi);
                    fmpz_mod_poly_get_coeff_fmpz(coeff, x_pi2, k);
                    if (!fmpz_is_zero(coeff))
                    {
                        fmpz_sub_ui(coeff, coeff, 1);
                        fmpz_mod(coeff, coeff, p);
                        fmpz_mod_poly_set_coeff_fmpz(x_pi2, k, coeff);
                    }
                    else
                        fmpz_mod_poly_set_coeff_fmpz(x_pi2, k, p_1);

                    fmpz_mod_poly_to_fmpz_mat_col(matrix, k, x_pi2);
                    fmpz_mod_poly_mulmod_preinv(x_pi, x_pi, x_p, phi, finv);
                }

                fmpz_mod_poly_clear(finv);
                fmpz_mod_poly_clear(x_pi);
                fmpz_mod_poly_clear(x_pi2);

                if (TALKATIVE >= 3)
                {
                    printf("\tmatrix");
                    tac(5);
                }

                /* Row reduce Q - I */
                if (TALKATIVE)
                    tic(5);

                rank = mat_rank_mod(matrix, p);
                *r = n / (n - rank);

                if (TALKATIVE >= 3)
                {
                    printf("rank");
                    tac(5);
                    printf("r = %d ", *r);
                }

                fmpz_clear(coeff);
                fmpz_clear(p_1);
                fmpz_mat_clear(matrix);

                fmpz_mod_poly_clear(x_p);
                fmpz_mod_poly_clear(polynomial);

                return ATKIN;
            }
        }

        case 1:
            fmpz_mod_poly_set(splittingPart, polynomial);
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return SINGLE_ROOT;

        case 2:
            fmpz_mod_poly_set(splittingPart, polynomial);
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return ELKIES;

        default:
            fmpz_mod_poly_set(splittingPart, polynomial);
            fmpz_mod_poly_clear(x_p);
            fmpz_mod_poly_clear(polynomial);
            return COMPLETELY_SPLIT;

    }

}


void findJfromF(fmpz_t j2,
                fmpz_mod_poly_t h,
                const fmpz_t j1,
                const fmpz_t f1,
                const fmpz_mod_poly_t psi_j1,
                const ModularPolynomial * modular,
                const EllipticCurve * curve,
                const int computeKernelPolynomialH)
{
    fmpz_mod_poly_t psi_f1, h1;
    fmpz_t new_j, A2,B2;
    fmpz *roots = NULL;
    int r, k, counter_j1, computed_h = 0, correct_j = 0;

    fmpz_init(new_j);
    fmpz_init(A2);
    fmpz_init(B2);
    fmpz_mod_poly_init(h1, curve->p);

    fmpz_mod_poly_init(psi_f1, curve->p);
    evaluate_ModularPolynomialAtX(psi_f1, modular, f1);

    /* Find j2 */
    r = rootsInFp(&roots, psi_f1, curve->p);

    /* For j2, we need to try ALL the roots of this polynomial... */
    k = 0;
    counter_j1 = 0;

    /* Try all the possible roots for j2 until the right one is found */
    while (1)
    {
        fmpz_set(new_j, roots + k);

        /* if j2 is just j1, try the next root... unless j1 has already been 
         * met: it is possible that j1 it its own neighbor
         */
        if (fmpz_equal(j1, new_j) && counter_j1 == 0)
        {
            ++counter_j1;
            ++k;
            continue;
        }

        /* if there are at most 2 roots, one of them is j1 and the other one is necessarily correct */
        if (r <= 2)
            break;
        
        /* if there are more than two roots, then they really need to be tested... the
         * following is suboptimal, but does not have a big impact anyway
         */

        correct_j = computeIsogeny_AtkinPolynomials(h1, A2, B2, curve, j1, new_j, f1,
                                             modular->l, psi_j1, psi_f1,
                                             modular, (r > 2) && (k != r - 1));
        if (!correct_j)
        {
            ++k;
            continue;
        }

        computed_h  = 1;

        break;
    }

    fmpz_set(j2, new_j);

    if (computeKernelPolynomialH)
    {
        if (!computed_h)
            computeIsogeny_AtkinPolynomials(h, A2, B2, curve, j1, j2, f1,
                                             modular->l, psi_j1, psi_f1,
                                             modular, 0);
        else
            fmpz_mod_poly_set(h, h1);
    }



    for (int i = 0; i < r; ++i)
        fmpz_clear(roots + i);
    free(roots);
    fmpz_mod_poly_clear(psi_f1);
    fmpz_mod_poly_clear(h1);
    fmpz_clear(new_j);
    fmpz_clear(A2);
    fmpz_clear(B2);
}

int
findJNeighbors_splittingDone(fmpz ** neighbors,
             const fmpz_mod_poly_t splittingPart,
             const ModularPolynomial * modular,
             const fmpz_mod_poly_t psi_j1,
             const fmpz_t j1,
             const EllipticCurve * curve)
{
    fmpz_t f1;
    fmpz_t j2;

    fmpz_mod_poly_factor_t fac;

    fmpz_init(f1);
    fmpz_init(j2);

    /* find f1 -- any of the possible f1's works (each of
     * them corresponds to an isogenous curve) */
    fmpz_mod_poly_factor_init(fac);

    fmpz_mod_poly_factor(fac, splittingPart);


    int d = fmpz_mod_poly_degree(splittingPart);

    if (d != 0)
    {
        fmpz *find_roots = malloc(sizeof(fmpz) * d);

        for (int i = 0; i < d; ++i)
        {
            fmpz_mod_poly_get_coeff_fmpz(f1, fac->poly + i, 0);
            fmpz_neg(f1, f1);
            fmpz_mod(f1, f1, curve->p);
            
            findJfromF(j2,NULL, j1, f1, psi_j1, modular, curve, 0);

            fmpz_init(find_roots + i);
            fmpz_set(find_roots + i, j2);

        }
        *neighbors = find_roots;
    }

    fmpz_clear(f1);
    fmpz_clear(j2);
    fmpz_mod_poly_factor_clear(fac);
    return d;
}

int
findJNeighbors(fmpz ** neighbors,
             const ModularPolynomial * modular,
             const fmpz_t j1,
             const EllipticCurve * curve)
{
    fmpz_mod_poly_t polynomial, x_p, psi_j1;

    if (TALKATIVE)
        tic(13);
    fmpz_mod_poly_init(psi_j1, curve->p);
    evaluate_ModularPolynomialAtY(psi_j1, modular, j1);

    fmpz_mod_poly_init(polynomial, curve->p);
    fmpz_mod_poly_zero(polynomial);
    fmpz_mod_poly_set_coeff_ui(polynomial, 1, 1);   /* polynomial = x */

    fmpz_mod_poly_init(x_p, curve->p);
    /* BIG EXPONENTIATION */


    powerXMod(x_p, curve->p, psi_j1);     /* x_p = x^p mod psi_j1 */


    fmpz_mod_poly_sub(polynomial, x_p, polynomial); /* polynomial = x^p - x */
    fmpz_mod_poly_gcd(polynomial, polynomial, psi_j1);
    if (TALKATIVE >= 3){
        printf("neighPre"); tac(13);
    }
    if (TALKATIVE)
        tic(13);
    int d = findJNeighbors_splittingDone(neighbors,polynomial,modular,psi_j1,j1,curve);
    if (TALKATIVE >= 3){
        printf("neighFind"); tac(13);
    }
    fmpz_mod_poly_clear(x_p);
    fmpz_mod_poly_clear(polynomial);
    fmpz_mod_poly_clear(psi_j1);
    return d;
}

int
findOneJNeighbor(fmpz_t j2,
                 fmpz_mod_poly_t h,
                 const fmpz_mod_poly_t splittingPart,
                 const ModularPolynomial * modular,
                 const fmpz_mod_poly_t psi_j1,
                 const fmpz_t j1,
                 const EllipticCurve * curve)
{
    fmpz_t f1;

    fmpz_mod_poly_factor_t fac;

    fmpz_init(f1);

    /* find f1 -- any of the possible f1's works (each of
     * them corresponds to an isogenous curve) */
    fmpz_mod_poly_factor_init(fac);

    fmpz_mod_poly_factor(fac, splittingPart);


    int d = fmpz_mod_poly_degree(splittingPart);

    if (d != 0)
    {   
        fmpz_mod_poly_get_coeff_fmpz(f1, fac->poly, 0);
        fmpz_neg(f1, f1);
        fmpz_mod(f1, f1, curve->p);
        
        findJfromF(j2, h,j1, f1, psi_j1, modular, curve,1);
    }

    fmpz_clear(f1);
    fmpz_mod_poly_factor_clear(fac);
    return d;
}

int
countJNeighbors(const ModularPolynomial * modular,
                const fmpz_t j1,
                const EllipticCurve * curve)
{
    fmpz_mod_poly_t polynomial, x_p, psi_j1;

    //if (TALKATIVE)
    //    tic(13);
    fmpz_mod_poly_init(psi_j1, curve->p);
    evaluate_ModularPolynomialAtY(psi_j1, modular, j1);

    fmpz_mod_poly_init(polynomial, curve->p);
    fmpz_mod_poly_zero(polynomial);
    fmpz_mod_poly_set_coeff_ui(polynomial, 1, 1);   /* polynomial = x */

    fmpz_mod_poly_init(x_p, curve->p);

    /* BIG EXPONENTIATION */
    powerXMod(x_p, curve->p, psi_j1);     /* x_p = x^p mod psi_j1 */


    fmpz_mod_poly_sub(polynomial, x_p, polynomial); /* polynomial = x^p - x */
    fmpz_mod_poly_gcd(polynomial, polynomial, psi_j1);

    int d = fmpz_mod_poly_degree(polynomial);

    fmpz_mod_poly_clear(x_p);
    fmpz_mod_poly_clear(polynomial);
    fmpz_mod_poly_clear(psi_j1);

    return d;
}

/*
 * Navigates in the isogeny volcano to determine how many
 * times l divides the discriminant of frobenius (actually,
 * just a lower bound for efficiency)
 */
int
valuationDiscriminant(const fmpz_mod_poly_t splittingPart,
             const ModularPolynomial * modular,
             const fmpz_mod_poly_t psi_j1,
             const fmpz_t j1,
             const EllipticCurve * curve)
{

    int depth = 0, ramifies = 0;
    /* if the isogeny volcano is a point or a cycle, depth 0*/
    if (fmpz_mod_poly_degree(splittingPart) == 0 || fmpz_mod_poly_degree(splittingPart) == 2)
    {
        depth = 0;
        ramifies = 0;
    }
    /* if there is a single neighbor, the curve could be at the bottom of
     * a non-trivial volcano, or the volcano could have depth 0 */
    else if (fmpz_mod_poly_degree(splittingPart) == 1)
    {

        if (TALKATIVE >= 3)
            tic(15);

        fmpz *neighbors1 = NULL;
        int n_neighbors1 = findJNeighbors_splittingDone(&neighbors1, splittingPart, modular, psi_j1, j1, curve);

        int ctr = countJNeighbors(modular, neighbors1, curve);
        
        if (ctr == 1)
        {
            depth = 0;
            ramifies = 1;
        }
        else
        {
            depth = 1;
            ramifies = 0;
        }

        for (int i = 0; i < n_neighbors1; ++i)
                fmpz_clear(neighbors1 + i);
        free(neighbors1);

        if (TALKATIVE >= 3){
            printf("depth"); tac(15);
        }
    }
    /* otherwise, there are l+1 neighbors, and the volcano has depth at least 1 */
    else
    {
        depth = 1;
        ramifies = 0;
    }

    return 2*depth + ramifies;
}

int
findTrace_Elkies_ramified(const int l,
                 const DivisionPolynomials * divPolys,
                 const fmpz_mod_poly_t h,
                 const EllipticCurve * curve, int twist_secure)
{
    int tl_guess;
    fmpz_t tmp;

    fmpz_init(tmp);

    int p_l = fmpz_mod_ui(tmp, curve->p, l);
    
    // TODO: PROBLEM: if p_l is not a square residue mod l, then tl_guess = 0, BUG

    tl_guess = n_sqrtmod((4*p_l)%l, l);

    int lambda = (tl_guess * n_invmod(2, l))%l;
    
    fmpz_clear(tmp);
    //printf("\nGuess2 [%d %d] data: %d %d\n",lambda,tl_guess,p_l,l);
    //tl =  findSignOfTrace(eigenvalue, divPolys,l,h,lambda,tl);
    tl_guess = findSignOfTrace(&lambda, divPolys, l, h, lambda, tl_guess);
    //printf("\nCorrect2 %d %d\n",lambda,tl_guess);
    return tl_guess;
}


int
findTrace_Elkies(const fmpz_mod_poly_t splittingPart,
                 const ModularPolynomial * modular, int n,
                 const DivisionPolynomials * divPolys,
                 const fmpz_mod_poly_t psi_j1, const fmpz_t j1,
                 const EllipticCurve * curve, int twist_secure)
{
    slong i,j;
    fmpz_t f1, tmp, A2, B2;
    fmpz_t j2;

    fmpz_mod_poly_t psi_f1, h;
    fmpz_mod_poly_factor_t fac;

    fmpz *roots = NULL;
    int r, lambda, tl, k, ctr_j1;

    fmpz_init(f1);
    fmpz_init(j2);
    fmpz_init(tmp);
    fmpz_init(A2);
    fmpz_init(B2);
    fmpz_mod_poly_init(psi_f1, curve->p);
    fmpz_mod_poly_init(h, curve->p);

    if (TALKATIVE >= 4)
        tic(19);

    /* find f1 -- any of the possible f1's works (each of
     * them corresponds to an isogenous curve) */

    //fmpz_mod_poly_factor_init(fac);
    // fmpz_mod_poly_factor(fac, splittingPart);

    // fmpz_mod_poly_get_coeff_fmpz(f1, fac->poly, 0);
    // fmpz_neg(f1, f1);
    // fmpz_mod(f1, f1, curve->p);
    //fmpz_mod_poly_factor_clear(fac);

    r = rootsInFp(&roots, splittingPart, curve->p);
    fmpz_set(f1, roots);

    for (i = 0; i < r; ++i)
        fmpz_clear(roots + i);
    if (roots) free(roots);

    evaluate_ModularPolynomialAtX(psi_f1, modular, f1);

    fmpz_mod_poly_factor_init(fac);

    /* Find j2 */
    if (fmpz_mod_poly_degree(psi_f1) == 2)
    {
        r = 1;
        roots = malloc(sizeof(fmpz));
        fmpz_init(roots);

        fmpz_mod_poly_get_coeff_fmpz(roots, psi_f1, 1);
        fmpz_add(roots, roots, j1);
        fmpz_neg(roots, roots);
        fmpz_mod(roots, roots, curve->p);
    }
    else
        r = rootsInFp(&roots, psi_f1, curve->p);

    /* For j2, we need to try ALL the roots of this polynomial... */
    lambda = -1;
    tl = -1;
    k = 0;
    ctr_j1 = 0;

    if (TALKATIVE >= 4)
    {
        printf("factor");tac(19);
    }

    /* Try all the possible roots for j2 until the right one is found */
    while (1)
    {

        fmpz_set(j2, roots + k);

        /* if j1 is a doule root, it could be the correct j-invariant of the
         * target (i.e., the isogeny is an endomorphism) */
        if (ctr_j1 == 0 && fmpz_equal(j1, j2))
        {
            ++ctr_j1;
            ++k;
            continue;
        }

        if (TALKATIVE >= 4)
            tic(18);
        if (!computeIsogeny_AtkinPolynomials(h, A2, B2, curve, j1, j2, f1,
                                             modular->l, psi_j1, psi_f1,
                                             modular, (r > 2) && (k != r - 1)))
        {
            ++k;
            continue;
            if (TALKATIVE >= 4)
            {
                printf("isog");tac(18);
            }
        }
        if (TALKATIVE >= 4)
        {
            printf("isog");tac(18);
        }

        /*
         * Might fail (i.e. returns tl = -2) if ((modular[i]->l > 31)
         * && (modular[i]->l != 41) && (modular[i]->l != 47)
         * && (modular[i]->l != 59) && (modular[i]->l != 71))
         * in case j2 is not the right root
         */

        tl = findTrace(&lambda, divPolys, modular->l, 1, 0, h, twist_secure);

        if (tl == -1)
        {
            /*
             * An early abort happenned in findTrace, meaning that
             * twist_secure == true and a small divisor was found for the
             * curve's order or its twist's
             */

            for (i = 0; i < r; ++i)
                fmpz_clear(roots + i);

            free(roots);
            fmpz_clear(f1);
            fmpz_clear(A2);
            fmpz_clear(B2);
            fmpz_clear(j2);
            fmpz_clear(tmp);
            fmpz_mod_poly_factor_clear(fac);
            fmpz_mod_poly_clear(psi_f1);
            fmpz_mod_poly_clear(h);

            return -1;
        }
        else if (tl < 0)
        {
            ++k;
            continue;
        }

        if (n > 1){

            /* Using isogeny cycles */
            /* Based on the paper "Isogeny cycles and the Schoof-Elkies-Atkin
             * algorithm", Couveignes, Dewaghe, Morain */

            fmpz_mod_poly_t psi_f2, psi_j2, h_i, phi, iso_psi, iso_phi, iso_psi_next, iso_phi_next;
            fmpz_t f2, j3, A3, B3;
            EllipticCurve *E2 = new_EllipticCurve(curve->p, A2, B2);

            fmpz_mod_poly_init(psi_f2, curve->p);
            fmpz_mod_poly_init(phi, curve->p);
            fmpz_mod_poly_init(h_i, curve->p);
            fmpz_mod_poly_init(iso_psi, curve->p);
            fmpz_mod_poly_init(iso_phi, curve->p);
            fmpz_mod_poly_init(iso_psi_next, curve->p);
            fmpz_mod_poly_init(iso_phi_next, curve->p);
            fmpz_mod_poly_init(psi_j2, curve->p);

            fmpz_init(f2);
            fmpz_init(j3);
            fmpz_init(A3);
            fmpz_init(B3);

            fmpz_mod_poly_set(iso_psi, h);

            isogenyPhiFromPsi(iso_phi, iso_psi, curve);

            for (i = 2; i <= n; ++i)
            {

                if (TALKATIVE >= 2)
                {
                    printf("\n [%d^%ld]\t", modular->l, i);
                }

                /* 1: Find the next curve */
                for (j = 0; j < r; ++j)
                    fmpz_clear(roots + j);
                free(roots);

                if (TALKATIVE >= 4)
                    tic(21);
                evaluate_ModularPolynomialAtY(psi_j2, modular, j2);

                r = rootsInFp(&roots, psi_j2, curve->p);

                /*
                 * Find the j-invariant of next curve in the chain (i.e., make sure
                 * not going backwards)
                 */

                if (fmpz_equal(f1, roots))
                    fmpz_set(f2, roots + 1);
                else
                    fmpz_set(f2, roots);

                for (j = 0; j < r; ++j)
                    fmpz_clear(roots + j);
                free(roots);


                evaluate_ModularPolynomialAtX(psi_f2, modular, f2);

                if (fmpz_mod_poly_degree(psi_f2) == 2)
                {
                    r = 1;
                    roots = malloc(sizeof(fmpz));
                    fmpz_init(roots);

                    fmpz_mod_poly_get_coeff_fmpz(roots, psi_f2, 1);
                    fmpz_add(roots, roots, j2);
                    fmpz_neg(roots, roots);
                    fmpz_mod(roots, roots, curve->p);
                }
                else
                    r = rootsInFp(&roots, psi_f2, curve->p);

                /*
                 * Isogeny cycles are only used for the smallest primes, in
                 * particular, there will be only two roots, no ambiguity !
                 * The first bad prime is 37, which is too big.
                 */

                if (fmpz_equal(j2, roots))
                    fmpz_set(j3, roots + 1);
                else
                    fmpz_set(j3, roots);

                if (TALKATIVE >= 4)
                {
                    printf("factor");tac(21);
                }

                if (TALKATIVE >= 4)
                    tic(22);
                /* 2: Compute factor h_i of kernel of the new l-isogeny */
                computeIsogeny_AtkinPolynomials(h_i, A3,
                                                B3, E2, j2, j3,
                                                f2, modular->l,
                                                psi_j2, psi_f2, modular, 0);

                /* 3: Set h to be the numerator of h_i o I_{n-1} o ... o I_1,
                 * where I_j is the j-th isogeny.
                 */
                
                fmpz_mod_poly_mul(h, iso_psi, iso_psi);
                numeratorRationalComposition(h, h_i, iso_phi, h);

                /*
                 * 4: find eigenvalue mod l^i, i.e., 0 <= tau < l such that
                 * lambda_n = lambda_{n-1} + tau*l^{n-1} is eigenvalue
                 */


                if (TALKATIVE >= 4)
                {
                    printf("isog");tac(22);
                }

                tl = findTrace(&lambda, divPolys, modular->l, i, lambda, h,
                            twist_secure);

                if (TALKATIVE && (tl < 0))
                {
                    printf("ERROR: Isogeny cycles failure\n");
                }

                if (i != n)
                {

                    if (TALKATIVE >= 4)
                        tic(23);
                    fmpz_mod_poly_mul(iso_psi_next, h, iso_psi);

                    isogenyPhiFromPsi(phi, h_i, E2);
                    fmpz_mod_poly_mul(h, iso_psi, iso_psi);
                    numeratorRationalComposition(iso_phi_next, phi, iso_phi, h);

                    fmpz_mod_poly_set(iso_psi, iso_psi_next);
                    fmpz_mod_poly_set(iso_phi, iso_phi_next);

                    /* update E2 */

                    fmpz_set(E2->A, A3);
                    fmpz_set(E2->B, B3);
                    fmpz_set(j2, j3);
                    fmpz_set(f1, f2);
                    if (TALKATIVE >= 4)
                    {
                        printf("next");tac(22);
                    }
                }
            }

            fmpz_clear(A3);
            fmpz_clear(B3);
            fmpz_clear(f2);
            fmpz_clear(j3);
            fmpz_mod_poly_clear(psi_f2);
            fmpz_mod_poly_clear(phi);
            fmpz_mod_poly_clear(psi_j2);
            fmpz_mod_poly_clear(h_i);
            fmpz_mod_poly_clear(iso_psi);
            fmpz_mod_poly_clear(iso_phi);
            fmpz_mod_poly_clear(iso_psi_next);
            fmpz_mod_poly_clear(iso_phi_next);
            free_EllipticCurve(E2);
        }
        break;
    }

    for (i = 0; i < r; ++i)
        fmpz_clear(roots + i);

    free(roots);
    fmpz_clear(f1);
    fmpz_clear(j2);
    fmpz_clear(A2);
    fmpz_clear(B2);
    fmpz_clear(tmp);
    fmpz_mod_poly_factor_clear(fac);
    fmpz_mod_poly_clear(psi_f1);
    fmpz_mod_poly_clear(h);

    return tl;
}

AtkinPrime *
findTraces_Atkin(const int r, const PrimeData * prime,
                 const fmpz_t p, AtkinPrime * next)
{
    slong k;
    fmpz_t l_fmpz;
    fq_ctx_t ctx;
    fmpz_mod_poly_t poly;
    int p2l;
    AtkinPrime *candidateTraces = new_AtkinPrime(prime->l, next);

    fq_t g, gpowk;

    fmpz_init_set_ui(l_fmpz, prime->l);
    fmpz_mod_poly_init(poly, l_fmpz);

    fmpz_clear(l_fmpz);

    fmpz_mod_poly_set_coeff_ui(poly, 0, prime->nonSquare);
    fmpz_mod_poly_neg(poly, poly);
    fmpz_mod_poly_set_coeff_ui(poly, 2, 1);

    fq_ctx_init_modulus(ctx, poly, "x");

    fq_init2(g, ctx);
    fq_init2(gpowk, ctx);

    fmpz_poly_set_coeff_ui(g, 0, prime->generator1);
    fmpz_poly_set_coeff_ui(g, 1, prime->generator2);

    fq_pow_ui(g, g, (prime->l * prime->l - 1) / r, ctx);

    p2l = (fmpz_fdiv_ui(p, prime->l) * ((prime->l + 1) / 2)) % prime->l;

    fmpz_poly_one(gpowk);


    for (k = 1; k <= (r >> 1); ++k)
    {
        fq_mul(gpowk, gpowk, g, ctx);

        if (n_gcd(r, k) == 1)
        {
            int z = (p2l * (fmpz_poly_get_coeff_ui(gpowk, 0) + 1)) % prime->l;

            if (z == 0)
                atkin_addCandidate(candidateTraces, 0);

            else if (n_jacobi(z, prime->l) == 1)
            {
                int x = n_sqrtmod(z, prime->l);
                x = (x * 2) % prime->l;
                atkin_addCandidate(candidateTraces, x);
                atkin_addCandidate(candidateTraces, prime->l - x);
            }
        }
    }

    fq_clear(g, ctx);
    fq_clear(gpowk, ctx);
    fq_ctx_clear(ctx);
    fmpz_mod_poly_clear(poly);

    return candidateTraces;
}

/*
 * Returns 0 if the curve of equation y^2 = square_y(x) has an even number of 
 * points over F_p, 1 otherwise
 */

int
order_mod_2(const fmpz_mod_poly_t square_y, const fmpz_t p)
{
    fmpz_mod_poly_t splittingPart, x_p;
    int deg;

    fmpz_mod_poly_init(splittingPart, p);
    fmpz_mod_poly_init(x_p, p);

    fmpz_mod_poly_zero(splittingPart);

    /* splittingPart = x */
    fmpz_mod_poly_set_coeff_ui(splittingPart, 1, 1);

    /* BIG EXPONENTIATION */
    powerXMod(x_p, p, square_y);    /* x_p = x^p */

    /* splittingPart = x^p - x */
    fmpz_mod_poly_sub(splittingPart, x_p, splittingPart);

    /* splittingPart = part of square_2y splitting in F_p */
    fmpz_mod_poly_gcd(splittingPart, splittingPart, square_y);

    deg = fmpz_mod_poly_degree(splittingPart);

    fmpz_mod_poly_clear(splittingPart);
    fmpz_mod_poly_clear(x_p);

    return (deg == 0);
}

/*
 * See paper "Elliptic Curves with the Montgomery-Form and Their Cryptographic
 * Applications" for details (Okeya, Kurumatani, Sakurai)
 * TODO: we could easily extend to check: is montgomery (so 4 is a cofactor),
 * but 8 is not a cofactor of the order... (Theorem 1 in paper)
 */

int
curveIsMontgomery(const fmpz_mod_poly_t square_y, const EllipticCurve * curve)
{
    slong i;
    fmpz *roots = NULL;
    int r = rootsInFp(&roots, square_y, curve->p);

    if (r == 0)
        return 0;
    else {
        int is_montgomery, index, not_divisible_by_8;
        fmpz_t tmp;
        fmpz_t A, B;                /* A and B for the Montgomery form (not Weierstrass!) */

        fmpz_init_set(tmp, roots);  /* tmp = alpha, a root of square_2y */

        fmpz_mul(tmp, tmp, tmp);
        fmpz_mul_ui(tmp, tmp, 3);
        fmpz_add(tmp, tmp, curve->A);   /* tmp = 3*alpha^2 + a */

        is_montgomery = (fmpz_jacobi(tmp, curve->p) == 1);

        if (!is_montgomery)
        {                           /* not Montgomery */

            for (i = 0; i < r; ++i)
                fmpz_clear(roots + i);

            fmpz_clear(tmp);

            return 0;
        }

        fmpz_init(A);
        fmpz_init(B);

        fmpz_invmod(B, tmp, curve->p);
        fmpz_sqrtmod(B, B, curve->p);
        fmpz_mul_ui(A, B, 3);
        fmpz_mul(A, A, roots);
        fmpz_mod(A, A, curve->p);

        index = (fmpz_mod_ui(tmp, curve->p, 4) == 3);
        fmpz_add_ui(tmp, A, 2);
        index = (index << 1) | (fmpz_jacobi(tmp, curve->p) == 1);
        fmpz_sub_ui(tmp, A, 2);
        index = (index << 1) | (fmpz_jacobi(tmp, curve->p) == 1);
        index = (index << 1) | (fmpz_jacobi(B, curve->p) == 1);
        not_divisible_by_8 = montgomery_mod8[index];

        for (i = 0; i < r; ++i)
            fmpz_clear(roots + i);

        fmpz_clear(tmp);
        fmpz_clear(A);
        fmpz_clear(B);

        return (not_divisible_by_8) ? 1 : -1;
    }
}

BSGSList *
multipleScalarMultiplications(BSGSList * bsgs, const ECPoint * P,
                              const ECPoint * PAdd,
                              const EllipticCurve * curve)
{
    fmpz_t tmp;
    fmpz_t tmp2, s;

    int word_size = 8;
    int L = fmpz_sizeinbase(bsgs->ma, 2);

    ECPoint *P4 = new_ECPoint();

    ECPoint ***Ps; /* precomputed data */

    BSGSList *saved_bsgs;

    fmpz_init(tmp);
    fmpz_init(tmp2);
    fmpz_init(s);

    if(TALKATIVE)
        tic(1);

    Ps = precomp_multipleScalarMultiplications(P, word_size, L, curve);

    if (TALKATIVE >= 3)
    {
        printf("precomp");
        tac(1);
    }
    if (TALKATIVE)
        tic(1);

    bsgs = orderBSGSListByRs(bsgs);

    if (TALKATIVE >= 3) {
        printf("sort");
        tac(1);
    }

    saved_bsgs = bsgs;

    scalarMultiplications_withPrecomp(bsgs->P, bsgs->r, word_size, Ps, curve);

    EllipticCurve_add(bsgs->P, bsgs->P, PAdd, curve);

    if (TALKATIVE)
    {
        tic(1);
    }

    while (bsgs->next)
    {

        fmpz_sub(tmp, bsgs->next->r, bsgs->r);

        scalarMultiplications_withPrecomp_optimized(bsgs->next->P, tmp,
                                                    word_size, Ps, curve, P4,
                                                    tmp2, s);

        EllipticCurve_add_memory(bsgs->next->P, bsgs->next->P, bsgs->P, curve,
                                 P4, tmp2, s);
        bsgs = bsgs->next;
    }

    if (TALKATIVE >= 2){
        printf("main");
        tac(1);
    }

    free_precomp_multipleScalarMultiplications(Ps, word_size, L);
    fmpz_clear(tmp);
    fmpz_clear(tmp2);
    fmpz_clear(s);
    free_ECPoint(P4);

    return saved_bsgs;
}

/* M = M_Elkies * M_Atkin */
void
solveBSGS(fmpz_t trace, AtkinPrime * atkinPrimes, const fmpz_t traceElkies,
          const fmpz_t MElkies, const fmpz_t M, const fmpz_t hasseBound,
          const EllipticCurve * curve, flint_rand_t state)
{
    slong k;
    fmpz_t tmp, bsgs2Bound, K;
    fmpz_t MA1, MA2;

    int shiftableBits;
    int not_found, word_size, L;

    AtkinPrime *L1, *L2;

    BSGSList *bsgs1, *bsgs2, *ptr;
    BSGSList *table[HASH_TABLE_SIZE];

    ECPoint *P, *P1, *P2, *deltaP, *minusDeltaP;

    ECPoint ***Ps;

    fmpz_t tmp2, s, rPlusAM2;
    ECPoint *P4, *rP;


    fmpz_init(tmp);
    fmpz_init(K);
    fmpz_init(bsgs2Bound);

    fmpz_cdiv_q(K, hasseBound, M);

    if (TALKATIVE >= 3)
    {
        printf("Initializing BSGS ");
    }
    if (TALKATIVE)
        tic(1);

    if (TALKATIVE >= 3)
    {
        printf("[K = ");
        fmpz_print(K);
        printf("] ");
    }

    evenlySplitAtkinPrimes(&shiftableBits, &L1, &L2, atkinPrimes, K);

    bsgs1 = buildBSGSListOfCandidateTraces(L1);
    ptr = NULL;

    bsgs2 = buildBSGSListOfCandidateTraces(L2);

    if (TALKATIVE >= 3)
    {
        printf("buildBSGS");
        tac(1);
    }

    fmpz_init_set(MA1, bsgs1->ma);
    fmpz_init_set(MA2, bsgs2->ma);

    /*
     * Given r2 mod m2, canditates for actual values of r2 are (r2 + a * m2)
     * such that |r2 + a * m2| <= bsgs2Bound
     */

    fmpz_add_ui(bsgs2Bound, K, 2);
    fmpz_fdiv_q_ui(bsgs2Bound, bsgs2Bound, 2);
    fmpz_mul(bsgs2Bound, bsgs2Bound, MA2);

    computeRInBSGSListCentered(bsgs1, traceElkies, MElkies, bsgs2->ma);
    computeRInBSGSListCentered(bsgs2, traceElkies, MElkies, bsgs1->ma);

    if (TALKATIVE >= 3)
    {
        printf("computeRs");
        tac(1);
    }
    if (TALKATIVE)
    {
        tic(1);
    }

    P = new_ECPoint();
    P1 = new_ECPoint();
    P2 = new_ECPoint();

    EllipticCurve_randomPoint(P, state, curve);

    // if (TALKATIVE >= 5)
    // {
    //     printf("\nP = ");
    //     print_ECPoint(P);
    //     printf("\n");
    // }

    deltaP = new_ECPoint();

    EllipticCurve_scalarMul(deltaP, P, bsgs1->ma, curve);
    EllipticCurve_scalarMul(deltaP, deltaP, bsgs2->ma, curve);
    EllipticCurve_scalarMul(deltaP, deltaP, MElkies, curve);

    minusDeltaP = new_ECPoint_copy(deltaP);

    fmpz_neg(minusDeltaP->y, minusDeltaP->y);

    for (k = 0; k < HASH_TABLE_SIZE; ++k)
        table[k] = NULL;

    if (TALKATIVE >= 3)
    {
        printf("post");
        tac(1);
        printf("\n");
        printf("Baby step: ");
    }

    /* P1 = -[bsgs2->ma * MElkies](xP,yP) */
    EllipticCurve_scalarMul(P1, P, bsgs2->ma, curve);
    EllipticCurve_scalarMul(P1, P1, MElkies, curve);
    fmpz_sub(P1->y, curve->p, P1->y);

    /* tmp = p + 1 - tE */
    fmpz_add_ui(tmp, curve->p, 1);
    fmpz_sub(tmp, tmp, traceElkies);

    /* P2 = [p + 1 - tE]P */
    EllipticCurve_scalarMul(P2, P, tmp, curve);

    bsgs1 = multipleScalarMultiplications(bsgs1, P1, P2, curve);

    if (TALKATIVE)
    {
        tic(1);
    }

    while (bsgs1)
    {
        ptr = bsgs1->next;
        bsgs1->next = table[HASH(bsgs1->P)];
        table[HASH(bsgs1->P)] = bsgs1;
        bsgs1 = ptr;
    }

    if (TALKATIVE >= 3)
    {
        printf("post");
        tac(1);
        printf("\n");
        printf("Giant step: ");
    }

    /* P1 = [bsgs1->ma * MElkies](xP,yP) */
    EllipticCurve_scalarMul(P1, P, MA1, curve);
    EllipticCurve_scalarMul(P1, P1, MElkies, curve);

    not_found = 1;

    if (TALKATIVE)
        tic(1);

    bsgs2 = orderBSGSListByRs(bsgs2);

    if (TALKATIVE >= 3)
    {
        printf("sort");
        tac(1);
    }

    bsgs1 = bsgs2;              /* Save list bsgs2 to free it at the end */
    ptr = bsgs2;

    word_size = 8;
    L = fmpz_sizeinbase(bsgs2->ma, 2);

    if (TALKATIVE)
    {
        tic(1);
    }

    Ps = precomp_multipleScalarMultiplications(P1, word_size, L, curve);

    if (TALKATIVE >= 3)
    {
        printf("precomp");
        tac(1);
    }
    if (TALKATIVE)
    {
        tic(1);
    }

    P4 = new_ECPoint();
    rP = new_ECPoint();
    fmpz_init(tmp2);
    fmpz_init(s);
          
    fmpz_init(rPlusAM2); /* r + a * m2 */

    while (not_found && bsgs2)
    {
        if (bsgs2 != ptr)
        {                       /* If not the first round in the while loop */
            fmpz_sub(tmp, bsgs2->r, ptr->r);
            scalarMultiplications_withPrecomp_optimized(P2, tmp, word_size, Ps,
                                                        curve, P4, tmp2, s);
            EllipticCurve_add_memory(bsgs2->P, ptr->P, P2, curve, P4, tmp2, s);
        }
        else
            scalarMultiplications_withPrecomp(bsgs2->P, bsgs2->r, word_size,
                                              Ps, curve);

        ptr = table[HASH(bsgs2->P)];
        while (not_found && ptr)
        {
            /*
             * If the points are both zero, or both non-zero and agree on
             * x-coordinate
             */

            if ((ptr->P->is_zero == bsgs2->P->is_zero) && (bsgs2->P->is_zero ||
                                                           fmpz_equal(ptr->P->
                                                                      x,
                                                                      bsgs2->
                                                                      P->x)))
            {
                not_found = 0;
                fmpz_mul(trace, bsgs2->r, ptr->ma);

                if (!fmpz_equal(ptr->P->y, bsgs2->P->y))
                    fmpz_neg(trace, trace);

                fmpz_mul(tmp, ptr->r, bsgs2->ma);
                fmpz_add(trace, trace, tmp);
                fmpz_mul(trace, trace, MElkies);
                fmpz_add(trace, trace, traceElkies);
            }

            ptr = ptr->next;
        }

        /*
         * Adding m2's to r
         * rPlusAM2 is now positive (in case it wasn't in the first place)
         */

        fmpz_add(rPlusAM2, bsgs2->r, bsgs2->ma);    
        EllipticCurve_add_memory(rP, bsgs2->P, deltaP, curve, P4, tmp2, s);
        while (not_found && fmpz_cmp(rPlusAM2, bsgs2Bound) <= 0)
        {
            ptr = table[HASH(rP)];
            while (not_found && ptr)
            {
                if ((ptr->P->is_zero == rP->is_zero) && (rP->is_zero ||
                                                         fmpz_equal(ptr->P->x,
                                                                    rP->x)))
                {
                    not_found = 0;

                    fmpz_mul(trace, rPlusAM2, ptr->ma);

                    if (!fmpz_equal(ptr->P->y, rP->y))
                        fmpz_neg(trace, trace);

                    fmpz_mul(tmp, ptr->r, bsgs2->ma);
                    fmpz_add(trace, trace, tmp);
                    fmpz_mul(trace, trace, MElkies);
                    fmpz_add(trace, trace, traceElkies);
                }

                ptr = ptr->next;
            }

            fmpz_add(rPlusAM2, rPlusAM2, bsgs2->ma);
            EllipticCurve_add_memory(rP, rP, deltaP, curve, P4, tmp2, s);
        }

        /*
         * Subtracting m2's to r :
         * rPlusAM2 is now suppose to be negative, but we consider the absolute
         * value
         */

        fmpz_sub(rPlusAM2, bsgs2->ma, bsgs2->r);
        
        EllipticCurve_add_memory(rP, bsgs2->P, minusDeltaP, curve, P4, tmp2, s);
        while (not_found && fmpz_cmp(rPlusAM2, bsgs2Bound) <= 0)
        {
            ptr = table[HASH(rP)];
            while (not_found && ptr)
            {
                if ((ptr->P->is_zero == rP->is_zero) && (rP->is_zero ||
                                                         fmpz_equal(ptr->P->x,
                                                                    rP->x)))
                {
                    not_found = 0;
                    fmpz_mul(trace, rPlusAM2, ptr->ma);

                    if (!fmpz_equal(ptr->P->y, rP->y))
                        fmpz_neg(trace, trace);

                    fmpz_mul(tmp, ptr->r, bsgs2->ma);
                    fmpz_sub(trace, tmp, trace);
                    fmpz_mul(trace, trace, MElkies);
                    fmpz_add(trace, trace, traceElkies);
                }
                ptr = ptr->next;
            }

            fmpz_add(rPlusAM2, rPlusAM2, bsgs2->ma);
            EllipticCurve_add_memory(rP, rP, minusDeltaP, curve, P4, tmp2, s);
        }

        ptr = bsgs2;
        bsgs2 = bsgs2->next;
    }

    if (TALKATIVE >= 3)
    {
        printf("main");
        tac(1);
    }

    if (not_found && (TALKATIVE > 0))
    {
        printf("ERROR: no match was found!\n");
        tac(1);
    }

    bsgs2 = bsgs1;

    free_precomp_multipleScalarMultiplications(Ps, word_size, L);

    fmpz_clear(rPlusAM2);
    fmpz_clear(tmp2);
    fmpz_clear(s);
    free_ECPoint(P4);
    free_ECPoint(rP);

    for (k = 0; k < HASH_TABLE_SIZE; ++k)
        free_BSGSList_deep(table[k]);

    free_BSGSList_deep(bsgs2);

    free_AtkinPrime_deep(L1);
    free_AtkinPrime_deep(L2);

    free_ECPoint(P);
    free_ECPoint(P1);
    free_ECPoint(P2);
    free_ECPoint(deltaP);
    free_ECPoint(minusDeltaP);
    fmpz_clear(K);
    fmpz_clear(bsgs2Bound);
    fmpz_clear(tmp);
    fmpz_clear(MA1);
    fmpz_clear(MA2);
}

/*
 * Adds the first element of A2 to the list A1, returns pointer to new list
 */

AtkinPrime *
addAtkinPrimeToList(AtkinPrime * A1, AtkinPrime * A2)
{
    slong i;
    AtkinPrime *newList = new_AtkinPrime(A2->l, A1);
    for (i = 0; i < A2->N; ++i)
        atkin_addCandidate(newList, A2->tl[i]);

    return newList;
}

/*
 * Goes through the collected data (processed primes and corresponding residues,
 * or candidate residues), and determines wether of not it is sufficient to
 * continue tu the BSGS step
 */

AtkinPrime *
analyseCollectedPrimes(int *enough, fmpz_t finalM, fmpz_t MElkies,
                       fmpz_t hasseBound, AtkinPrime * atkinPrimes,
                       fmpz_t bsgsBound, int buildFinalList)
{
    fmpz_t tmp, M, atkin_candidates_virtual, sqrt_bound_bsgs;
    AtkinPrime *atkin_iterator, *finalListOfAtkinPrimes = NULL;
    int exit_loop;

    fmpz_init(tmp);
    fmpz_init(M);
    fmpz_init(atkin_candidates_virtual);
    fmpz_init(sqrt_bound_bsgs);

    fmpz_sqrt(sqrt_bound_bsgs, bsgsBound);
    fmpz_mul_ui(sqrt_bound_bsgs, sqrt_bound_bsgs, 2);

    fmpz_set(M, MElkies);

    exit_loop = (fmpz_cmp(M, hasseBound) >= 0);
    *enough = exit_loop;

    atkin_iterator = atkinPrimes;
    fmpz_set_ui(atkin_candidates_virtual, 1);

    while (!exit_loop)
    {
        if (atkin_iterator)
        {                       /* If there are more atkin primes available */
            fmpz_mul_ui(tmp, M, atkin_iterator->l);
            if (fmpz_cmp(tmp, hasseBound) <= 0)
            {
                /* atkin_iterator->l is not too big, it is useful! */
                fmpz_set(M, tmp);
                fmpz_mul_ui(atkin_candidates_virtual, atkin_candidates_virtual,
                            atkin_iterator->N);
                if (buildFinalList)
                    finalListOfAtkinPrimes =
                        addAtkinPrimeToList(finalListOfAtkinPrimes,
                                            atkin_iterator);
                if (fmpz_cmp(atkin_candidates_virtual, bsgsBound) > 0)
                {
                    /*
                     * Not going out of the loop yet, unless there are already
                     * too many candidates
                     */
                    exit_loop = 1;
                    *enough = 0;
                }
            }
            else
            {                   /* Getting close to the end, check whats the best move ! */

                /*
                 * Estimate the remaining number of classes if we don't add this
                 * prime
                 */
                fmpz_cdiv_q(tmp, hasseBound, M);

                if (fmpz_cmp_ui(tmp, atkin_iterator->N) <= 0)
                {
                    /*
                     * It is better not to add atkin_iterator->l,
                     * just go to the next prime... maybe we'll be luckier
                     */
                }
                else
                {
                    /*
                     * Getting tricky. Adding atkin_iterator->l is better than
                     * not adding any other atkin prime... but might be better
                     * not to add atkin_iterator->l and to add other primes...
                     *
                     * TODO : possible strategy : remove from atkin_iterator all
                     * the primes with more candidates than atkin_iterator->l,
                     * consider all the remaining primes l such that
                     * l * M > Hasse, and keep only the one with the least
                     * number of candidates.
                     * We are left with that l vs. all the others
                     *
                     * For simplicity, for the moment, just add it...
                     */
                    if (buildFinalList)
                    {
                        finalListOfAtkinPrimes =
                            addAtkinPrimeToList(finalListOfAtkinPrimes,
                                                atkin_iterator);
                    }

                    fmpz_mul_ui(M, M, atkin_iterator->l);
                    fmpz_mul_ui(atkin_candidates_virtual,
                                atkin_candidates_virtual, atkin_iterator->N);
                    exit_loop = 1;
                    *enough =
                        (fmpz_cmp(atkin_candidates_virtual, bsgsBound) <= 0);
                }

            }
            atkin_iterator = atkin_iterator->next;
        }
        else
        {                       /* No Atkin primes left, check if we are close enough to the end */

            /* Estimate the remaining number of classes */
            fmpz_cdiv_q(tmp, hasseBound, M);

            /*
             * TODO: with a proper BSGS for the remaining candidates, the
             * following splitting will not be necessary!
             */
            if (fmpz_cmp(tmp, sqrt_bound_bsgs) < 0)
            {                   /* we miss very little */

                /* Info: we are getting there ! */
                fmpz_mul(atkin_candidates_virtual,
                         atkin_candidates_virtual, tmp);

                *enough = (fmpz_cmp(atkin_candidates_virtual, bsgsBound) <= 0);
            }
            else                /* Missing to much information... look for more */
                *enough = 0;

            exit_loop = 1;      /* We exit the loop since there is no more Atkin prime */
        }
    }

    fmpz_set(finalM, M);

    fmpz_clear(sqrt_bound_bsgs);
    fmpz_clear(tmp);
    fmpz_clear(M);
    fmpz_clear(atkin_candidates_virtual);

    return finalListOfAtkinPrimes;
}

void
computeBSGSBound(fmpz_t bsgsBound, const fmpz_t p)
{
    long logp = fmpz_sizeinbase(p, 2);

    if (logp < 100)
        fmpz_set_ui(bsgsBound, 70000000);

    else if (logp < 255)
    {
        fmpz_set_ui(bsgsBound, logp);
        fmpz_pow_ui(bsgsBound, bsgsBound, 5);
        fmpz_cdiv_q_ui(bsgsBound, bsgsBound, 140);
    }
    else
    {
        fmpz_set_ui(bsgsBound, logp);
        fmpz_mul_ui(bsgsBound, bsgsBound, 35156250);
    }

    /*
     * Good values :
     *   160 - 1'000'000'000 (pari 200'000'000)
     *   192 - (pari 1'000'000'000)
     *   256 - 9'000'000'000 (pari 9'000'000'000)
     *   512 - 18'000'000'000 (pari 52'000'000'000)
     */
}

/*
 * Returns -2 in case of error, -1 in case of early abort, 0 if the cardinality
 * has been computed
 * TODO: algorithm not implemented for j-invariants 0 and 1728!!
 */

int
SEA_AtkinPolynomials(fmpz_t cardinality, const EllipticCurve * curve,
                     PrimeData primes[], flint_rand_t state,
                     AbortStrategy strategy, int twist_secure)
{

    if (EllipticCurve_isSupersingular_proba(curve, state))
    {
       fmpz_add_ui(cardinality, curve->p, 1UL);
       return 0;
    }
    else {
        DivisionPolynomials *divPolys = new_DivisionPolynomials(curve);
        /* TODO: deal with j = 0,1728 */

        /* First, taking care of l = 2 */
        int feedback_mod_2 = 0;
        int enough;
        int i = 0;
        int r;
        fmpz_t hasseBound, M, MElkies, traceElkies, trace, tmp, bsgsBound;

        fmpz_mod_poly_t psi_j1, splittingPart;

        AtkinPrime *atkinPrimes = NULL, *atkin_iterator;
        fmpz_t atkin_candidates_virtual;

        fmpz_t j1;

        ECPoint *P, *P_random;

        if ((strategy == MONTGOMERY_4) || (strategy == MONTGOMERY_8))
        {
            /* If we want a Montgomery curve */
            int feedback_montgomery = curveIsMontgomery(divPolys->square_y, curve);

            if (feedback_montgomery == 0)
            {                       /* if not Montgomery */
                /* Early abort! */
                if (TALKATIVE >= 2)
                    printf("EARLY ABORT: curve is not Montgomery\n");

                free_DivisionPolynomials(divPolys);
                return -1;
            }
            else if (feedback_montgomery == -1)
            {
                feedback_mod_2 = 8; /* 8 is cofactor of the curve */
                if (strategy == MONTGOMERY_4)
                {                   /* cofactor 8 is not allowed */
                    /* Early abort! */
                    if (TALKATIVE >= 2 && (feedback_mod_2 == -1))
                        printf("EARLY ABORT: order divisible by 8\n");

                    free_DivisionPolynomials(divPolys);
                    return -1;
                }
            }
            else
            {
                feedback_mod_2 = 4; /* 4 is cofactor of the curve, but not 8 */
            }
        }
        else
        {                           /* Not looking for Montgomery */
            feedback_mod_2 = order_mod_2(divPolys->square_y, curve->p);
            /* feedback_mod_2 is #E mod 2 */
            if ((feedback_mod_2 == 0) && (strategy == PRIME))
            {
                /* Early abort! */
                if (TALKATIVE >= 2)
                    printf("EARLY ABORT: order divisible by 2\n");

                free_DivisionPolynomials(divPolys);
                return -1;
            }
        }

        /* Compute the j-invariant */
        fmpz_init(j1);
        if (jInvariant(j1, curve) == 0)
        {
            if (TALKATIVE)
                flint_printf("j: UNDEFINED\n");

            free_DivisionPolynomials(divPolys);
            fmpz_clear(j1);
            return -2;
        }

        fmpz_mod(j1, j1, curve->p);
        fmpz_sub_ui(tmp, j1, 1728);

        if (fmpz_is_zero(j1) || fmpz_is_zero(tmp))
            return -2; // TODO: algorithm not implemented for j-invariants 0 and 1728!!

        fmpz_init(traceElkies);
        fmpz_init(tmp);
        fmpz_init(trace);
        fmpz_init_set_ui(M, 1);
        fmpz_init_set_ui(MElkies, 1);
        fmpz_init(bsgsBound);

        computeBSGSBound(bsgsBound, curve->p);

        if (TALKATIVE >= 3)
        {
            printf("BSGS bound: ");
            fmpz_print(bsgsBound);
            printf("\n");
        }

        /* Compute the Hasse bound */
        fmpz_init_set(hasseBound, curve->p);
        fmpz_mul_ui(hasseBound, hasseBound, 16);
        fmpz_sqrt(hasseBound, hasseBound);
        fmpz_add_ui(hasseBound, hasseBound, 1);

        if (TALKATIVE >= 3)
        {
            printf("Hasse Bound: ");
            fmpz_print(hasseBound);
            printf("\n");
        }

        fmpz_mod_poly_init(psi_j1, curve->p);
        fmpz_mod_poly_init(splittingPart, curve->p);

        if (feedback_mod_2 < 2)
        {
            /* We already know the order of the curve is feedback_mod_2 mod 2 */
            fmpz_CRT_ui(traceElkies, traceElkies, MElkies, feedback_mod_2, 2, 1);

            fmpz_mul_ui(MElkies, MElkies, 2);
            fmpz_mul_ui(M, M, 2);
        }
        else
        {                           /* We know the order is 0 mod feedback_mod_2 (which is either 4 or 8) */
            fmpz_CRT_ui(traceElkies, traceElkies, MElkies, 0, feedback_mod_2, 1);
            fmpz_mul_ui(MElkies, MElkies, feedback_mod_2);
            fmpz_mul_ui(M, M, feedback_mod_2);
        }


        /*
         * A virtual atkin prime is an atkin prime, or a prime which has not been
         * processed yet and is considered as Atkin with all the possible traces
         * in its candidate list.
         */

        fmpz_init_set_ui(atkin_candidates_virtual, 1);

        fmpz_set(M, MElkies);

        i = 0;

        enough = 0;

        if (TALKATIVE)
            tic(3);

        while (!enough)
        {

            int earlyAbort = 0;
            int exponent = 1;
            int ln = primes[i].l;

            if (i >= MAX_NBR_SMALL_PRIMES)
            {       /* not enough small primes are available */
                    
                    if (TALKATIVE)
                        printf("ERROR: not enough modular polynomials\n");

                    free_AtkinPrime_deep(atkinPrimes);

                    free_DivisionPolynomials(divPolys);
                    fmpz_mod_poly_clear(splittingPart);
                    fmpz_mod_poly_clear(psi_j1);
                    fmpz_clear(hasseBound);
                    fmpz_clear(M);
                    fmpz_clear(MElkies);
                    fmpz_clear(j1);
                    fmpz_clear(traceElkies);
                    fmpz_clear(trace);
                    fmpz_clear(atkin_candidates_virtual);
                    fmpz_clear(tmp);
                    fmpz_clear(bsgsBound);

                    return -2;
            }
            else if (primes[i].modular == NULL) {
                primes[i].modular = read_ModularPolynomialAsymmetric(primes[i].l, PATH_ATKIN);
            }

            if (TALKATIVE >= 2)
            {
                printf("[%d] \t", primes[i].l);
                fflush(stdout);
            }
            if (TALKATIVE) {
                tic(2);
            }

            evaluate_ModularPolynomialAtY(psi_j1, primes[i].modular, j1);

            prime_type type = splittingModularPolynomial(splittingPart, &r, psi_j1, curve->p);
            long depth_volcano = 0;
            int tl;
            switch (type)
            {

                case ELKIES:
                case SINGLE_ROOT:
                case COMPLETELY_SPLIT:
                /* primes[i].l is an Elkies prime */

                if (TALKATIVE >= 2)
                {
                    printf("Elkies-%ld\t",
                           fmpz_mod_poly_degree(splittingPart));
                    fflush(stdout);
                }

                if (ln == 3)
                {
                    /*
                     * For the fastest possible early abort, treat 3 as a 
                     * special case
                     */
                    if (fmpz_tdiv_ui(curve->p, 3) == 2)
                    {
                        /* then t == #E == 0 mod 3 */
                        if (strategy != NONE)
                        {
                            earlyAbort = 1;
                            break;
                        }
                        tl = 0;

                    }
                    else if (twist_secure)
                    {
                        /*
                         * If p = 1 mod 3, then t = 1 or -1 mod 3, so either
                         * the curve or its twist is 0 mod 3!
                         */
                        earlyAbort = 1;
                        break;
                    }
                    else
                    {

                        /*
                         * in this case, no special treatement for 3: let the general case handle it
                         */

                    }
                }

                switch (type)
                {

                    case ELKIES:
                    {
                        depth_volcano = 0;
                        if (strategy == NONE) 
                        { /* prime powers do not help for early aborts */
                            if (primes[i].l < 14) {
                                fmpz_set_ui(tmp,(long)fmpz_dlog(curve->p));
                                double a = fmpz_dlog(tmp);
                                fmpz_set_ui(tmp,primes[i].l);
                                exponent = (long)(a/fmpz_dlog(tmp));

                                if (exponent < 1)
                                    exponent = 1;

                                for (int j = 1; j < exponent; ++j)
                                    ln *= primes[i].l;
                            }
                        }
                        tl = findTrace_Elkies(splittingPart, primes[i].modular,
                                                  exponent, divPolys, psi_j1, j1,
                                                  curve, twist_secure);
                    }
                    break;
                    /* either one root or l+1 roots: l is ramified in Z[frob] */
                    case SINGLE_ROOT:
                    {
                        fmpz_t j2;
                        fmpz_init(j2);
                        fmpz_mod_poly_t h;
                        fmpz_mod_poly_init(h, curve->p);

                        if (TALKATIVE >= 3)
                            tic(15);

                        findOneJNeighbor(j2, h, splittingPart, primes[i].modular, psi_j1, j1, curve);

                        if (countJNeighbors(primes[i].modular, j2, curve) == 1)
                            depth_volcano = 0; /* the depth of the volcano is 0 */
                        else
                            depth_volcano = 1; /* the depth of the volcano is AT LEAST 1 */

                        if (TALKATIVE >= 3){
                            printf("depth"); tac(15);
                        }

                        tl = findTrace_Elkies_ramified(primes[i].l,divPolys,h,curve,twist_secure);

                        fmpz_mod_poly_clear(h);
                        fmpz_clear(j2);
                    }
                    break;

                    case COMPLETELY_SPLIT:
                    {                   
                        depth_volcano = 1; /* the depth of the volcano is AT LEAST 1 */

                        fmpz_t j2;
                        fmpz_init(j2);
                        fmpz_mod_poly_t h;
                        fmpz_mod_poly_init(h, curve->p);
                        
                        findOneJNeighbor(j2, h, splittingPart, primes[i].modular, psi_j1, j1, curve);

                        tl = findTrace_Elkies_ramified(primes[i].l,divPolys,h,curve,twist_secure);

                        fmpz_mod_poly_clear(h);
                        fmpz_clear(j2);
                    }
                    break;
                    default: /* never reached */
                    break;
                }
     
                /*
                 * findTrace_Elkies returns -1 in case of early abort, but does not
                 * necessarily detect every occasions to abort
                 */

                if (tl == -1)
                {
                    earlyAbort = 1;
                    break;
                }
                if (strategy != NONE)
                {
                    int pl = fmpz_mod_ui(tmp, curve->p, primes[i].l);
                    if ((pl + 1 - tl) % primes[i].l == 0)
                    {
                        earlyAbort = 1;
                        break;
                    }
                    if (twist_secure && ((pl + 1 + tl) % primes[i].l == 0))
                    {   /* check the twist */
                        earlyAbort = 1;
                        break;
                    }
                }

                if (depth_volcano > 0)
                {
                    if (TALKATIVE >= 3)
                        printf("VOLCANO DETECTED\t");
                    int l = primes[i].l;
                    int n_roots = 0;
                    mp_limb_t *sqrt = NULL;

                    ln = 1;
                    for (int i = 0; i < depth_volcano; ++i)
                    {
                        ln *= l*l;
                    }

                    int p_ln = fmpz_mod_ui(tmp, curve->p, ln);
                    n_roots = n_sqrtmod_primepow(&sqrt, 4*p_ln, l, 2*depth_volcano);

                    for (int i = 0; i < n_roots; ++i)
                    {
                        if ((tl % l) == (*(sqrt + i) % l))
                        {
                            tl = *(sqrt + i);
                            break;
                        }
                    }

                    free(sqrt);
                }


                fmpz_CRT_ui(traceElkies, traceElkies, MElkies, tl, ln, 1);
                fmpz_mul_ui(MElkies, MElkies, ln);
                if (TALKATIVE >= 2)
                    printf("t mod %d = %d\t", ln, tl);


                break;

                case ATKIN:
                {                   /* primes[i].l is an Atkin prime */
                    /*
                     * TODO: For small Atkin primes, and big field, use the original Schoof algorithm
                     */

                    atkinPrimes = findTraces_Atkin(r, primes + i, curve->p,
                                                   atkinPrimes);

                    if (atkinPrimes == NULL)
                    {
                        if (TALKATIVE) printf("ERROR: no traces for Atkin prime\n");
                        //earlyAbort = 1;
                        //break;

                        free_AtkinPrime_deep(atkinPrimes);
                        free_DivisionPolynomials(divPolys);
                        fmpz_mod_poly_clear(splittingPart);
                        fmpz_mod_poly_clear(psi_j1);
                        fmpz_clear(hasseBound);
                        fmpz_clear(M);
                        fmpz_clear(MElkies);
                        fmpz_clear(j1);
                        fmpz_clear(traceElkies);
                        fmpz_clear(trace);
                        fmpz_clear(atkin_candidates_virtual);
                        fmpz_clear(tmp);
                        fmpz_clear(bsgsBound);

                        return -2;
                    }

                    if (TALKATIVE >= 2)
                    {
                        printf("%d candidates\t", atkinPrimes->N);
                        printf("score %.2f\t",
                               (float) (atkinPrimes->l) / (float) atkinPrimes->N);
                    }

                    /*
                     * Can the following ever happen?
                     * if (atkinPrimes->N == 2)
                     * { 
                     *     Only two candidates, can detect early abort.
                     *     int pl = fmpz_mod_ui(tmp, curve->p, primes[i].l);
                     *
                     *     If l divides pl + 1 - tl, then either the curve or its 
                     *     twist has order divisible by l.
                     *     Does this ever happen ?
                     *
                     *      if (((pl + 1 - atkinPrimes->tl[0]) % primes[i].l == 0)
                     *          || ((pl + 1 + atkinPrimes->tl[0]) % primes[i].l == 0))
                     *          printf("ATKIN EARLY ABORT\n");
                     * }
                     */


                    if (atkinPrimes->N == 1)
                    {               /* Only one candidate, treat it as Elkies */

                        fmpz_CRT_ui(traceElkies, traceElkies, MElkies,
                                    atkinPrimes->tl[0], primes[i].l, 1);
                        fmpz_mul_ui(MElkies, MElkies, primes[i].l);

                        atkin_iterator = atkinPrimes->next;
                        free_AtkinPrime(atkinPrimes);
                        atkinPrimes = atkin_iterator;

                    }

                }
                    break;

                case DROP:
                    /* Happens for Atkin primes which are too large to be interesting */
                    break;
            }

            if (TALKATIVE >= 2)
            {
                tac(2);
                printf("\n");
            }

            if (earlyAbort)
            {
                /*
                 * Early abort! Small divisor found for the curve's order or its
                 * twist's
                 */

                if (TALKATIVE >= 2)
                    printf("EARLY ABORT: order or twist divisible by %d\n",
                           primes[i].l);

                free_AtkinPrime_deep(atkinPrimes);
                free_DivisionPolynomials(divPolys);
                fmpz_mod_poly_clear(splittingPart);
                fmpz_mod_poly_clear(psi_j1);
                fmpz_clear(hasseBound);
                fmpz_clear(M);
                fmpz_clear(MElkies);
                fmpz_clear(j1);
                fmpz_clear(traceElkies);
                fmpz_clear(trace);
                fmpz_clear(atkin_candidates_virtual);
                fmpz_clear(tmp);
                fmpz_clear(bsgsBound);

                return -1;
            }

            fmpz_set(M, MElkies);

            atkinPrimes = insertAtkinPrimes(atkinPrimes);
            analyseCollectedPrimes(&enough, M, MElkies, hasseBound, atkinPrimes,
                                   bsgsBound, 0);

            ++i;

        }

        if (TALKATIVE >= 2)
        {
            printf("Total collection time ");
            tac(3);
            printf("\n");
        }

        /* Find the atkin primes which will be used in BSGS */
        atkin_iterator = analyseCollectedPrimes(&enough, M, MElkies, hasseBound,
                                                atkinPrimes, bsgsBound, 1);

        free_AtkinPrime_deep(atkinPrimes);

        solveBSGS(trace, atkin_iterator, traceElkies, MElkies, M, hasseBound,
                  curve, state);

        fmpz_add_ui(cardinality, curve->p, 1);
        fmpz_sub(cardinality, cardinality, trace);

        /* Sanity check */
        P = new_ECPoint();
        P_random = new_ECPoint();
        
        EllipticCurve_randomPoint(P_random, state, curve);

        tic(11);
        EllipticCurve_scalarMul(P, P_random, cardinality, curve);

        if (TALKATIVE >= 3)
        {
            printf("\n");
            printf("sanity");
            tac(11);
        }
        if (TALKATIVE >= 2) {
            if (P->is_zero)
                printf("Correct!\n");
        }

        if (!P->is_zero)
            printf("ERROR! Wrong cardinality!\n");

        if (TALKATIVE >= 2)
        {
            printf("N = ");
            fmpz_print(cardinality);
            printf("\n");
        }

        free_ECPoint(P);
        free_ECPoint(P_random);

        free_DivisionPolynomials(divPolys);
        fmpz_mod_poly_clear(splittingPart);
        fmpz_mod_poly_clear(psi_j1);
        fmpz_clear(hasseBound);
        fmpz_clear(M);
        fmpz_clear(MElkies);
        fmpz_clear(j1);
        fmpz_clear(traceElkies);
        fmpz_clear(trace);
        fmpz_clear(atkin_candidates_virtual);
        fmpz_clear(tmp);
        fmpz_clear(bsgsBound);

        return 0;
    }
}
