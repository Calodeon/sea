
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "SEA/EllipticCurvePointCounting.h"
#include "SEA/toolbox.h"

#include "SEA/params.h"

#define SECURITY_BITS 256
#define PATH_AS "./build/data/A256"
#define PATH_BS "./build/data/B256"


int
main(int _, char **v)
{

    fmpz_t p, p2, n, A, B;
    slong i;

    fmpz_init(p);
    fmpz_init(p2);
    fmpz_init(n);
    fmpz_init(A);
    fmpz_init(B);


    //fmpz_set_str(p2, "fa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bf", 16); // example with a volcano at ell = 3
    //fmpz_set_str(p2, "79bfa2079bfa20779bfa20779bf20779bfa20779bfa20779bfa20779bf207", 16);
    fmpz_set_str(p2, "79bfa2079bfa20779bfa20779bf20779bfa2", 16);
    //fmpz_set_str(p2, "fa20779b", 16); // a very small example

    nextPrime(p, p2);

    flint_printf("p = ");
    fmpz_print(p);
    flint_printf("\n");
    printf("Initialization...\n");

    SEA_init();

    fmpz As[250], Bs[250];
    FILE *fileA = fopen(PATH_AS, "r"), *fileB = fopen(PATH_BS, "r");
    char *line = malloc(sizeof(char) * 1000);
    size_t len = 0;

    if (!fileA || !fileB)
    {
        printf("error opening files \"A256\" and \"B256\"\n");
        return 0;
    }

    for (i = 0; i < 250; ++i)
    {

        if (getline(&line, &len, fileA) == -1)
        {
            printf("Not enough entries in \"A256\"\n");
            return 0;
        }
        fmpz_init(As + i);
        fmpz_set_str(As + i, line, 16);
        if (getline(&line, &len, fileB) == -1)
        {
            printf("Not enough entries in \"B256\"\n");
            return 0;
        }
        fmpz_init(Bs + i);
        fmpz_set_str(Bs + i, line, 16);

    }

    fclose(fileA);
    fclose(fileB);
    if (line)
        free(line);


    printf("Computing...\n");

    int a = 1;
    int b = 61;

    for (a = 2; a <= 10; ++a)
    {
        for (b = 1; b <= 10; ++b)
        {

            fmpz_set(A, As + a - 1);
            fmpz_set(B, Bs + b - 1);

            tic(0);

            //if (SEA_simple(n, p, A, B) == 0)
            if (SEA(n, p, A, B, NONE, 0) == 0)
            {

                printf("N = ");
                fmpz_print(n);

                printf("\n");


                flint_printf("A = ");
                fmpz_print(A);
                flint_printf("\n");
                flint_printf("B = ");
                fmpz_print(B);
                flint_printf("\n");



                if (fmpz_is_probabprime(n))
                {
                    printf("\n(%d,\t%d)\t", a, b);
                    fflush(stdout);
                    printf("N = ");
                    fmpz_print(n);

                    fmpz_neg(n, n);
                    fmpz_add_ui(n, n, 2);
                    fmpz_add(n, n, p);
                    fmpz_add(n, n, p);

                    if (fmpz_is_probabprime(n))
                    {
                        printf("\tTwist-secure curve");
                    }
                    else
                    {
                        printf("\tPrime curve");
                    }
                    printf("\n");

                }

            }
            else
            {
                /* printf("bad\n"); */
            }

            printf("SEA ");
            tac(0);
            printf("\n");

    
        }
    }

    SEA_clear();

    fmpz_clear(p);
    fmpz_clear(p2);
    fmpz_clear(n);
    fmpz_clear(A);
    fmpz_clear(B);

    for (i = 0; i < 250; ++i)
    {
        fmpz_clear(As + i);
        fmpz_clear(Bs + i);
    }


    return 0;
}
