
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

    fmpz_t p, p1, p2, n, A, B;
    slong i;

    fmpz_init(p);
    fmpz_init(p1);
    fmpz_init(p2);
    fmpz_init(n);
    fmpz_init(A);
    fmpz_init(B);

/*
    // p1 = 2^256-2^224+2^192+2^96-1
    //fmpz_set_str(p1, "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16);
    //fmpz_set_str(p2, "8000000b", 16); // 30
    //fmpz_set_str(p2, "ffffffffffffffffffbf", 16); // 80
    //fmpz_set_str(p2, "ffffffffffffffffffffffffffff89", 16); // 120
    //fmpz_set_str(p2, "ffffffffffffffffffffffffffffffffffffffb5", 16);
    //fmpz_set_str(p2, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16);
    //fmpz_set_str(p2, "143294146727694541181438963548158737568808714913643145551794573457318826785209", 10);
    //fmpz_set_str(p2, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec3", 16);
    
    
    fmpz_set_str(p2, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc31", 16);
    */

    //fmpz_set_str(p2, "222481008190394217244519595537473922147738241476467368881208896887570023433541", 10);
    //fmpz_set_str(p2, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec3", 16);

    fmpz_set_str(p2, "fa20779b", 16); // 32
    fmpz_set_str(p2, "fa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bf", 16); // example with a volcano at ell = 3
    //fmpz_set_str(p2, "fa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa20779bfa2077", 16); // a good example
    //fmpz_set_str(p2, "cfa20779bfa2079bfa20779bfa20779bf20779bfa20779bfa20779bfa20779bfa20779bfa2077", 16); // a example of totally split polynomial
    fmpz_set_str(p2, "79bfa2079bfa20779bfa20779bf20779bfa20779bfa20779bfa20779bfa20779bfa2077", 16); // lots of small Elkies-2 and a volcano
    fmpz_set_str(p2, "79bfa2079bfa20779bfa20779bf20779bfa20779bfa20779bfa20779bfa207", 16);// BSGS is too slow!

    fmpz_set_str(p2, "79bfa2079bfa20779bfa20779bf20779bfa20779bfa20779bfa20779bf207", 16);

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

    for (a = 2; a <= 2; ++a)
    {
        for (b = 1; b <= 1; ++b)
        {

    /*

    // (1,  30)     for the random prime

    // (1,  191)    N = 1461501637330902918203684802757851899579874130883

    // (1,  61)     N = 39402006196394479212279040100143613805079739270465446667954644381951029102592974749194850288263425714245031242752969

    // (2,173)      N = 13407807929942597099574024998205846127479365820592393377723561443721764030073716727470690379682046336162977123878693911073393884217346248783955285983522527

    // (75, 139)    Interesting, large K

    // (45, 32)     Twist-secure

    // (1,  21)     N = 115792089237316195423570985008687907853678737181469263003689564268589498588517

    // (1,  155)    N = 115792089237316195423570985008687907852790345417262834154816357227779229390721

    // (1,  199)    N = 115792089237316195423570985008687907852873215582086565804271915348399548582469

    // (36, 126)    N = 115792089237316195423570985008687907852922043384143762988465437130772777567057 IS WRONG !!
    */

    fmpz_set(A, As + a - 1);
    fmpz_set(B, Bs + b - 1);


    fmpz_set_str(A, "2", 16);
    fmpz_set_str(B, "0", 16);


    /*
    //fmpz_set_str(A, "12262006897685374802459563478209403269601397690715000682733368792197394659564763981711177157072342628291477231626062217776282472947907019062743220907697790", 10);

    //fmpz_set_str(B, "6234153715992723126717925414482162937811959033417435472435973370761990827444042321082187458434297662217189326884751629526112335681076712624996207473588958", 10);
    */

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
    if (TALKATIVE)
    {
        printf("SEA ");
        tac(0);
        printf("\n");
    }

    
        }
    }

    SEA_clear();

    fmpz_clear(p);
    fmpz_clear(p1);
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
