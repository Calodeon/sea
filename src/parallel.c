// gcc main.c -o gen -Wall -I/usr/local/include/flint/ -I/usr/local/include -lgmp -lflint -L/usr/local/lib/

#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "SEA/EllipticCurvePointCounting.h"
#include "SEA/toolbox.h"

#include "SEA/params.h"

#define NTHREADS 4
#define SECURITY_BITS 256

static pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;
static int twist_counter = 0;

static clock_t timer;

typedef struct
{
    int id;
    fmpz_t p;
    flint_rand_t state;
    int *isInterrupted;
} ThreadData;

void *
testing(void *x)
{
    ThreadData *data = ((ThreadData *) x);
    fmpz_t A, B, n;


    flint_randinit(data->state);
    gmp_randinit_default(data->state->gmp_state);
    data->state->gmp_init = 1;
    gmp_randseed_ui(data->state->gmp_state, 12345 + 1 + data->id);


    fmpz_init(n);
    int i = data->id;

    while (!(*(data->isInterrupted)))
    {
        i += NTHREADS;
        fmpz_randm(A, data->state, data->p);
        fmpz_randm(B, data->state, data->p);

        if (SEA_simple(n, data->p, A, B) == 0)
        {

            if (fmpz_is_probabprime(n))
            {

                fmpz_neg(n, n);
                fmpz_sub_ui(n, n, 2);
                fmpz_add(n, n, data->p);
                fmpz_add(n, n, data->p);

                if (fmpz_is_probabprime(n))
                {

                    pthread_mutex_lock(&file_mutex);
                    printf("Twist-secure curve %d: ", i);
                    flint_printf("A = ");
                    fmpz_print(A);
                    flint_printf(", ");
                    flint_printf("B = ");
                    fmpz_print(B);
                    flint_printf("\n");


                    FILE *twist_secure_list = fopen("twist_secure.txt", "a");

                    if (twist_secure_list == NULL)
                    {
                        printf("ERROR opening file twist_secure.txt\n");
                    }
                    else
                    {
                        // format is "[index, time, p, A, B, N]\n"
                        twist_counter++;

                        fprintf(twist_secure_list, "[%d, ", twist_counter);

                        fprintf(twist_secure_list, "%d, ",
                                (int) ((float) (clock() - timer) /
                                       CLOCKS_PER_SEC));

                        char *str = fmpz_get_str(NULL, 16, data->p);
                        fprintf(twist_secure_list, "%s, ", str);
                        free(str);

                        str = fmpz_get_str(NULL, 16, A);
                        fprintf(twist_secure_list, "%s, ", str);
                        free(str);

                        str = fmpz_get_str(NULL, 16, B);
                        fprintf(twist_secure_list, "%s, ", str);
                        free(str);

                        str = fmpz_get_str(NULL, 16, n);
                        fprintf(twist_secure_list, "%s]\n", str);
                        free(str);

                        fclose(twist_secure_list);
                    }

                    pthread_mutex_unlock(&file_mutex);

                }
                else
                {
                    printf("Prime curve #%d, time %ds: ", i,
                           (int) ((float) (clock() - timer) / CLOCKS_PER_SEC));
                    flint_printf("A = ");
                    fmpz_print(A);
                    flint_printf(", ");
                    flint_printf("B = ");
                    fmpz_print(B);
                    flint_printf("\n");
                }

            }
        }

    }
    flint_randclear(data->state);
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(n);
    return NULL;
}

int
main(int argc, char *argv[])
{

    slong i;
    fmpz_t p, n, A, B;

    fmpz_init(p);
    fmpz_init(n);
    fmpz_init(A);
    fmpz_init(B);

    fmpz_set_ui(n, 2);
    fmpz_pow_ui(n, n, SECURITY_BITS);

    nextPrime(p, n);

    fmpz_set_ui(n, 2);
    fmpz_pow_ui(n, n, 100);

    //fmpz_set_str(p, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffc31", 16);    // 512
    fmpz_set_str(p, "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16); // 256
    //fmpz_set_str(p, "ffffffffffffffffffffffffffffffffffffffb5", 16); // 160
    //fmpz_set_str(p, "ffffffffffffffffffffffffffff89", 16); // 120
    //fmpz_set_str(p, "ffffffffffffffffffbf", 16); // 80


    flint_printf("Prime: ");
    fmpz_print(p);
    flint_printf("\n");
    printf("Initialization...\n");


    SEA_init();
    timer = clock();

    printf("Research...\n");

    pthread_t threads[NTHREADS];
    int isInterrupted = 0;

    int feedback;

    ThreadData data[NTHREADS];

    for (i = 0; i < NTHREADS; ++i)
    {
        fmpz_init_set(data[i].p, p);
        data[i].id = i;
        data[i].isInterrupted = &isInterrupted;
        //testing(&(data[i]));
        feedback =
            pthread_create(&threads[i], NULL, testing, (void *) &data[i]);
        if (feedback != 0)
        {
            printf("ERROR pthread creation failed with error code %i\n",
                   feedback);
            return feedback;
        }
    }

    for (i = 0; i < NTHREADS; ++i)
    {
        feedback = pthread_join(threads[i], NULL);
        if (feedback != 0)
        {
            printf("ERROR pthread join failed with error code %i\n", feedback);
            return feedback;
        }
    }

    SEA_clear();

    for (i = 0; i < NTHREADS; ++i)
    {
        fmpz_clear(data[i].p);
    }

    fmpz_clear(p);
    fmpz_clear(n);
    fmpz_clear(A);
    fmpz_clear(B);


    return 0;
}
