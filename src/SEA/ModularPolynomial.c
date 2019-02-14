/* reading modular polynomials */

#include "ModularPolynomial.h"
#include <stdio.h>
#include <string.h>

ModularPolynomial *
read_ModularPolynomialSymmetric(const int l, char path_model[])
{
    slong i;
    fmpz_t coef, tmp;
    char path[1024];
    FILE *file;
    int a = 0, b = 0, status, success;

    ModularPolynomial *poly = malloc(sizeof(ModularPolynomial));

    poly->l = l;
    poly->maxDeg = l + 1;

    poly->coef = malloc(sizeof(fmpz_poly_t) * (poly->maxDeg + 1));
    for (i = 0; i < poly->maxDeg + 1; ++i)
        fmpz_poly_init(poly->coef[i]);

    sprintf(path, path_model, l);

    file = fopen(path, "r");
    if (!file)
    {
        printf("error opening file \"%s\"\n", path);
        free_ModularPolynomial(poly);
        return 0;
    }

    fmpz_init(coef);
    fmpz_init(tmp);

    success = fscanf(file, "[%d,%d]", &a, &b);
    if (success < 2)
    {
        printf("A\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free_ModularPolynomial(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    status = fmpz_fread(file, coef);
    success = fscanf(file, "\n");
    if (success < 0)
    {
        printf("B\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free_ModularPolynomial(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    while (status > 0)
    {

        fmpz_poly_get_coeff_fmpz(tmp, poly->coef[a], b);
        fmpz_add(coef, coef, tmp);

        fmpz_poly_set_coeff_fmpz(poly->coef[a], b, coef);

        if (a != b)
        {
            fmpz_poly_get_coeff_fmpz(tmp, poly->coef[b], a);
            fmpz_add(coef, coef, tmp);

            fmpz_poly_set_coeff_fmpz(poly->coef[b], a, coef);
        }

        success = fscanf(file, "[%d,%d]", &a, &b);
        if (success < 2)
        {
            printf("C\n");
            printf("error scanning file \"%s\" : success : %i\n", path, success);

            free_ModularPolynomial(poly);
            fmpz_clear(coef);
            fmpz_clear(tmp);

            fclose(file);

            return 0;
        }

        status = fmpz_fread(file, coef);
        success = fscanf(file, "\n");
        if (success < 0)
        {
            printf("D\n");
            printf("error scanning file \"%s\" : success : %i\n", path, success);

            free_ModularPolynomial(poly);
            fmpz_clear(coef);
            fmpz_clear(tmp);

            fclose(file);

            return 0;
        }
    }

    fmpz_clear(coef);
    fmpz_clear(tmp);

    fclose(file);

    return poly;
}

ModularPolynomial *
read_ModularPolynomialAsymmetric(const int l, char path_model[])
{
    slong i;
    fmpz_t coef, tmp;
    char path[1024];
    FILE *file;
    int a = 0, b = 0, status, success;
    ModularPolynomial *poly = malloc(sizeof(ModularPolynomial));

    if(TALKATIVE >= 4) { printf("> Loading modular polynomial %d\n", l); }

    poly->l = l;

    sprintf(path, path_model, l);

    file = fopen(path, "r");
    if (!file)
    {
        printf("error opening file \"%s\"\n", path);
        free(poly);
        return 0;
    }

    fmpz_init(coef);
    fmpz_init(tmp);

    poly->maxDeg = 0;
    success = fscanf(file, "[%d,%d]", &b, &(poly->maxDeg));
    if (success < 2)
    {
        printf("E\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    status = fmpz_fread(file, coef);
    success = fscanf(file, "\n");
    if (success < 0)
    {
        printf("F\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    while (status > 0 && success >= 0)
    {
        success = fscanf(file, "[%d,%d]", &b, &a);

        if (success >= 0) {
            if (a > poly->maxDeg)
                poly->maxDeg = a;

            status = fmpz_fread(file, coef);
            success = fscanf(file, "\n");
        }
    }

    poly->coef = malloc(sizeof(fmpz_mod_poly_t) * (poly->maxDeg + 1));
    for (i = 0; i < poly->maxDeg + 1; ++i)
        fmpz_poly_init(poly->coef[i]);

    rewind(file);

    success = fscanf(file, "[%d,%d]", &b, &a);
    if (success < 2)
    {
        printf("I\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free_ModularPolynomial(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    status = fmpz_fread(file, coef);
    success = fscanf(file, "\n");
    if (success < 0)
    {
        printf("J\n");
        printf("error scanning file \"%s\" : success : %i\n", path, success);

        free_ModularPolynomial(poly);
        fmpz_clear(coef);
        fmpz_clear(tmp);

        fclose(file);

        return 0;
    }

    while (status > 0 && success >= 0)
    {

        fmpz_poly_get_coeff_fmpz(tmp, poly->coef[a], b);
        fmpz_add(coef, coef, tmp);

        fmpz_poly_set_coeff_fmpz(poly->coef[a], b, coef);

        success = fscanf(file, "[%d,%d]", &b, &a);
        if (success >= 2)
        {
            status = fmpz_fread(file, coef);
            success = fscanf(file, "\n");
        }
    }

    fmpz_clear(coef);
    fmpz_clear(tmp);

    fclose(file);

    return poly;
}

/* evaluate at X = f */
void
evaluate_ModularPolynomialAtX(fmpz_mod_poly_t result,
                              const ModularPolynomial * poly, const fmpz_t f)
{
    slong i;
    fmpz_t coef;
    fmpz_mod_poly_t tmp;

    fmpz_init(coef);
    fmpz_mod_poly_init(tmp, fmpz_mod_poly_modulus(result));

    fmpz_mod_poly_zero(result);

    for (i = 0; i < poly->maxDeg + 1; ++i)
    {
        fmpz_mod_poly_set_fmpz_poly(tmp, poly->coef[i]);
        fmpz_mod_poly_evaluate_fmpz(coef, tmp, f);
        fmpz_mod_poly_set_coeff_fmpz(result, i, coef);
    }

    fmpz_clear(coef);
    fmpz_mod_poly_clear(tmp);
}

/* evaluate at Y = f */
void
evaluate_ModularPolynomialAtY(fmpz_mod_poly_t result,
                              const ModularPolynomial * poly, const fmpz_t f)
{
    slong i;
    fmpz_t coef;
    fmpz_mod_poly_t tmp;

    fmpz_init(coef);
    fmpz_mod_poly_init(tmp, fmpz_mod_poly_modulus(result));

    fmpz_mod_poly_zero(result);

    for (i = 0; i < poly->maxDeg + 1; ++i)
    {
        fmpz_powm_ui(coef, f, i, fmpz_mod_poly_modulus(result));
        fmpz_mod_poly_set_fmpz_poly(tmp, poly->coef[i]);
        fmpz_mod_poly_scalar_mul_fmpz(tmp, tmp, coef);
        fmpz_mod_poly_add(result, result, tmp);
    }

    fmpz_clear(coef);
    fmpz_mod_poly_clear(tmp);
}

/* evaluate the X-derivative at X = f */
void
evaluate_ModularPolynomialDXAtX(fmpz_mod_poly_t result,
                                const ModularPolynomial * poly, const fmpz_t f)
{
    slong i;
    fmpz_t coef;
    fmpz_mod_poly_t derivative;

    fmpz_init(coef);
    fmpz_mod_poly_init(derivative, fmpz_mod_poly_modulus(result));

    fmpz_mod_poly_zero(result);

    for (i = 0; i < poly->maxDeg + 1; ++i)
    {
        fmpz_mod_poly_set_fmpz_poly(derivative, poly->coef[i]);
        fmpz_mod_poly_derivative(derivative, derivative);
        fmpz_mod_poly_evaluate_fmpz(coef, derivative, f);
        fmpz_mod_poly_set_coeff_fmpz(result, i, coef);
    }

    fmpz_mod_poly_clear(derivative);
    fmpz_clear(coef);
}

/* computes F'(f) (j = F(f), where f is the weber f function) */
void
weberToJDerivative(fmpz_t result, const fmpz_t f, const fmpz_t p)
{
    /* F1 = F'(f1) (j = F(f), where f is the weber f function) */
    fmpz_t tmp;

    fmpz_init(tmp);

    fmpz_powm_ui(tmp, f, 24, p);
    fmpz_sub_ui(result, tmp, 16);
    fmpz_mul(result, result, result);
    fmpz_mod(result, result, p);
    fmpz_mul_ui(result, result, 48);
    fmpz_add_ui(tmp, tmp, 8);
    fmpz_mul(result, result, tmp);
    fmpz_mod(result, result, p);
    fmpz_sub_ui(tmp, tmp, 8);
    fmpz_mul(tmp, tmp, f);
    fmpz_invmod(tmp, tmp, p);
    fmpz_mul(result, result, tmp);
    fmpz_mod(result, result, p);

    fmpz_clear(tmp);
}

/* computes F'(f) and F''(f) (j = F(f), where f is the weber f function) */
void
weberToJDerivativeDouble(fmpz_t FP, fmpz_t FPP, const fmpz_t f, const fmpz_t p)
{
    fmpz_t tmp, f24;

    fmpz_init(tmp);
    fmpz_init(f24);

    /* computing FP = F'(f) */
    fmpz_powm_ui(f24, f, 24, p);
    fmpz_sub_ui(FP, f24, 16);
    fmpz_mul(FP, FP, FP);
    fmpz_mod(FP, FP, p);
    fmpz_mul_ui(FP, FP, 48);
    fmpz_add_ui(tmp, f24, 8);
    fmpz_mul(FP, FP, tmp);
    fmpz_mod(FP, FP, p);
    fmpz_mul(tmp, f24, f);
    fmpz_invmod(tmp, tmp, p);
    fmpz_mul(FP, FP, tmp);
    fmpz_mod(FP, FP, p);

    /* computing FPP = F''(f) */
    fmpz_mul(tmp, f24, f24);
    fmpz_mod(tmp, tmp, p);
    fmpz_mul_ui(FPP, tmp, 552);
    fmpz_add_ui(FPP, FPP, 51200);
    fmpz_mul(tmp, tmp, f24);
    fmpz_mod(tmp, tmp, p);
    fmpz_mul_ui(tmp, tmp, 47);
    fmpz_sub(FPP, tmp, FPP);
    fmpz_mul_ui(FPP, FPP, 48);
    fmpz_mul(tmp, f24, f);
    fmpz_mod(tmp, tmp, p);
    fmpz_mul(tmp, tmp, f);
    fmpz_mod(tmp, tmp, p);
    fmpz_invmod(tmp, tmp, p);
    fmpz_mul(FPP, FPP, tmp);
    fmpz_mod(FPP, FPP, p);

    fmpz_clear(tmp);
    fmpz_clear(f24);
}

void
print_ModularPolynomial(const ModularPolynomial * poly)
{
    slong i;

    for (i = 0; i < poly->maxDeg; ++i)
    {
        flint_printf("Y^%wd(", i);
        fmpz_poly_print_pretty(poly->coef[i], "X");
        flint_printf(") + ");
    }
    printf("Y^%d(", poly->maxDeg);
    fmpz_poly_print_pretty(poly->coef[poly->maxDeg], "X");
    printf(")\n");
}

void
free_ModularPolynomial(ModularPolynomial * poly)
{
    slong i;

    if (poly)
    {
        for (i = 0; i < poly->maxDeg + 1; ++i)
            fmpz_poly_clear(poly->coef[i]);

        free(poly->coef);
        free(poly);
    }
}
