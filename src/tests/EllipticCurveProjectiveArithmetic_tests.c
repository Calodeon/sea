#include <stdio.h>
#include <stdlib.h>
#include "fmpz.h"
#include "../tests.h"
#include "../SEA/EllipticCurveArithmetic.h"
#include "../SEA/EllipticCurveProjectiveArithmetic.h"
#include "../SEA/EllipticCurve.h"

EllipticCurve *E;
fmpz_t *lambda;
flint_rand_t state;

#define print_lambdas() do { \
	for (int q = 0; q < ECPROJ_LAMBDA_SIZE; ++q) { \
		printf("lambda[%d] = ", q); fmpz_print(lambda[q]); printf("\n"); \
	} \
} while (0)

void test_ECProjPoint_init() {
	fmpz_t p;
	fmpz_t A;
	fmpz_t B;

	fmpz_init(A);
	fmpz_init(B);
	fmpz_init(p),
	fmpz_set_ui(A, 663);
	fmpz_set_ui(B, 368);
	fmpz_set_ui(p, 1009);

	E = new_EllipticCurve(p, A, B);

	lambda = malloc(ECPROJ_LAMBDA_SIZE * sizeof(fmpz_t));
	for (int i = 0; i < ECPROJ_LAMBDA_SIZE; ++i) {
		fmpz_init(lambda[i]);
	}

	flint_randinit(state);

	fmpz_clear(A);
	fmpz_clear(B);
	fmpz_clear(p);
}

void test_ECProjPoint_clear() {
	free_EllipticCurve(E);
	for (int i = 0; i < ECPROJ_LAMBDA_SIZE; ++i) {
		fmpz_clear(lambda[i]);
	}
	free(lambda);
	flint_randclear(state);
}

char * test_ECProjPoint_newPoint() {
	ECProjPoint* P = new_ECProjPoint();


	t_assert(fmpz_equal_ui(P->x, 0), "P is not init at 0.");
	t_assert(fmpz_equal_ui(P->y, 0), "P is not init at 0.");
	t_assert(fmpz_equal_ui(P->z, 0), "P is not init at 0.");
	t_assert(P->is_zero, "P is not init at 0.");

	free_ECProjPoint(P);

	return NULL;
}

char * test_ECProjPoint_copyPoint() {
	ECProjPoint* P = new_ECProjPoint();

	fmpz_set_ui(P->x, 42);
	fmpz_set_ui(P->y, 39);
	fmpz_set_ui(P->z, 1);
	P->is_zero = 0;

	ECProjPoint* Q = new_ECProjPoint_copy(P);


	t_assert(fmpz_equal(P->x, Q->x), "Q is not equal to P.");
	t_assert(fmpz_equal(P->y, Q->y), "Q is not equal to P.");
	t_assert(fmpz_equal(P->z, Q->z), "Q is not equal to P.");
	t_assert(P->is_zero == Q->is_zero, "Q is not equal to P.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);

	return NULL;
}

char * test_ECProjPoint_equals() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* Q = new_ECProjPoint();

	// P = (42 : 39 : 1)
	fmpz_set_ui(P->x, 42);
	fmpz_set_ui(P->y, 39);
	fmpz_set_ui(P->z, 1);
	P->is_zero = 0;

	// Q = 2 * (42 : 39 : 1)
	fmpz_set_ui(Q->x, 84);
	fmpz_set_ui(Q->y, 78);
	fmpz_set_ui(Q->z, 2);
	Q->is_zero = 0;


	t_assert(EllipticCurve_projective_equals(P, Q, E), "Q is not equal to P.");

	free_ECProjPoint(P);
	P = new_ECProjPoint_copy(Q);

	t_assert(EllipticCurve_projective_equals(P, Q, E), "Q is not equal to P.");

	fmpz_set_ui(P->x, 7);

	t_assert(!EllipticCurve_projective_equals(P, Q, E), "Q is equal to P.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);

	return NULL;
}

char * test_ECProjPoint_isEquivalentTo() {
	ECProjPoint* P = new_ECProjPoint();
	ECPoint* a = new_ECPoint();

	// P = (420 : 390 : 10)
	fmpz_set_ui(P->x, 420);
	fmpz_set_ui(P->y, 390);
	fmpz_set_ui(P->z, 10);
	P->is_zero = 0;

	// a = (42 : 39)
	fmpz_set_ui(a->x, 42);
	fmpz_set_ui(a->y, 39);
	a->is_zero = 0;


	t_assert(EllipticCurve_projective_isEquivalentTo(P, a, E), "P is not equivalent to its equivalent.");

	fmpz_set_ui(a->y, 38);

	t_assert(!EllipticCurve_projective_isEquivalentTo(P, a, E), "P is equivalent to something not equivalent.");


	free_ECProjPoint(P);
	free_ECPoint(a);

	return NULL;

}

char * test_ECProjPoint_isOnCurve() {
	ECProjPoint* P = new_ECProjPoint();
	ECPoint* a = new_ECPoint();

	EllipticCurve_randomPoint(a, E, state);

	// P = a
	fmpz_set(P->x, a->x);
	fmpz_set(P->y, a->y);
	fmpz_set_ui(P->z, 1);
	P->is_zero = 0;


	t_assert(EllipticCurve_projective_isOnCurve(P, E), "P is not on the curve.");

	fmpz_set_ui(P->z, 3);

	t_assert(!EllipticCurve_projective_isOnCurve(P, E), "invalid P is on the curve.");


	free_ECPoint(a);
	free_ECProjPoint(P);

	return NULL;
}

char * test_ECProjPoint_randomPoint() {
	ECProjPoint* P = new_ECProjPoint();

	EllipticCurve_projective_randomPoint(P, E, state);


	t_assert(!P->is_zero, "Random point is zero.");
	t_assert(!fmpz_equal_ui(P->z, 0), "Random point has z = 0.");
	t_assert(!fmpz_equal_ui(P->x, 0) && !fmpz_equal_ui(P->y, 0), "Random point is (0, 0).");
	t_assert(EllipticCurve_projective_isOnCurve(P, E), "Random point is not on curve.");


	free_ECProjPoint(P);

	return NULL;
}

char * test_ECProjPoint_add() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* Q = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	ECProjPoint* tmp = new_ECProjPoint();


	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	// Q = (709 : 157 : 901) random point on E
	fmpz_set_ui(Q->x, 709);
	fmpz_set_ui(Q->y, 157);
	fmpz_set_ui(Q->z, 901);
	Q->is_zero = 0;


	EllipticCurve_projective_add_memory(R, P, Q, E, tmp, lambda);

	t_assert(fmpz_equal_ui(R->x, 514), "Added point R != P + Q.");
	t_assert(fmpz_equal_ui(R->y, 892), "Added point R != P + Q.");
	t_assert(fmpz_equal_ui(R->z, 989), "Added point R != P + Q.");
	t_assert(!(R->is_zero), "Added point is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(R, E), "Added point is not on curve.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);
	free_ECProjPoint(R);
	free_ECProjPoint(tmp);

	return NULL;
}

char * test_ECProjPoint_double() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	ECProjPoint* tmp = new_ECProjPoint();

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;


	EllipticCurve_projective_double_memory(R, P, E, tmp, lambda);

	t_assert(fmpz_equal_ui(R->x, 277), "Doubled point R != [2] P.");
	t_assert(fmpz_equal_ui(R->y, 385), "Doubled point R != [2] P.");
	t_assert(fmpz_equal_ui(R->z, 611), "Doubled point R != [2] P.");
	t_assert(!(R->is_zero), "Doubled point is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(R, E), "Doubled point is not on curve.");


	free_ECProjPoint(P);
	free_ECProjPoint(R);
	free_ECProjPoint(tmp);

	return NULL;
}

char * test_ECProjPoint_scalarMul() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	fmpz_t m;

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	fmpz_init(m);
	fmpz_set_ui(m, 593);


	EllipticCurve_projective_scalarMul(R, P, m, E);

	t_assert(fmpz_equal_ui(R->x, 718), "Multiplied point R != [593] P.");
	t_assert(fmpz_equal_ui(R->y, 160), "Multiplied point R != [593] P.");
	t_assert(fmpz_equal_ui(R->z, 808), "Multiplied point R != [593] P.");
	t_assert(!(R->is_zero), "Multiplied point is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(R, E), "Multiplied point is not on curve.");


	free_ECProjPoint(P);
	free_ECProjPoint(R);
	fmpz_clear(m);

	return NULL;
}

char * test_ECProjPoint_precomp_scalarMul() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	fmpz_t m;

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	fmpz_init(m);
	fmpz_set_ui(m, 593);

	ECProjPoint*** Ps = precomp_projective_multipleScalarMultiplications(P, 8, 10, E);


	projective_scalarMultiplications_withPrecomp(R, m, 8, Ps, E);

	t_assert(fmpz_equal_ui(R->x, 297), "Precomputed multiplied point R != [593] P.");
	t_assert(fmpz_equal_ui(R->y, 336), "Precomputed multiplied point R != [593] P.");
	t_assert(fmpz_equal_ui(R->z, 486), "Precomputed multiplied point R != [593] P.");
	t_assert(!(R->is_zero), "Precomputed multiplied point is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(R, E), "Precomputed multiplied point is not on curve.");


	free_precomp_projective_multipleScalarMultiplications(Ps, 8, 10);

	free_ECProjPoint(P);
	free_ECProjPoint(R);
	fmpz_clear(m);

	return NULL;
}

char * test_ECProjPoint_conversion() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* Q = new_ECProjPoint();
	ECPoint* a = new_ECPoint();

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;


	EllipticCurve_projective_to_affine(a, P, E);
	t_assert(fmpz_equal_ui(a->x, 577), "Conversion from projective to affine failed.");
	t_assert(fmpz_equal_ui(a->y, 803), "Conversion from projective to affine failed.");
	t_assert(!(a->is_zero), "Conversion from projective to affine led to 0.");

	EllipticCurve_affine_to_projective(Q, a, E);
	t_assert(EllipticCurve_projective_equals(P, Q, E), "Conversion from affine to projective failed.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);
	free_ECPoint(a);

	return NULL;
}


/****************************
 *     elaborated cases     *
 ****************************/


char * test_ECProjPoint_inPlaceAddition() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* Q = new_ECProjPoint();

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	// Q = (709 : 157 : 901) random point on E
	fmpz_set_ui(Q->x, 709);
	fmpz_set_ui(Q->y, 157);
	fmpz_set_ui(Q->z, 901);
	Q->is_zero = 0;


	// We add in place P += Q
	EllipticCurve_projective_add(P, P, Q, E);

	t_assert(fmpz_equal_ui(P->x, 514), "Added point P is not in place.");
	t_assert(fmpz_equal_ui(P->y, 892), "Added point P is not in place.");
	t_assert(fmpz_equal_ui(P->z, 989), "Added point P is not in place.");
	t_assert(!(P->is_zero), "Added point in place is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(P, E), "Added point in place is not on curve.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);

	return NULL;
}

char * test_ECProjPoint_inPlaceDoubling() {
	ECProjPoint* P = new_ECProjPoint();

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;


	// We double in place P = [2] P
	EllipticCurve_projective_double(P, P, E);

	t_assert(fmpz_equal_ui(P->x, 277), "Doubled point P is not in place.");
	t_assert(fmpz_equal_ui(P->y, 385), "Doubled point P is not in place.");
	t_assert(fmpz_equal_ui(P->z, 611), "Doubled point P is not in place.");
	t_assert(!(P->is_zero), "Doubled point in place is zero.");
	t_assert(EllipticCurve_projective_isOnCurve(P, E), "Doubled point in place is not on curve.");


	free_ECProjPoint(P);

	return NULL;
}

char * test_ECProjPoint_convertedAddition() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* Q = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	ECProjPoint* RR = new_ECProjPoint();

	ECPoint* a = new_ECPoint();
	ECPoint* b = new_ECPoint();
	ECPoint* c = new_ECPoint();
	ECPoint* cc = new_ECPoint();


	// we take random points in projective
	EllipticCurve_projective_randomPoint(P, E, state);
	EllipticCurve_projective_randomPoint(Q, E, state);

	// we convert them to affine
	EllipticCurve_projective_to_affine(a, P, E);
	EllipticCurve_projective_to_affine(b, Q, E);

	// we add the projective points, and convert the result in affine
	EllipticCurve_projective_add(R, P, Q, E);
	EllipticCurve_projective_to_affine(cc, R, E);

	// we add the affine points, and converts the result in projective
	EllipticCurve_add(c, a, b, E);
	EllipticCurve_affine_to_projective(RR, c, E);

	// we see if the results match
	t_assert(EllipticCurve_projective_equals(R, RR, E), "Projective addition converted failed.");
	t_assert(fmpz_equal(c->x, cc->x) && fmpz_equal(c->y, cc->y), "Projective addition converted failed.");


	free_ECProjPoint(P);
	free_ECProjPoint(Q);
	free_ECProjPoint(R);
	free_ECProjPoint(RR);

	free_ECPoint(a);
	free_ECPoint(b);
	free_ECPoint(c);
	free_ECPoint(cc);

	return NULL;
}

char * test_ECProjPoint_convertedDoubling() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	ECProjPoint* RR = new_ECProjPoint();

	ECPoint* a = new_ECPoint();
	ECPoint* c = new_ECPoint();
	ECPoint* cc = new_ECPoint();


	// we take one random point in projective
	EllipticCurve_projective_randomPoint(P, E, state);

	// we convert it to affine
	EllipticCurve_projective_to_affine(a, P, E);

	// we double the projective point, and convert the result in affine
	EllipticCurve_projective_double(R, P, E);
	EllipticCurve_projective_to_affine(cc, R, E);

	// we double the affine point, and converts the result in projective
	EllipticCurve_double(c, a, E);
	EllipticCurve_affine_to_projective(RR, c, E);

	// we see if the results match
	t_assert(EllipticCurve_projective_equals(R, RR, E), "Projective doubling converted failed.");
	t_assert(fmpz_equal(c->x, cc->x) && fmpz_equal(c->y, cc->y), "Projective doubling converted failed.");


	free_ECProjPoint(P);
	free_ECProjPoint(R);
	free_ECProjPoint(RR);

	free_ECPoint(a);
	free_ECPoint(c);
	free_ECPoint(cc);

	return NULL;
}

char * test_ECProjPoint_matching_precomp() {
	ECProjPoint* P = new_ECProjPoint();
	ECPoint* a = new_ECPoint();

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	// a = (577, 803) equivalent to P
	fmpz_set_ui(a->x, 577);
	fmpz_set_ui(a->y, 803);
	a->is_zero = 0;

	int word_size = 8;
	int exp_size = 10;


	ECProjPoint*** Pps = precomp_projective_multipleScalarMultiplications(P, word_size, exp_size, E);
	ECPoint*** Pas = precomp_multipleScalarMultiplications(a, word_size, exp_size, E);

	int n_words = (exp_size - 1)/word_size + 1;
	int n_patterns = (1 << word_size) - 1;

	for (int j = 0; j < n_patterns; ++j) {
		for (int i = 0; i < n_words; ++i) {
			t_assert(EllipticCurve_projective_isEquivalentTo(Pps[j][i], Pas[j][i], E), "A precomputed projective point does not match its equivalent in affine coordinates.");
		}
	}


	free_precomp_projective_multipleScalarMultiplications(Pps, 8, 10);
	free_precomp_multipleScalarMultiplications(Pas, 8, 10);

	return NULL;
}

char * test_ECProjPoint_convertedMultiplication() {
	ECProjPoint* P = new_ECProjPoint();
	ECProjPoint* R = new_ECProjPoint();
	ECPoint* a = new_ECPoint();
	ECPoint* b = new_ECPoint();
	fmpz_t m;

	// P = (369 : 218 : 861) random point on E
	fmpz_set_ui(P->x, 369);
	fmpz_set_ui(P->y, 218);
	fmpz_set_ui(P->z, 861);
	P->is_zero = 0;

	// a = (577, 803) equivalent to P
	fmpz_set_ui(a->x, 577);
	fmpz_set_ui(a->y, 803);
	a->is_zero = 0;

	fmpz_init(m);
	fmpz_set_ui(m, 593);

	ECProjPoint*** Ps = precomp_projective_multipleScalarMultiplications(P, 8, 10, E);
	ECPoint*** Pa = precomp_multipleScalarMultiplications(a, 8, 10, E);


	projective_scalarMultiplications_withPrecomp(R, m, 8, Ps, E);
	scalarMultiplications_withPrecomp(b, m, 8, Pa, E);

	t_assert(EllipticCurve_projective_isEquivalentTo(R, b, E), "Projective multiplication converted failed.");


	free_precomp_projective_multipleScalarMultiplications(Ps, 8, 10);
	free_precomp_multipleScalarMultiplications(Pa, 8, 10);

	free_ECProjPoint(P);
	free_ECProjPoint(R);
	free_ECPoint(a);
	free_ECPoint(b);
	fmpz_clear(m);

	return NULL;
}