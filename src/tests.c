#include <stdio.h>
#include <stdlib.h>
#include "tests.h"

#include "tests/EllipticCurveProjectiveArithmetic_tests.h"

/***************************
 *      Test functions     *
 ***************************/

// run tests/EllipticCurveProjectiveArithmetic_tests.c
int test_ECProjPoint() {

	int error = 0;

	test_ECProjPoint_init();

	t_run_unit(test_ECProjPoint_newPoint);
	t_run_unit(test_ECProjPoint_copyPoint);
	t_run_unit(test_ECProjPoint_equals);
	t_run_unit(test_ECProjPoint_isEquivalentTo);
	t_run_unit(test_ECProjPoint_isOnCurve);
	t_run_unit(test_ECProjPoint_randomPoint);
	t_run_unit(test_ECProjPoint_add);
	t_run_unit(test_ECProjPoint_double);
	t_run_unit(test_ECProjPoint_scalarMul);
	t_run_unit(test_ECProjPoint_precomp_scalarMul);
	t_run_unit(test_ECProjPoint_conversion);
	t_run_unit(test_ECProjPoint_inPlaceAddition);
	t_run_unit(test_ECProjPoint_inPlaceDoubling);
	t_run_unit(test_ECProjPoint_convertedAddition);
	t_run_unit(test_ECProjPoint_convertedDoubling);
	t_run_unit(test_ECProjPoint_matching_precomp);
	t_run_unit(test_ECProjPoint_convertedMultiplication);

	test_ECProjPoint_clear();

	return error;
}


// All tests
int all_tests() {

	int error = 0;

	t_run_test(test_ECProjPoint, "ECProjPoint");

	return error;
}

/***************************
 *           main          *
 ***************************/

int main(void) {

	int error = all_tests();

	if (!error) {
		printf("\nALL TESTS PASSED !\n");
		return 0;
	}

	printf("\n\nSOME TESTS FAILED.\n");
	return -1;
}