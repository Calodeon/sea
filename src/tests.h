#include <stdio.h>

/***************************
 *      Test framework     *
 ***************************/

#ifndef TESTS_H_
#define TESTS_H_

#define t_assert(cond, error) do { \
	if (!(cond)) { \
		return error; \
	} \
} while (0)

extern int error;

#define t_run_unit(unit) do { \
	char * result = unit(); \
	if (result) { \
		printf("\n\t[ERROR] %s", result); \
		error = 1; \
	} \
} while (0)

#define t_run_test(test, name) do { \
	printf("Testing %s .... ", name); \
	error = test(); \
	if (!error) { \
		printf("PASS\n"); \
	} \
} while (0)

#endif /* TESTS_H_ */