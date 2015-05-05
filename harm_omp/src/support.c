/******************************************************************************
 *cr
 *cr         (C) Copyright 2010-2013 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "support.h"

void verify_d(double* input, double* output, unsigned num_elements) {

	const double relativeTolerance = 1e-12;
	const double absoluteTolerance = 1e-12;
	double absErr, relErr;
	bool fail = false;
	int fail_count = 0;
	int i;

	for(i = 0; i < num_elements; ++i) {
		absErr = fabs(input[i] - output[i]);
		relErr = fabs(absErr/input[i]);
		if (relErr > relativeTolerance && absErr > absoluteTolerance) {
			printf("TEST FAILED at i = %d, OMP = %e, orig = %e, rel_err = %e, abs_err = %e\n\n", i, input[i], output[i], relErr, absErr);
//			exit(0);
		fail = true; 
		fail_count++;
		} else {
//			printf("TEST FAILED at i = %d, cpu = %e, gpu = %e, rel_err = %e, abs_err = %e\n\n", i, input[i], output[i], relErr, absErr);
		}
	}
	if (fail) { printf("TEST FAILED on %d/%d elements.\n", fail_count, num_elements); }
//	else { printf("TEST PASSED\n"); }

}


void startTime(Timer* timer) {
	gettimeofday(&(timer->startTime), NULL);
}

void stopTime(Timer* timer) {
	gettimeofday(&(timer->endTime), NULL);
}

float elapsedTime(Timer timer) {
	return ((float) ((timer.endTime.tv_sec - timer.startTime.tv_sec) \
		+ (timer.endTime.tv_usec - timer.startTime.tv_usec)/1.0e6));
}

