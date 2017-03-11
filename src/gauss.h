/**
 * gauss.h
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 */

#ifndef R_IDENTIFIABILITY_GAUSS_HEADER
#define R_IDENTIFIABILITY_GAUSS_HEADER

#include "common.h"

double
c_lgauss_s(double x, double mu, double s2);

void
c_lgauss(double const *restrict x, R_xlen_t n, double mu, double s2, double *restrict r);

#endif
