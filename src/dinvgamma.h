/**
 * dinvgamma.h
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 */

#ifndef R_IDENTIFIABILITY_INVGAMMA_HEADER
#define R_IDENTIFIABILITY_INVGAMMA_HEADER

#include "common.h"

double
c_dinvgamma_s(double x, R_xlen_t n, double a, double b);

void
c_dinvgamma(double const *restrict x, R_xlen_t n, double a, double b, double *restrict r);

#endif
