/**
 * linvgamma.h
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 */

#ifndef R_IDENTIFIABILITY_LINVGAMMA_HEADER
#define R_IDENTIFIABILITY_LINVGAMMA_HEADER

#include "common.h"

double
c_linvgamma_s(double x, double a, double b);

void
c_linvgamma(double const *restrict x, R_xlen_t n, double a, double b, double *restrict r);

#endif
