/**
 * common.h
 * Copyright (c) 2017 Victhor S. Sartório
 * This file and its contents are licensed under the terms of the MIT License
 *
 * This is the common header which includes and defines basic facilities that shall be
 * included throughout every/most other C files in this project.
 */

#ifndef R_IDENTIFIABILITY_COMMON_HEADER
#define R_IDENTIFIABILITY_COMMON_HEADER

// ----------------------------------------------
// Includes
// ----------------------------------------------

#include <R.h>
#include <Rinternals.h>

// ----------------------------------------------
// Macros
// ----------------------------------------------

#if defined(__GNUC__) || defined(__clang__)
#define HOTFUNC __attribute__((hot))
#else
#define HOTFUNC
#endif

#endif
