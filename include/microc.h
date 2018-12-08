/*
 *  This source code is part of MicroC: a finite element code
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MICROC_H
#define MICROC_H

#include <stdio.h>

#include "petscksp.h"
#include "petscdm.h"
#include "petscdmda.h"

#define NGP            8
#define NPE            8
#define NVOI           6
#define DIM            3
#define NEWTON_MIN_TOL 1.0e-1
#define NEWTON_MAX_ITS 2

#define FINAL_TIME     10.0
#define TIME_STEPS     1
#define VTU_FREQ       -1
#define DT             0.001
#define NX             5
#define NY             5
#define NZ             5
#define LX             10.0
#define LY             1.0
#define LZ             1.0

#define U_MAX          -0.8

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";

#define CONSTXG        0.577350269189626

static double xg[8][3] = {
	{ -CONSTXG, -CONSTXG, -CONSTXG },
	{ +CONSTXG, -CONSTXG, -CONSTXG },
	{ +CONSTXG, +CONSTXG, -CONSTXG },
	{ -CONSTXG, +CONSTXG, -CONSTXG },
	{ -CONSTXG, -CONSTXG, +CONSTXG },
	{ +CONSTXG, -CONSTXG, +CONSTXG },
	{ +CONSTXG, +CONSTXG, +CONSTXG },
	{ -CONSTXG, +CONSTXG, +CONSTXG } };

double lx, ly, lz, dx, dy, dz;
double wg;

PetscReal dt, final_time;
PetscReal newton_min_tol;
PetscInt ts;
PetscInt vtu_freq;
PetscInt newton_max_its;

DM da;
PC pc;
KSP ksp;
Mat A;
Vec u, du, b;

PetscErrorCode microc_init();
PetscErrorCode microc_finish();

#endif
