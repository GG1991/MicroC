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
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "petscksp.h"
#include "petscdm.h"
#include "petscdmda.h"

#include "gp.h"
#include "material.h"

#define nod_index(i,j,k)   ((k) * nx * ny + (j) * nx + (i))

#define NGP            8
#define NPE            8
#define NVOI           6
#define DIM            3
#define NEWTON_MIN_TOL 1.0e-1
#define NEWTON_REL_TOL 1.0e-3
#define NEWTON_MAX_ITS 2

#define LX             1.0
#define LY             1.0
#define LZ             1.0

enum {
       	MIC_SPHERE,
       	MIC_LAYER_Y,
       	MIC_CILI_FIB_Z,
       	MIC_CILI_FIB_XZ,
       	MIC_QUAD_FIB_XYZ,
       	MIC_QUAD_FIB_XZ,
       	MIC_QUAD_FIB_XZ_BROKEN_X
};

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

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";

int first_init; // variable to know if init was called

double lx, ly, lz, dx, dy, dz;
double wg;

double newton_min_tol;
int newton_max_its;

int nx, ny, nz, nn;
int nndim;
int nelem;
int nex, ney, nez;
int ngp;
int nvars;

int *elem_type;
int micro_type;
double special_param;

int nbcs;
int *bc_index;
double *bc_value;

const int *eix;

gp_t *gp_list;
material_t material_list[2];

DM da;
PC pc;
KSP ksp;
Mat A, A0;
Vec u, du, b;

// init.c
int microc_init(const int ngp, const int size[3], const int micro_type,	const double *micro_params);
int microc_finish(void);

// homogenize.c
void microc_set_macro_strain(const int gp_id, const double *macro_strain);
void microc_get_macro_strain(const int gp_id, double *macro_strain);
void microc_get_macro_stress(const int gp_id, double *macro_stress);
void microc_get_macro_ctan(const int gp_id, double *macro_ctan);
void microc_homogenize();

//general.c
void calc_strain(const double u_e[NPE * DIM], const double B[NVOI][NPE * DIM], double strain[6]);
void calc_stress(const int ie, const double strain[NVOI], const double *vars_old, double stress[NVOI]);
void calc_ctan(const int ie, const double strain[NVOI], const double *vars_old, double ctan[NVOI][NVOI]);
void calc_B(int gp, double B[6][NPE * DIM]);
int get_elem_type(const int ex, const int ey, const int ez);
material_t *get_material(const int ie);
void calc_elemental_displacements_with_ie(const int ie, const double *u_global, double u_e[NPE * DIM]);

// assembly.c
void get_elem_rhs_with_ie(const int ie, const double *u, const double *varsold, double be[NPE * DIM]);
void get_elem_mat_with_ie(const int ie, const double *u_global, const double *varsold, double be[NPE * DIM]);
double assembly_res(Vec b, Vec u, double *vars_old);
int assembly_jac(Mat A, Vec u, double *vars_old);

// bcs.c
void mat_vec(const double strain_mat[3][3], const double coor[3], double disp[3]);
int bc_apply_on_u(Vec u, const double strain[6]);

#endif
