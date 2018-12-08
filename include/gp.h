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

#ifndef GP_H
#define GP_H

#include "stdbool.h"

#define GP_NR_MAX_ITS 10

typedef struct {

	double macro_strain[6];
	double macro_stress[6];
	double *macro_ctan;

	bool allocated; // flag for memory optimization

	double *int_vars_n; // vectors for calculations
	double *int_vars_k;
	double *u_n;
	double *u_k;

	int sigma_solver_its[GP_NR_MAX_ITS];
	int sigma_newton_its;
	double sigma_solver_err[GP_NR_MAX_ITS];
	double sigma_newton_err[GP_NR_MAX_ITS];
	long int sigma_cost;
} gp_t;

#endif
