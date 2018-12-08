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

#include "microc.h"

void get_elem_rhs_with_ie(const int ie, const double *u_global, const double *varsold, double be[NPE * DIM])
{
	double stress[NVOI], strain[NVOI];
	double B[NVOI][NPE * DIM];
	double u_e[NPE * DIM];

	memset(be, 0, NPE * DIM * sizeof(double));

	calc_elemental_displacements_with_ie(ie, u_global, u_e);

	int gp;
	for (gp = 0; gp < NPE; ++gp) {

		calc_B(gp, B);

		calc_strain(u_e, B, strain);
		calc_stress(ie, strain, varsold, stress);

		for (int i = 0; i < NPE * DIM; ++i)
			for (int j = 0; j < NVOI; ++j)
				be[i] += B[j][i] * stress[j] * wg;
	}
}


double assembly_res(Vec u, Vec b, double *vars_old)
{
	int ierr;
	int i, j, k;
	int n, d, gp;
	int ix[NPE * DIM];
	double u_e[NPE * DIM];
	double stress[NVOI];
	double be[NPE * DIM * NPE * DIM];
	double B[NVOI][NPE * DIM];
	int ie;
	double *u_arr, *b_arr;

	Vec u_loc;

	ierr = VecZeroEntries(b); CHKERRQ(ierr);

	ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
	ierr = VecGetArray(u, &u_arr);CHKERRQ(ierr);

	for(ie = 0; ie < nelem; ++ie) {

		// Calculates the elemental residue <be>
		get_elem_rhs_with_ie(ie, u_arr, vars_old, be);

		// Calculates the indeces where <be> should be assembly
		for(n = 0; n < NPE; ++n)
			for(d = 0; d < DIM; ++d)
				ix[n * DIM + d] = eix[ie * NPE + n] * DIM + d;

		// Sum <be> over <b_arr>
		for(n = 0; n < NPE * DIM; ++n)
			b_arr[ix[n]] += be[n];
	}

	ierr = VecRestoreArray(b, &b_arr); CHKERRQ(ierr);
	ierr = VecRestoreArray(u, &u_arr); CHKERRQ(ierr);

	//VecView(b, PETSC_VIEWER_STDOUT_WORLD);

//	ierr = apply_bc_on_res(b);

	ierr = VecScale(b, -1.); CHKERRQ(ierr);

	double norm;
	VecNorm(b, NORM_2, &norm);

	return norm;
}
