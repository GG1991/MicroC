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

void get_elem_rhs_with_ie(const int ie, const double *u_arr, const double *vars_old, double be[NPE * DIM])
{
	double stress[NVOI], strain[NVOI];
	double B[NVOI][NPE * DIM];
	double u_e[NPE * DIM];

	memset(be, 0, NPE * DIM * sizeof(double));

	calc_elemental_displacements_with_ie(ie, u_arr, u_e);

	int gp;
	for (gp = 0; gp < NPE; ++gp) {

		calc_B(gp, B);

		calc_strain(u_e, B, strain);
		calc_stress(ie, strain, vars_old, stress);

		for (int i = 0; i < NPE * DIM; ++i)
			for (int j = 0; j < NVOI; ++j)
				be[i] += B[j][i] * stress[j] * wg;
	}
}


void get_elem_mat_with_ie(const int ie, const double *u_arr, const double *vars_old, double Ae[NPE * DIM * NPE * DIM])
{
	int i, j, k, m;
	double strain[NVOI], ctan[NVOI][NVOI];
	double B[NVOI][NPE * DIM];
	double cxb[NVOI][NPE * DIM];
	double u_e[NPE * DIM];

	memset(Ae, 0, NPE * DIM * NPE * DIM * sizeof(double));

	calc_elemental_displacements_with_ie(ie, u_arr, u_e);

	int gp;
	for (gp = 0; gp < NPE; ++gp) {

		calc_B(gp, B);

		calc_strain(u_e, B, strain);
		calc_ctan(ie, strain, vars_old, ctan);

		for (i = 0; i < NVOI; ++i) {
			for (j = 0; j < NPE * DIM; ++j) {
				double tmp = 0.0;
				for (k = 0; k < NVOI; ++k)
					tmp += ctan[i][k] * B[k][j];
				cxb[i][j] = tmp;
			}
		}

		for (i = 0; i < NPE * DIM; ++i) {
			for (m = 0; m < NVOI; ++m) {
				for (j = 0; j < NPE * DIM; ++j) {
					Ae[i * (NPE * DIM) + j] += B[m][i] * cxb[m][j] * wg;
				}
			}
		}
	}
}


double assembly_res(Vec b, Vec u, const double *vars_old)
{
	int ierr;
	int i, ie, n, d;
	int ix[NPE * DIM];
	double be[NPE * DIM * NPE * DIM];
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

	for(i = 0; i < nbcs; ++i)
		b_arr[bc_index[i]] = 0.0;

	ierr = VecRestoreArray(b, &b_arr); CHKERRQ(ierr);
	ierr = VecRestoreArray(u, &u_arr); CHKERRQ(ierr);

	ierr = VecScale(b, -1.); CHKERRQ(ierr);
	//VecView(b, PETSC_VIEWER_STDOUT_WORLD);

	double norm;
	VecNorm(b, NORM_2, &norm);

	return norm;
}


int assembly_jac(Mat A, Vec u, const double *vars_old)
{
	int ierr;
	int ie, n, d;
	int ix[NPE * DIM];
	double Ae[NPE * DIM * NPE * DIM];
	double *u_arr;

	ierr = VecGetArray(u, &u_arr);CHKERRQ(ierr);
	ierr = MatZeroEntries(A); CHKERRQ(ierr);

	for(ie = 0; ie < nelem; ++ie) {

		// Calculates the elemental matrix <Ae>
		get_elem_mat_with_ie(ie, u_arr, vars_old, Ae);

		// Calculates the indeces where <Ae> should be assembly
		for(n = 0; n < NPE; ++n)
			for(d = 0; d < DIM; ++d)
				ix[n * DIM + d] = eix[ie * NPE + n] * DIM + d;

		// Sum <Ae> over <A>
		ierr = MatSetValuesLocal(A, NPE * DIM, ix, NPE * DIM, ix, Ae, ADD_VALUES); CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = MatZeroRowsColumns(A, nbcs, bc_index, 1.0, NULL, NULL); CHKERRQ(ierr);
	ierr = VecRestoreArray(u, &u_arr); CHKERRQ(ierr);
	MatView(A, PETSC_VIEWER_DRAW_WORLD);

	return ierr;
}
