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


void mat_vec(const double strain_mat[3][3], const double coor[3], double disp[3])
{
	int i, j;
	for (i = 0; i < 3; ++i) {
		disp[i] = 0.0;
		for (j = 0; j < 3; ++j) {
			disp[i] += strain_mat[i][j] * coor[j];
		}
	}
	return;
}


int bc_apply_on_u(Vec u, const double strain[6])
{
	int ierr;
	int i, j, k, d;
	int count = 0;
	double disp[3];
	double strain_mat[3][3] = {
		{      strain[0], .5 * strain[3], .5 * strain[4] },
		{ .5 * strain[3],      strain[1], .5 * strain[5] },
		{ .5 * strain[4], .5 * strain[5],      strain[2] } };

	/* Y = 0 */
	for (i = 0; i < nx; ++i) {
		for (k = 0; k < nz; ++k) {
			const double coor[3] = { i * dx, 0 * dy, k * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	/* Y = LY */
	for (i = 0; i < nx; ++i) {
		for (k = 0; k < nz; ++k) {
			const double coor[3] = { i * dx, (ny - 1) * dy, k * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	/* Z = 0 */
	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			const double coor[3] = { i * dx, j * dy, 0 * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	/* Z = LZ */
	for (i = 0; i < nx; ++i) {
		for (k = 0; k < ny; ++k) {
			const double coor[3] = { i * dx, j * dy, (nz - 1) * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	/* X = 0 */
	for (j = 0; j < ny; ++j) {
		for (k = 0; k < nz; ++k) {
			const double coor[3] = { 0 * dx, j * dy, k * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	/* X = LX */
	for (j = 0; j < ny; ++j) {
		for (k = 0; k < nz; ++k) {
			const double coor[3] = { (nx - 1) * dx, j * dy, k * dz };
			mat_vec(strain_mat, coor, disp);
			for (d = 0; d < DIM; ++d)
				bc_value[count++] = disp[d];
		}
	}

	ierr = VecSetValues(u, nbcs, bc_index, bc_value, INSERT_VALUES); CHKERRQ(ierr);

	return ierr;
}
