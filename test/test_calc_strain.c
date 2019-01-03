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
 *
 */


#include "microc.h"

int main(void)
{
	int ierr;
	int size[3] = { 2, 2, 2 };
	int ngp = 100;
	int type = 2;
	double params[1] = { 0.2 };

	material_t materials[NMATERIALS];
	material_set(&materials[0], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], 1.0e8, 0.25, 1.0e8, 1.0e4, 1);

	ierr = microc_init(ngp, size, type, params, materials);

	const double strain_macro[6][6] = {
	       	{1., 0., 0., 0., 0., 0. },
	       	{0., 1., 0., 0., 0., 0. },
	       	{0., 0., 1., 0., 0., 0. },
	       	{0., 0., 0., 1., 0., 0. },
	       	{0., 0., 0., 0., 1., 0. },
	       	{0., 0., 0., 0., 0., 1. }
	};

	double strain[6];

	double B[NVOI][NPE * DIM];
	double u_e[NPE * DIM];
	double *u_arr;

	int i, j, e = 0, gp = 0;

	ierr = VecGetArray(u[0], &u_arr); CHKERRQ(ierr);

	for (i = 0; i < 6; ++i) {

		ierr = bc_apply_on_u(u[0], strain_macro[i]);
		calc_elemental_displacements_with_ie(e, u_arr, u_e);

		calc_B(gp, B);
		calc_strain(u_e, B, strain);
		for (j = 0; j < 6; ++j)
			assert(fabs(strain[j] - strain_macro[i][j]) < 1.0e-5);
	}

	ierr = VecRestoreArray(u[0], &u_arr); CHKERRQ(ierr);

	microc_finish();

	return 0;
}
