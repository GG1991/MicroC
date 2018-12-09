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

int main(void)
{
	int ierr;
	int size[3] = { 2, 2, 2 };
	int ngp = 100;
	int type = 2;
	double params[1] = { 0.2 };

	ierr = microc_init(ngp, size, type, params);

	const double strain[6] = { 1., 2., 3., 1., 1., 1. };
	const double u_exact[24] = {
		0.00000, 0.00000, 0.00000,
		1.00000, 0.50000, 0.50000,
		0.50000, 2.00000, 0.50000,
		1.50000, 2.50000, 1.00000,
		0.50000, 0.50000, 3.00000,
		1.50000, 1.00000, 3.50000,
		1.00000, 2.50000, 3.50000,
		2.00000, 3.00000, 4.00000 };

	ierr = bc_apply_on_u(u, strain);
	VecView(u, PETSC_VIEWER_STDOUT_WORLD);

	int nn = size[0] * size[1] * size[2];
	double *u_arr;
	ierr = VecGetArray(u, &u_arr);CHKERRQ(ierr);
	int i;
	for (i = 0; i < nn * DIM; ++i)
		assert(fabs(u_arr[i] - u_exact[i]) < 1.0e-5);
	ierr = VecRestoreArray(u, &u_arr); CHKERRQ(ierr);

	microc_finish();

	ierr = microc_init(ngp, size, type, params);
	microc_finish();

	return 0;
}
