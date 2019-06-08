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
#include <time.h>


const double strain[6] = { 1., 2., 3., 1., 1., 1. };


int main(int argc, char **argv)
{
	int ierr;
	int ngp = 1;
	int type = MIC_CILI_FIB_Z;
	double params[1] = { 0.2 };
	if (argc < 3) {
		printf("Usage: %s <n> <a>\n", argv[0]);
		return 1;
	}
	const int n = atoi(argv[1]);
	const double a = atof(argv[2]);
	int size[3] = { n, n, n };

	double Em = 1.0e8;

	material_t materials[NMATERIALS];
	material_set(&materials[0], Em, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], a * Em, 0.25, 1.0e8, 1.0e4, 0);
	material_print(&materials[0]);
	material_print(&materials[1]);

	microc_initv(ngp, size, type, params, materials, argc, argv);

	double err;
	int its;
	VecZeroEntries(u[0]);
	bc_apply_on_u(u[0], strain);

//	VecView(u[0], PETSC_VIEWER_STDOUT_WORLD);

	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);
	printf("SOLVER_START\n");

	clock_t start = clock();

	solve(A0[0], b[0], du[0], &its, &err);

	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf("\nn = %d\t time = %lf\t iter = %d\n", n * n * n, seconds, its);

	printf("SOLVER_END\n");
	VecAXPY(u[0], 1., du[0]);

	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);

	microc_finish();

	return 0;
}
