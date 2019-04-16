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


#define REPETITIONS 1

const double strain[6] = { 1., 2., 3., 1., 1., 1. };

int SOLVER_0(void)
{

	double err;
	int its;
	VecZeroEntries(u[0]);
	bc_apply_on_u(u[0], strain);
	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);
	MICROC_INST_START
	printf("SOLVER_START\n");
	solve(A0[0], b[0], du[0], &its, &err);
	printf("SOLVER_END\n");
	MICROC_INST_END
	VecAXPY(u[0], 1., du[0]);
	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);
	return its;
}

int SOLVER(void)
{

	double err;
	int its;
	VecZeroEntries(u[0]);
	bc_apply_on_u(u[0], strain);

//	VecView(u[0], PETSC_VIEWER_STDOUT_WORLD);

	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);
	MICROC_INST_START
	printf("SOLVER_START\n");
	solve(A0[0], b[0], du[0], &its, &err);
	printf("SOLVER_END\n");
	MICROC_INST_END
	VecAXPY(u[0], 1., du[0]);

	err = assembly_res(b[0], u[0], NULL);
	printf("|NR err| = %lf\n", err);
	return its;
}


int main(int argc, char **argv)
{
	int ierr;
	int ngp = 1;
	int type = MIC_CILI_FIB_Z;
	double params[1] = { 0.2 };
	if (argc <= 1) {
		printf("Usage: %s <n> <a>\n", argv[0]);
		return 1;
	}
	const int n = (argc >= 2) ? atoi(argv[1]) : 10;
	const double a = (argc >= 3) ? atof(argv[2]) : 1.0;
	int size[3] = { n, n, n };

	MICROC_INST_START;

	double Em = 1.0e8;

	material_t materials[NMATERIALS];
	material_set(&materials[0], Em, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], a * Em, 0.25, 1.0e8, 1.0e4, 0);
	material_print(&materials[0]);
	material_print(&materials[1]);

	microc_initv(ngp, size, type, params, materials, argc, argv);

	int its;
//	its = SOLVER_0();

	int i;
	for (i = 0; i < REPETITIONS; ++i)
		its = SOLVER();

	double total = get_total_time(7);
	int calls = get_total_calls(7);

	printf("\nn = %d\t mean = %lf\t iter = %d\n", n * n * n, total / calls, its);

	microc_finish();

	MICROC_INST_END
	MICROC_INST_PRINT

	return 0;
}