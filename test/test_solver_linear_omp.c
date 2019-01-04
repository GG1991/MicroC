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


#define REPETITIONS 10


int main(int argc, char **argv)
{
	int ierr;
	int ngp = 1;
	int type = 2;
	double params[1] = { 0.2 };
	if (argc <= 1) {
		printf("Usage: %s <n>\n", argv[0]);
		return 1;
	}
	int n = atoi(argv[1]);
	int size[3] = { n, n, n };

	MICROC_INST_START

	material_t materials[NMATERIALS];
	material_set(&materials[0], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);

	microc_initv(ngp, size, type, params, materials, argc, argv);

	double time = omp_get_wtime();

	int i;
#pragma omp parallel for
	for (i = 0; i < REPETITIONS; ++i) {

		int thread_id = omp_get_thread_num();
		double err;
		int its;
		const double strain[6] = { 1., 2., 3., 1., 1., 1. };

		VecZeroEntries(u[thread_id]);
		bc_apply_on_u(u[thread_id], strain);

		err = assembly_res(b[thread_id], u[thread_id], NULL);
		printf("Thread : %d |NR err| = %lf\n", thread_id, err);
		solve_ts(A0[thread_id], b[thread_id], du[thread_id], PCJACOBI, KSPCG, &its, &err);
		VecAXPY(u[thread_id], 1., du[thread_id]);

		err = assembly_res(b[thread_id], u[thread_id], NULL);
		printf("Thread : %d |NR err| = %lf\n", thread_id, err);
	}

	time = omp_get_wtime() - time;
	printf("time = %lf\n", time);
	microc_finish();

	MICROC_INST_END
	MICROC_INST_PRINT

	return 0;
}
