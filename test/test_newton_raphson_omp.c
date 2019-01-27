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


#define REPETITIONS 5


int main(void)
{
	MICROC_INST_START

	int ierr;
	int size[3] = { 50, 50, 50 };
	int ngp = 100;
	int type = 2;
	double params[1] = { 0.2 };

	material_t materials[NMATERIALS];
	material_set(&materials[0], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);

	ierr = microc_init(ngp, size, type, params, materials);

	double time = omp_get_wtime();

	int i;
#pragma omp parallel for
	for (i = 0; i < REPETITIONS; ++i) {
		const double strain[6] = { 1., 2., 3., 1., 1., 1. };
		int thread_id = omp_get_thread_num();
		Mat _A = A[thread_id];
		Mat _A0 = A0[thread_id];
		Vec _b = b[thread_id];
		Vec _u = u[thread_id];
		Vec _du = du[thread_id];
		VecZeroEntries(_u);
		newton_t newton;
		newton_raphson_v(_A, _A0, _b, _u, _du, 
				 false, strain, NULL,
				 PCJACOBI, KSPCG,
				 &newton, true);
	}
	time = omp_get_wtime() - time;
	printf("time = %lf\n", time);

	microc_finish();

	MICROC_INST_END
	MICROC_INST_PRINT

	return 0;
}
