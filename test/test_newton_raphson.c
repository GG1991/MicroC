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

const double strain[6] = { 1., 2., 3., 1., 1., 1. };


int main(void)
{
	MICROC_INST_START

	int ierr;
	int size[3] = { 10, 10, 10 };
	int ngp = 100;
	int type = 2;
	double params[1] = { 0.2 };

	material_t materials[NMATERIALS];
	material_set(&materials[0], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);
	material_set(&materials[1], 1.0e8, 0.25, 1.0e8, 1.0e4, 0);

	ierr = microc_init(ngp, size, type, params, materials);

	int i;
	for (i = 0; i < REPETITIONS; ++i) {
		VecZeroEntries(u[0]);
		newton_t newton;
		newton_raphson_v(A[0], A0[0], b[0], u[0], du[0], false, strain, NULL, PCJACOBI, KSPCG, &newton);
	}

	microc_finish();

	MICROC_INST_END
	MICROC_INST_PRINT

	return 0;
}
