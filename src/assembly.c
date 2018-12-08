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

void get_elem_rhs_with_ie(const int ie, const double *u,
		  const double *varsold, double be[npe * dim])
{
	double stress[NVOI], strain[NVOI];
	double B[NVOI][NPE * DIM];

	memset(be, 0, NPE * DIM * sizeof(double));

	int gp;
	for (gp = 0; gp < npe; ++gp) {

		calc_B(gp, B);

		calc_strain(u_e, B, strain);
		get_stress(ie, strain, varsold, stress);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += B[j][i] * stress[j] * wg;
	}
}
