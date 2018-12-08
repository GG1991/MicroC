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


int get_elem_type(const int ex, const int ey, const int ez)
{
	int i;
	const double coor[3] = {
		ex * dx + dx / 2.,
		ey * dy + dy / 2.,
		ez * dz + dz / 2. }; // 2D -> dz = 0

	if (micro_type == MIC_SPHERE) { // sphere in the center

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (i = 0; i < DIM; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == MIC_LAYER_Y) { // 2 flat layers in y dir

		const double width = special_param;
		return (coor[1] < width);

	} else if (micro_type == MIC_CILI_FIB_Z) { // a cilindrical fiber in z dir

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (i = 0; i < 2; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == MIC_CILI_FIB_XZ) { // 2 cilindrical fibers one in x and z dirs

		const double rad = special_param;
		const double cen_1[3] = { lx / 2., ly * .75, lz / 2. };
		double tmp_1 = 0.;
		for (i = 0; i < 2; ++i)
			tmp_1 += (cen_1[i] - coor[i]) * (cen_1[i] - coor[i]);

		const double cen_2[3] = { lx / 2., ly * .25, lz / 2. };
		double tmp_2 = 0.;
		for (i = 1; i < 3; ++i)
			tmp_2 += (cen_2[i] - coor[i]) * (cen_2[i] - coor[i]);

		return ((tmp_1 < rad * rad) || (tmp_2 < rad * rad));

	} else if (micro_type == MIC_QUAD_FIB_XYZ) {

	       	/* 3 quad fibers in x, y and z dirs */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XZ) {

	       	/* 2 quad fibers in x and z dirs */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XZ_BROKEN_X) {

	       	/* 2 quad fibers in x and z dirs and the one in x is broken */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width &&
		    (coor[0] < lx * .8 || coor[0] > lx * .9))
			return 1;

		return 0;
	}

	return -1;
}


material_t *get_material(const int ie)
{
	return &material_list[elem_type[ie]];
}


void calc_stress(const int ie, const double strain[NVOI], const double *vars_old, double stress[NVOI])
{
	material_t *mat = get_material(ie);

	if (mat->plasticity == true) {

		const double *eps_p_old = &vars_old[0];
		const double alpha_old = vars_old[6];

		plastic_get_stress(mat, strain, eps_p_old, alpha_old, stress);

	} else {

		isolin_get_stress(mat, strain, stress);
	}

}


void calc_elemental_displacements_with_ie(const int ie, const double *u_global, double u_e[NPE * DIM])
{
	// Calculates the elemental displacement <u_e> using the global element
	// numeration <ie>

	int n, d;
	for(n = 0; n < NPE; ++n)
		for(d = 0; d < DIM; ++d)
			u_e[n * DIM + d] = u_global[eix[ie * NPE + n] * DIM + d];
	return;
}


void calc_strain(const double u_e[NPE * DIM], const double B[NVOI][NPE * DIM], double strain[6])
{
	/* Calculates the strain using the elemental displacements <u_e> and <B>
	 * strain = B * u_e
	*/
	int i, j;
	for(i = 0; i < NVOI; ++i)
		for(j = 0; j < NPE * DIM; ++j)
			strain[i] += B[i][j] * u_e[j];
	return;
}


void calc_B(int gp, double B[6][NPE * DIM])
{
	int i;

	const double dsh[NPE][DIM] = {
		{
			-(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz } };

	for (i = 0; i < NPE; ++i) {
		B[0][i * DIM    ] = dsh[i][0];
		B[0][i * DIM + 1] = 0;
		B[0][i * DIM + 2] = 0;
		B[1][i * DIM    ] = 0;
		B[1][i * DIM + 1] = dsh[i][1];
		B[1][i * DIM + 2] = 0;
		B[2][i * DIM    ] = 0;
		B[2][i * DIM + 1] = 0;
		B[2][i * DIM + 2] = dsh[i][2];
		B[3][i * DIM    ] = dsh[i][1];
		B[3][i * DIM + 1] = dsh[i][0];
		B[3][i * DIM + 2] = 0;
		B[4][i * DIM    ] = dsh[i][2];
		B[4][i * DIM + 1] = 0;
		B[4][i * DIM + 2] = dsh[i][0];
		B[5][i * DIM    ] = 0;
		B[5][i * DIM + 1] = dsh[i][2];
		B[5][i * DIM + 2] = dsh[i][1];
	}
}
