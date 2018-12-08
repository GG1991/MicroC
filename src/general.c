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

void calc_stress(int nmat, const double strain[NPE * DIM])

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


void calc_strain(const double u_e[NPE * DIM], int gp, double strain[6])
{
	// Calculates the strain using the elemental displacements <u_e>

	double B[NVOI][NPE * DIM];
	calc_B(gp, B);

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
