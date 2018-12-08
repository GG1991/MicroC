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
	int size[3] = { 10, 10, 10 };
	int ngp = 100;
	int type = 2;
	double params[1] = { 0.2 };

	microc_init(ngp, size, type, params);

	const double strain[6] = { 0.2, 0.1, -0.3, 5.7, 2.0, -4.9 };
	double strain_get[6];

	int gp, i;
	for (gp = 0; gp < ngp; ++gp) {
		microc_set_macro_strain(gp, strain);
		microc_get_macro_strain(gp, strain_get);
		for (i = 0; i < 6; ++i) {
			assert(fabs(strain[i] - strain_get[i]) < 1.0e-5);
		}
	}

	microc_finish();

	return 0;
}
