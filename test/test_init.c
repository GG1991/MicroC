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
	int size1[3] = { 10, 10, 10 };
	int ngp1 = 1;
	int type1 = 2;
	double params1[1] = { 0.2 };

	microc_init(ngp1, size1, type1, params1);
	microc_finish();

	int size2[3] = { 20, 20, 20 };
	int ngp2 = 10;
	int type2 = 5;
	double params2[1] = { 0.3 };

	microc_init(ngp2, size2, type2, params2);
	microc_finish();

	return 0;
}
