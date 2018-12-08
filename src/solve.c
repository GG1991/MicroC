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


int solve(Mat A, Vec b, Vec x, double *_err)
{
}


int newton_raphson(bool non_linear, double strain[NVOI], double *vars_old, 
		   Vec u, double _newton_err[NEWTON_MAX_ITS],
		   int _solver_its[NEWTON_MAX_ITS], double _solver_err[NEWTON_MAX_ITS])
{
	//set_displ_bc(strain, u);

	int ierr, lits = 0;
	double lerr, lerr0;
	int solver_its;
	double solver_err;

	if (_solver_its != NULL) memset(_solver_its, 0, NEWTON_MAX_ITS * sizeof(int));
	if (_solver_err != NULL) memset(_solver_err, 0, NEWTON_MAX_ITS * sizeof(double));
	if (_newton_err != NULL) memset(_newton_err, 0, NEWTON_MAX_ITS * sizeof(double));

	while (lits < NEWTON_MAX_ITS) {

		lerr = assembly_res(u, b, vars_old);

		if (lits == 0)
			lerr0 = lerr;

		if (_newton_err != NULL)
		       	_newton_err[lits] = lerr;

		if (lerr < NEWTON_MIN_TOL || lerr < lerr0 * NEWTON_REL_TOL)
			break;

		if (non_linear || lits > 0) {

			assembly_jac(A, u, vars_old);
			solver_its = solve(A, b, du, &solver_err);

		} else {

			solver_its = solve(A0, b, du, &solver_err);

		}

		if (_solver_its != NULL) _solver_its[lits] = solver_its;
		if (_solver_err != NULL) _solver_err[lits] = solver_err;

		ierr = VecAXPY(u, 1., du); CHKERRQ(ierr);

		lits++;

	}

	return lits;
}
