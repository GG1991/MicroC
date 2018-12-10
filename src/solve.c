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
	int ierr;
	int its = -1;
	double norm = -1;
			PetscViewer viewer;

	switch (solver) {
		case PETSC_CG:

			ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
			ierr = KSPSetUp(ksp); CHKERRQ(ierr);

			ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

			ierr = KSPGetIterationNumber(ksp, &its);
			ierr = KSPGetResidualNorm(ksp, &norm);
			break;

		case (CGPD):
			// MicroC solver implementation of the CG with
			// Diagonal Preconditioner
			break;
		default:
			ierr = 1;
			break;
	}
	printf("KSP : |Ax - b|/|Ax| = %e\tIts = %d\n", norm, its);
	return its;
}


int newton_raphson(const bool non_linear, const double strain[NVOI], const double *vars_old, Vec u, newton_t *newton)
{
	bc_apply_on_u(u, strain);

	int ierr, lits = 0;
	double lerr, lerr0;
	int solver_its;
	double solver_err;

	if (!newton) return 1;

	clean_newton(newton);

	while (lits < NEWTON_MAX_ITS) {

		lerr = assembly_res(b, u, vars_old);

		if (lits == 0)
			lerr0 = lerr;

		newton->err[lits] = lerr;

		if (lerr < NEWTON_MIN_TOL || lerr < lerr0 * NEWTON_REL_TOL)
			break;

		if (non_linear || lits > 0) {

			assembly_jac(A, u, vars_old);
			solver_its = solve(A, b, du, &solver_err);

		} else {

			solver_its = solve(A0, b, du, &solver_err);

		}

		newton->solver_its[lits] = solver_its;
		newton->solver_err[lits] = solver_err;

		ierr = VecAXPY(u, 1., du); CHKERRQ(ierr);

		lits++;

	}

	return lits;
}


void clean_newton(newton_t *newton)
{
	newton->its = 0;
	memset(newton->err, 0.0, NEWTON_MAX_ITS * sizeof(double));
	memset(newton->solver_its, 0, NEWTON_MAX_ITS * sizeof(int));
	memset(newton->solver_err, 0.0, NEWTON_MAX_ITS * sizeof(double));
	return;
}
