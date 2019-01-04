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


int solve(Mat _A, Vec _b, Vec _x, int *_its, double *_err)
{
	return solve_v(_A, _b, _x, PCJACOBI, KSPCG, _its, _err);
}


int solve_v(Mat _A, Vec _b, Vec _x, PCType _PC, KSPType _KSP,
	    int *_its, double *_err)
{
	MICROC_INST_START

	int ierr, its;
	double norm;

	ierr = KSPSetOperators(ksp, _A, _A); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, _KSP); CHKERRQ(ierr);
	ierr = PCSetType(pc, _PC); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	ierr = KSPSolve(ksp, _b, _x); CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp, &its);
	ierr = KSPGetResidualNorm(ksp, &norm);

	*_err = norm;
	*_its = its;

	MICROC_INST_END
	return its;
}


int solve_ts(Mat _A, Vec _b, Vec _x, PCType _PC, KSPType _KSP,
	     int *_its, double *_err)
{
	/* Thread safe routine for OMP */

	MICROC_INST_START

	int ierr, its;
	double norm;

	KSP ksp;
	PC pc;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp, _A, _A); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, _KSP); CHKERRQ(ierr);
	ierr = PCSetType(pc, _PC); CHKERRQ(ierr);
	//ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	ierr = KSPSolve(ksp, _b, _x); CHKERRQ(ierr);

	ierr = KSPGetIterationNumber(ksp, &its);
	ierr = KSPGetResidualNorm(ksp, &norm);

	*_err = norm;
	*_its = its;

	MICROC_INST_END
	return its;
}


//int newton_raphson(const bool non_linear, const double strain[NVOI],
//		   const double *vars_old, Vec u,
//		   newton_t *newton)
//{
//	return newton_raphson_v(non_linear, strain, vars_old, u,
//				PCJACOBI, KSPCG,
//				newton);
//}


int newton_raphson_v(Mat _A,
		     Mat _A0,
		     Vec _b,
		     Vec _u,
		     Vec _du,
		     const bool non_linear,
		     const double strain[NVOI],
		     const double *_vars_old,
		     PCType _PC, KSPType _KSP,
		     newton_t *newton,
		     bool print)
{
	MICROC_INST_START

	bc_apply_on_u(_u, strain);

	int ierr, lits = 0;
	double lerr, lerr0;
	int solver_its;
	double solver_err;

	if (!newton) return 1;

	clean_newton(newton);

	while (lits < NEWTON_MAX_ITS) {

		lerr = assembly_res(_b, _u, _vars_old);

		if (lits == 0)
			lerr0 = lerr;

		newton->err[lits] = lerr;
		if (print) {
			int thread_id = omp_get_thread_num();
			printf("Thread : %d |b| = %e\n", thread_id, lerr);
		}

		if (lerr < NEWTON_MIN_TOL || lerr < lerr0 * NEWTON_REL_TOL)
			break;

		if (non_linear || lits > 0) {

			assembly_jac(_A, _u, _vars_old);
			solve_ts(_A, _b, _du, _PC, _KSP, &solver_its, &solver_err);

		} else {

			solve_ts(_A0, _b, _du, _PC, _KSP, &solver_its, &solver_err);

		}

		newton->solver_its[lits] = solver_its;
		newton->solver_err[lits] = solver_err;

		ierr = VecAXPY(_u, 1., _du); CHKERRQ(ierr);

		lits++;

	}

	MICROC_INST_END
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
