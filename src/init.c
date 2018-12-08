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


PetscErrorCode microc_init()
{

	PetscErrorCode ierr;
	ierr = PetscInitialize(NULL, NULL, (char*)0, help);

       	if(ierr)
	       	return ierr;

	PetscPrintf(PETSC_COMM_WORLD,
		    "\nMicroC : A HPC for FE2 Multi-scale Simulations\n\n");

	lx = LX;
	ly = LY;
	lz = LZ;

	vtu_freq = VTU_FREQ;
	newton_max_its = NEWTON_MAX_ITS;
	newton_min_tol = NEWTON_MIN_TOL;

	PetscOptionsGetReal(NULL, NULL, "-lx", &lx, NULL);
	PetscOptionsGetReal(NULL, NULL, "-ly", &ly, NULL);
	PetscOptionsGetReal(NULL, NULL, "-lz", &lz, NULL);
	PetscOptionsGetReal(NULL, NULL, "-new_tol", &newton_min_tol, NULL);
	PetscOptionsGetInt(NULL, NULL, "-vtu_freq", &vtu_freq, NULL);
	PetscOptionsGetInt(NULL, NULL, "-new_its", &newton_max_its, NULL);

	DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE,
		       bz = DM_BOUNDARY_NONE;
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX,
			    NX, NY, NZ,
			    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			    DIM, 1, NULL, NULL, NULL, &da);

	ierr = DMSetMatType(da, MATAIJ); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, &A); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &u); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &b); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &du); CHKERRQ(ierr);

	ierr = VecSetOption(u, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
	ierr = VecSetOption(b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

	ierr = VecZeroEntries(u); CHKERRQ(ierr);
	ierr = VecZeroEntries(b); CHKERRQ(ierr);
	ierr = VecZeroEntries(du); CHKERRQ(ierr);

	PetscInt M, N, P;
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0,
			   0, 0, 0, 0, 0); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of Elements : %ld\n", (M - 1) * (N - 1) * (P - 1));
	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of Nodes    : %ld\n", M * N * P);
	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of DOFs     : %ld\n\n", (M * N * P) * DIM);

	dx = lx / (M - 1);
	dy = ly / (N - 1);
	dz = lz / (P - 1);
	wg = dx * dy * dz / NPE;

	KSPType ksptype;
	PetscReal rtol, abstol, dtol;
	PetscInt maxits;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
	ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	ierr = KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxits);
	ierr = KSPGetType(ksp, &ksptype);
	PetscPrintf(PETSC_COMM_WORLD,
		    "KSP Info: type = %s\trtol = %e\t\
		    abstol = %e\tdtol = %e\tmaxits = %d\n",
		    ksptype, rtol, abstol, dtol, maxits);

	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);

	return ierr;
}


PetscErrorCode finish()
{
	PetscErrorCode ierr;
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&u); CHKERRQ(ierr);
	ierr = VecDestroy(&b); CHKERRQ(ierr);
	ierr = VecDestroy(&du); CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

	ierr = PetscFinalize();
	return ierr;
}
