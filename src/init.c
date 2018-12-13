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


int microc_init(const int _ngp, const int _size[3], const int _type,
		const double *_params, const material_t *_materials)
{
	return microc_initv(_ngp, _size, _type, _params, _materials, 0, NULL);
}


int microc_initv(const int _ngp, const int _size[3], const int _type,
		const double *_params, const material_t *_materials,
		int argc, char **argv)
{
	MICROC_INST_START

	int ierr;
	if (!first_init) {
		// This is done to call init multiple time without causing
		// errors
		ierr = MPI_Init(NULL, NULL);
		first_init = 1;
	}

	ierr = PetscInitialize(&argc, &argv, (char*)0, help);

	if(ierr)
		return ierr;

	printf("\nMicroC : A HPC for FE2 Multi-scale Simulations\n\n");

	int m;
	for (m = 0; m < NMATERIALS; ++m)
		memcpy(&material_list[m], &_materials[m], sizeof(material_t));

	lx = LX;
	ly = LY;
	lz = LZ;

	nx = _size[0];
	ny = _size[1];
	nz = _size[2];

	newton_max_its = NEWTON_MAX_ITS;
	newton_min_tol = NEWTON_MIN_TOL;

	ierr = DMDACreate3d(PETSC_COMM_WORLD,
			    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
			    DMDA_STENCIL_BOX,
			    nx, ny, nz,
			    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			    DIM, 1, NULL, NULL, NULL, &da);

	ierr = DMSetMatType(da, MATAIJ); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, &A); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, &A0); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &u); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &b); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &du); CHKERRQ(ierr);

	int npe;
	ierr = DMDAGetElements(da, &nelem, &npe, &eix); CHKERRQ(ierr);
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);
	assert(npe == NPE);

	int ex, ey, ez;
	elem_type = malloc(nelem * sizeof(int));
	for (ex = 0; ex < nex; ++ex)
		for (ey = 0; ey < ney; ++ey)
			for (ez = 0; ez < nez; ++ez)
				elem_type[elm_index(ex,ey,ez)] =  get_elem_type(ex, ey, ez);


	PetscInt M, N, P;
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	nn = M * N * P;
	nndim = nn * DIM;

	printf("Number of Elements : %ld\n", nelem);
	printf("Number of Nodes    : %ld\n", M * N * P);
	printf("Number of DOFs     : %ld\n\n", (M * N * P) * DIM);

	dx = lx / (M - 1);
	dy = ly / (N - 1);
	dz = lz / (P - 1);
	wg = dx * dy * dz / NPE;

	KSPType ksptype;
	PetscReal rtol, abstol, dtol;
	PetscInt maxits;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
	ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);


	// Init <struct gp_t> list
	ngp = _ngp;
        gp_list = malloc(ngp * sizeof(gp_t));

	int gp;
	for (gp = 0; gp < ngp; ++gp) {
		gp_list[gp].u_n = (double *) calloc(nndim, sizeof(double));
		gp_list[gp].u_k = (double *) calloc(nndim, sizeof(double));
	}

	// Set BC structures
	nbcs = ((nx * nz) + (nx * ny) + (ny * nz)) * 2 * DIM;
	bc_index = malloc(nbcs * sizeof(int));
	bc_value = malloc(nbcs * sizeof(double));

	int i, j, k, d;
	int count = 0;

	/* Y = 0 */
	for (i = 0; i < nx; ++i) {
		for (k = 0; k < nz; ++k) {
			const int ix = nod_index(i, 0, k);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}

	/* Y = LY */
	for (i = 0; i < nx; ++i) {
		for (k = 0; k < nz; ++k) {
			const int ix = nod_index(i, ny - 1, k);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}

	/* Z = 0 */
	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			const int ix = nod_index(i, j, 0);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}

	/* Z = LZ */
	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			const int ix = nod_index(i, j, nz - 1);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}
				
	/* X = 0 */
	for (j = 0; j < ny; ++j) {
		for (k = 0; k < nz; ++k) {
			const int ix = nod_index(0, j, k);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}

	/* X = LX */
	for (j = 0; j < ny; ++j) {
		for (k = 0; k < nz; ++k) {
			const int ix = nod_index(nx - 1, j, k);
			for (d = 0; d < DIM; ++d)
				bc_index[count++] = ix * DIM + d;
		}
	}

	ierr = VecZeroEntries(u);
	ierr = assembly_jac(A0, u, NULL);

	MICROC_INST_END
	return ierr;
}


int microc_finish(void)
{
	int ierr;
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = MatDestroy(&A0); CHKERRQ(ierr);
	ierr = VecDestroy(&u); CHKERRQ(ierr);
	ierr = VecDestroy(&b); CHKERRQ(ierr);
	ierr = VecDestroy(&du); CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
	ierr = DMDestroy(&da); CHKERRQ(ierr);

	free(bc_index);
	free(bc_value);


	int gp;
	for (gp = 0; gp < ngp; ++gp) {
		free(gp_list[gp].u_n);
		free(gp_list[gp].u_k);
	}
	free(gp_list);
	free(elem_type);

	ierr = PetscFinalize(); CHKERRQ(ierr);
	return ierr;
}
