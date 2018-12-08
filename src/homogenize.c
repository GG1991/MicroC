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

#define D_EPS 1.0e-5


void microc_set_macro_strain(const int gp_id, const double *macro_strain)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(gp_list[gp_id].macro_strain, macro_strain, NVOI * sizeof(double));
}


void microc_get_macro_strain(const int gp_id, double *macro_strain)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(macro_strain, gp_list[gp_id].macro_strain, NVOI * sizeof(double));
}


void microc_get_macro_stress(const int gp_id, double *macro_stress)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(macro_stress, gp_list[gp_id].macro_stress, NVOI * sizeof(double));
}


void microc_get_macro_ctan(const int gp_id, double *macro_ctan)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);
	memcpy(macro_ctan, gp_list[gp_id].macro_ctan, NVOI * NVOI * sizeof(double));
}


void microc_homogenize()
{
	int newton_its, solver_its[NEWTON_MAX_ITS];
	double newton_err[NEWTON_MAX_ITS], solver_err[NEWTON_MAX_ITS];
	bool nl_flag;

	int igp;
	for (igp = 0; igp < ngp; ++igp) {

		gp_t *gp_ptr = &gp_list[igp];

		gp_ptr->sigma_cost = 0;

		double *vars_new;
		if (gp_ptr->allocated) {
			vars_new = gp_ptr->int_vars_k;
		} else {
//			vars_new = vars_new_aux;
		}

		// SIGMA 1 Newton-Raphson
		memcpy(gp_ptr->u_k, gp_ptr->u_n, nndim * sizeof(double));

//		newton_its = newton_raphson(gp_ptr->allocated, gp_ptr->macro_strain,
//					    gp_ptr->int_vars_n, gp_ptr->u_k,
//					    newton_err, solver_its, solver_err);

		gp_ptr->sigma_newton_its = newton_its;
		memcpy(gp_ptr->sigma_newton_err, newton_err, NEWTON_MAX_ITS * sizeof(double));
		memcpy(gp_ptr->sigma_solver_its, solver_its, NEWTON_MAX_ITS * sizeof(int));
		memcpy(gp_ptr->sigma_solver_err, solver_err, NEWTON_MAX_ITS * sizeof(double));

		for (int i = 0; i < newton_its; ++i)
			gp_ptr->sigma_cost += solver_its[i];

//		calc_ave_stress(gp_ptr->u_k, gp_ptr->int_vars_n, gp_ptr->macro_stress);


		/* Updates <vars_new> and <f_trial_max> */
//		nl_flag = calc_vars_new(gp_ptr->u_k, gp_ptr->int_vars_n, vars_new, &f_trial_max);

		if (nl_flag == true) {
			if (gp_ptr->allocated == false) {
//				gp_ptr->allocate(num_int_vars);
				memcpy(gp_ptr->int_vars_k, vars_new, nvars * sizeof(double));
			}
		}

		if (gp_ptr->allocated) {

			// CTAN 3/6 Newton-Raphsons in 2D/3D
			double eps_1[6], sig_0[6], sig_1[6];

//			memcpy(u_aux, gp_ptr->u_k, nndim * sizeof(double));
			memcpy(sig_0, gp_ptr->macro_stress, NVOI * sizeof(double));

			for (int i = 0; i < NVOI; ++i) {

				memcpy(eps_1, gp_ptr->macro_strain, NVOI * sizeof(double));
				eps_1[i] += D_EPS;

//				newton_its = newton_raphson(true, eps_1, gp_ptr->int_vars_n,
//							    u_aux,
//							    newton_err, solver_its, solver_err);

				for (int i = 0; i < newton_its; ++i)
					gp_ptr->sigma_cost += solver_its[i];

//				calc_ave_stress(u_aux, gp_ptr->int_vars_n, sig_1);

				for (int v = 0; v < NVOI; ++v)
					gp_ptr->macro_ctan[v * NVOI + i] =
						(sig_1[v] - sig_0[v]) / D_EPS;

			}
		}
	}
}


//void update_vars()
//{
//
//	for (igp = 0; igp < ngp; ++igp)
//		gp_list[igp].update_vars();
//}
