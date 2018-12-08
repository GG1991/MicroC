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

#include "material.h"


void material_set(material_t *mat, double _E, double _nu, double _Ka, double _Sy, int _type)
{
	mat->E = _E;
	mat->nu = _nu;
	mat->Ka = _Ka;
	mat->Sy = _Sy;
	mat->type = _type;

	mat->k = _E / (3. * (1. - 2. * _nu));
	mat->mu = _E / (2. * (1. + _nu));
	mat->lambda = _nu * _E / ((1. + _nu) * (1. - 2. * _nu));

	if (_type == 0) {
		// lineal
		mat->plasticity = false;
		mat->damage = false;
	} else if (_type == 1) {
		// con plasticidad
		mat->plasticity = true;
		mat->damage = false;
	} else if (_type == 2) {
		// con daño
		mat->plasticity = false;
		mat->damage = true;
	}
}


void material_print(material_t *mat)
{
	printf("E = %e, nu = %e, Ka = %e, Sy = %e, type = %1d\n",
	       mat->E, mat->nu , mat->Ka, mat->Sy, mat->type);
}

void get_dev_tensor(const double tensor[6], double tensor_dev[6])
{
	memcpy(tensor_dev, tensor, 6 * sizeof(double));

	int i;
	for (i = 0; i < 3; i++)
		tensor_dev[i] -= 0.333333333333 * (tensor[0] + tensor[1] + tensor[2]);
}


bool plastic_law(const material_t *mat, const double eps[6], const double eps_p_old[6], const double alpha_old,
		 double *_dl, double _normal[6], double _s_trial[6], double *_f_trial)
{
	/* Calculates _dl, _normal and _s_trial
	 * Optionally  returns f_trial value
	 */

	int i;
	double eps_dev[6], eps_p_dev_1[6];

	get_dev_tensor(eps_p_old, eps_p_dev_1);
	get_dev_tensor(eps, eps_dev);

	for (i = 0; i < 3; ++i)
		_s_trial[i] = 2 * mat->mu * (eps_dev[i] - eps_p_dev_1[i]);

	for (i = 3; i < 6; ++i)
		_s_trial[i] = mat->mu * (eps_dev[i] - eps_p_dev_1[i]);

	double tmp = 0.0;
	for (i = 0; i < 6; ++i)
		tmp += _s_trial[i] * _s_trial[i];
	double s_norm = sqrt(tmp);

	double f_trial = s_norm - SQRT_2DIV3 * (mat->Sy + mat->Ka * alpha_old);

	if (f_trial > 0 && mat->plasticity) {

		for (i = 0; i < 6; ++i)
			_normal[i] = _s_trial[i] / s_norm;
		*_dl = f_trial / (2. * mat->mu * (1. + mat->Ka / (3. * mat->mu)));

		if (_f_trial != NULL)
			*_f_trial = 0.;

		return true;

	} else {

		memset(_normal, 0, 6 * sizeof(double));
		*_dl = 0;

		if (_f_trial != NULL)
			*_f_trial = f_trial;

		return false;
	}
}


void plastic_get_stress(const material_t *mat, const double eps[6], const double eps_p_old[6],
			const double alpha_old, double stress[6])
{
	int i;
	double dl, normal[6], s_trial[6];
	bool nl_flag = plastic_law(mat, eps, eps_p_old, alpha_old, &dl, normal, s_trial, NULL);

	//sig_2 = s_trial + K * tr(eps) * 1 - 2 * mu * dl * normal;
	memcpy(stress, s_trial, 6 * sizeof(double));

	for (i = 0; i < 3; ++i)
		stress[i] += mat->k * (eps[0] + eps[1] + eps[2]);

	for (i = 0; i < 6; ++i)
		stress[i] -= 2 * mat->mu * dl * normal[i];
}


void plastic_get_ctan(const material_t *mat, const double eps[6], const double eps_p_old[6],
		      const double alpha_old, double ctan[6][6])
{
	int i, j;
	double stress_0[6];
	plastic_get_stress(mat, eps, eps_p_old, alpha_old, stress_0);

	for (i = 0; i < 6; ++i) {

		double eps_1[6];
		memcpy(eps_1, eps, 6 * sizeof(double));
		eps_1[i] += D_EPS_CTAN;

		double stress_1[6];
		plastic_get_stress(mat, eps_1, eps_p_old, alpha_old, stress_1);

		for (j = 0; j < 6; ++j)
			ctan[j][i] = (stress_1[j] - stress_0[j]) / D_EPS_CTAN;
	}
}


bool plastic_evolute(const material_t *mat, const double eps[6], const double eps_p_old[6], const double alpha_old,
		     double *eps_p_new, double *alpha_new, double *f_trial)
{
	double dl, normal[6], s_trial[6];
	bool nl_flag = plastic_law(mat, eps, eps_p_old, alpha_old, &dl, normal, s_trial, f_trial);

	int i;
	for (i = 0; i < 6; ++i)
		eps_p_new[i] = eps_p_old[i] + dl * normal[i];

	*alpha_new = alpha_old + SQRT_2DIV3 * dl;

	return nl_flag;
}

void isolin_get_ctan(const material_t *material, double ctan[6][6])
{
	// C = lambda * (1x1) + 2 mu I
	memset(ctan, 0, 6 * 6 * sizeof(double));

	int i;
	for (i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ctan[i][j] += material->lambda;

	for (i = 0; i < 3; ++i)
		ctan[i][i] += 2 * material->mu;

	for (i = 3; i < 6; ++i)
		ctan[i][i] = material->mu;
}


void isolin_get_stress(const material_t *material, const double eps[6], double stress[6])
{
	int i;
	for (i = 0; i < 3; ++i)
		stress[i] = material->lambda * (eps[0] + eps[1] + eps[2]) \
			    + 2 * material->mu * eps[i];

	for (i = 3; i < 6; ++i)
		stress[i] = material->mu * eps[i];
}
