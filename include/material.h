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

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <stdbool.h>

typedef struct {

	double E, nu, Ka, Sy;
	double k, mu, lambda;
	int type;
	bool plasticity, damage;


} material_t;

// material.c
void material_set(material_t *mat, double _E, double _nu, double _Ka, double _Sy, int _type);

void get_dev_tensor(const double tensor[6], double tensor_dev[6]);

bool plastic_law(const material_t *mat, const double eps[6], const double eps_p_old[6], const double alpha_old,
		 double *_dl, double _normal[6], double _s_trial[6], double *_f_trial);
void plastic_get_stress(const material_t *mat, const double eps[6], const double eps_p_old[6],
			const double alpha_old, double stress[6]);
void plastic_get_ctan(const material_t *mat, const double eps[6], const double eps_p_old[6],
		      const double alpha_old, double ctan[6][6]);
bool plastic_evolute(const material_t *mat, const double eps[6], const double eps_p_old[6], const double alpha_old,
		     double *eps_p_new, double *alpha_new, double *f_trial);
void isolin_get_ctan(const material_t *material, double ctan[6][6]);
void isolin_get_stress(const material_t *material, const double eps[6], double stress[6]);

void material_print(void);

#endif
