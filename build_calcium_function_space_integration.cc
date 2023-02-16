#include "build_calcium_function_space_integration.h"

using namespace std;

void BuildCalciumFunctionSpaceIntegration::LHS(const VectorN& segX, TriMatN& phi_out) {
	double tmp = (segX[1] - segX[0]) / 6.0;
	phi_out[0] = 2 * tmp;
	phi_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		phi_out[k - 1] = tmp;
		tmp = (segX[i + 1] - segX[i - 1]) / 3.0;
		phi_out[k] = tmp;
		tmp = (segX[i + 1] - segX[i]) / 6.0;
		phi_out[k + 1] = tmp;
		k += 3;
	}
	phi_out[k - 1] = tmp;
	phi_out[k] = 2 * tmp;
}

void BuildCalciumFunctionSpaceIntegration::RHS1(const VectorN& segX, const VectorN& a, const VectorN& v, TriMatN& chi_out) {
	double tmp1, tmp2 = (segX[1] - segX[0]) / 60.0;
	chi_out[0] = tmp2 * (12.0 * a[0] / v[0] + 3.0 * a[1] / v[0] + 3.0 * a[0] / v[1] + 2.0 * a[1] / v[1]);
	chi_out[1] = tmp2 * ( 3.0 * a[0] / v[0] + 2.0 * a[1] / v[0] + 2.0 * a[0] / v[1] + 3.0 * a[1] / v[1]);
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp1 = (segX[i] - segX[i - 1]) / 60.0;
		tmp2 = (segX[i + 1] - segX[i]) / 60.0;
		chi_out[k - 1] = tmp1 * ( 3.0 * a[i - 1] / v[i - 1] + 2.0 * a[i] / v[i - 1] + 2.0 * a[i - 1] / v[i] +  3.0 * a[i] / v[i]);
		chi_out[k]	   = tmp1 * ( 2.0 * a[i - 1] / v[i - 1] + 3.0 * a[i] / v[i - 1] + 3.0 * a[i - 1] / v[i] + 12.0 * a[i] / v[i]) +
						 tmp2 * (12.0 * a[i] / v[i] + 3.0 * a[i + 1] / v[i] + 3.0 * a[i] / v[i + 1] + 2.0 * a[i + 1] / v[i + 1]);
		chi_out[k + 1] = tmp2 * ( 3.0 * a[i] / v[i] + 2.0 * a[i + 1] / v[i] + 2.0 * a[i] / v[i + 1] + 3.0 * a[i + 1] / v[i + 1]);
		k += 3;
	}
	tmp1 = (segX[Params::N - 1] - segX[Params::N - 2]) / 60.0;
	chi_out[k - 1] = tmp1 * (3.0 * a[Params::N - 2] / v[Params::N - 2] + 2.0 * a[Params::N - 1] / v[Params::N - 2] + 2.0 * a[Params::N - 2] / v[Params::N - 1] +  3.0 * a[Params::N - 1] / v[Params::N - 1]);
	chi_out[k]	   = tmp1 * (2.0 * a[Params::N - 2] / v[Params::N - 2] + 3.0 * a[Params::N - 1] / v[Params::N - 2] + 3.0 * a[Params::N - 2] / v[Params::N - 1] + 12.0 * a[Params::N - 1] / v[Params::N - 1]);
}

void BuildCalciumFunctionSpaceIntegration::RHS2(const VectorN& segX, VectorN& rho_out) {
	rho_out[0] = (segX[1] - segX[0]) / 2.0;
	for (int k = 1; k < Params::N - 1; k++) {
		rho_out[k] = (segX[k + 1] - segX[k - 1]) / 2.0;
	}
	rho_out[Params::N - 1] = (segX[Params::N - 1] - segX[Params::N - 2]) / 2.0;
}

void BuildCalciumFunctionSpaceIntegration::RHS3_left(TriMatN& nu_left_out) {
	nu_left_out[0] = 0.0;
	nu_left_out[1] = 0.0;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		nu_left_out[k - 1] = - 0.5;
		nu_left_out[k]     =   0.5;
		nu_left_out[k + 1] =   0.0;
		k += 3;
	}
	nu_left_out[k - 1] = - 0.5;
	nu_left_out[k]     =   0.5;
}

void BuildCalciumFunctionSpaceIntegration::RHS3_right(TriMatN& nu_right_out) {
	nu_right_out[0] = - 0.5;
	nu_right_out[1] =   0.5;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		nu_right_out[k - 1] =   0.0;
		nu_right_out[k]     = - 0.5;
		nu_right_out[k + 1] =   0.5;
		k += 3;
	}
	nu_right_out[k - 1] = 0.0;
	nu_right_out[k]     = 0.0;
}

void BuildCalciumFunctionSpaceIntegration::RHS4_left_minus_right(const VectorN& a_l, const VectorN& a_r, const VectorN& v, TriMatN& ksi_out) {
	double tmp1, tmp2 = 1.0 / 12.0 * (3.0 * a_r[0] / v[0] + a_r[1] / v[0] + a_r[0] / v[1] + a_r[1] / v[1]);
	ksi_out[0] =   tmp2;
	ksi_out[1] = - tmp2;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp1 = 1.0 / 12.0 * (a_l[i - 1] / v[i - 1] + a_l[i] / v[i - 1] + a_l[i - 1] / v[i] + 3.0 * a_l[i] / v[i]);
		tmp2 = 1.0 / 12.0 * (3.0 * a_r[i] / v[i] + a_r[i + 1] / v[i] + a_r[i] / v[i + 1] + a_r[i + 1] / v[i + 1]);
		ksi_out[k - 1] = - tmp1;
		ksi_out[k]     =   tmp1 + tmp2;
		ksi_out[k + 1] = - tmp2;
		k += 3;
	}
	tmp1 = 1.0 / 12.0 * (a_l[Params::N - 2] / v[Params::N - 2] + a_l[Params::N - 1] / v[Params::N - 2] + a_l[Params::N - 2] / v[Params::N - 1] + 3.0 * a_l[Params::N - 1] / v[Params::N - 1]);
	ksi_out[k - 1] = - tmp1;
	ksi_out[k]     =   tmp1;
}

void BuildCalciumFunctionSpaceIntegration::RHS5(const VectorN& segX, const VectorN& a, const VectorN& v, TriMatN& mu_shell_out) {
	double tmp1, tmp2 = (segX[1] - segX[0]) / 60.0;
	mu_shell_out[0] = tmp2 * (12.0 * a[0] / v[0] + 3.0 * a[1] / v[0] + 3.0 * a[0] / v[1] + 2.0 * a[1] / v[1]);
	mu_shell_out[1] = tmp2 * ( 3.0 * a[0] / v[0] + 2.0 * a[1] / v[0] + 2.0 * a[0] / v[1] + 3.0 * a[1] / v[1]);
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp1 = (segX[i] - segX[i - 1]) / 60.0;
		tmp2 = (segX[i + 1] - segX[i]) / 60.0;
		mu_shell_out[k - 1] = tmp1 * ( 3.0 * a[i - 1] / v[i - 1] + 2.0 * a[i] / v[i - 1] + 2.0 * a[i - 1] / v[i] +  3.0 * a[i] / v[i]);
		mu_shell_out[k]     = tmp1 * ( 2.0 * a[i - 1] / v[i - 1] + 3.0 * a[i] / v[i - 1] + 3.0 * a[i - 1] / v[i] + 12.0 * a[i] / v[i]) +
							  tmp2 * (12.0 * a[i] / v[i] + 3.0 * a[i + 1] / v[i] + 3.0 * a[i] / v[i + 1] + 2.0 * a[i + 1] / v[i + 1]);
		mu_shell_out[k + 1] = tmp2 * ( 3.0 * a[i] / v[i] + 2.0 * a[i + 1] / v[i] + 2.0 * a[i] / v[i + 1] + 3.0 * a[i + 1] / v[i + 1]);
		k += 3;
	}
	tmp1 = (segX[Params::N - 1] - segX[Params::N - 2]) / 60.0;
	mu_shell_out[k - 1] = tmp1 * (3.0 * a[Params::N - 2] / v[Params::N - 2] + 2.0 * a[Params::N - 1] / v[Params::N - 2] + 2.0 * a[Params::N - 2] / v[Params::N - 1] +  3.0 * a[Params::N - 1] / v[Params::N - 1]);
	mu_shell_out[k]     = tmp1 * (2.0 * a[Params::N - 2] / v[Params::N - 2] + 3.0 * a[Params::N - 1] / v[Params::N - 2] + 3.0 * a[Params::N - 2] / v[Params::N - 1] + 12.0 * a[Params::N - 1] / v[Params::N - 1]);
}

void BuildCalciumFunctionSpaceIntegration::RHS6(const VectorN& segX, const VectorN& a, TriMatN& mu_cell_out) {
	double tmp1, tmp2 = (segX[1] - segX[0]) / 12.0;
	mu_cell_out[0] = tmp2 * (3.0 * a[0] + a[1]);
	mu_cell_out[1] = tmp2 * (	   a[0] + a[1]);
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp1 = (segX[i] - segX[i - 1]) / 12.0;
		tmp2 = (segX[i + 1] - segX[i]) / 12.0;
		mu_cell_out[k - 1] = tmp1 * (a[i - 1] + a[i]);
		mu_cell_out[k]     = tmp1 * a[i - 1] + 3 * (tmp1 + tmp2) * a[i] + tmp2 * a[i + 1];
		mu_cell_out[k + 1] = tmp2 * (a[i] + a[i + 1]);
		k += 3;
	}
	tmp1 = (segX[Params::N - 1] - segX[Params::N - 2]) / 12.0;
	mu_cell_out[k - 1] = tmp1 * (a[Params::N - 2] +		  a[Params::N - 1]);
	mu_cell_out[k]     = tmp1 * (a[Params::N - 2] + 3.0 * a[Params::N - 1]);
}
