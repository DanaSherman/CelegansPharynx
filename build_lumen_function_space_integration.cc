#include "build_lumen_function_space_integration.h"

void BuildLumenFunctionSpaceIntegration::LHS(const VectorN& segX, const VectorN& C, TriMatN& mat_out) {
	mat_out[0] = (3.0 * C[0] + C[1]) * (segX[1] - segX[0]) / 12.0;
	double tmp = (      C[0] + C[1]) * (segX[1] - segX[0]) / 12.0;
	mat_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		mat_out[k - 1] = tmp;
		tmp = (C[i - 1] * (segX[i] - segX[i - 1]) + 3.0 * C[i] * (segX[i + 1] - segX[i - 1]) + C[i + 1] * (segX[i + 1] - segX[i])) / 12.0;
		mat_out[k] = tmp;
		tmp = (C[i] + C[i + 1]) * (segX[i + 1] - segX[i]) / 12.0;
		mat_out[k + 1] = tmp;
		k += 3;
	}
	mat_out[k - 1] = tmp;
	mat_out[k] = (C[Params::N - 2] + 3 * C[Params::N - 1]) * (segX[Params::N - 1] - segX[Params::N - 2]) / 12.0;
}

void BuildLumenFunctionSpaceIntegration::RHS1(const VectorN& segX, TriMatN& phi_out) {
	double tmp = (segX[1] - segX[0]) / 6.0;
	phi_out[0] = 2 * tmp;
	phi_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		phi_out[k - 1] = tmp;
		tmp = (segX[i + 1] - segX[i - 1]) / 3.0;
		phi_out[k]     = tmp;
		tmp = (segX[i + 1] - segX[i]) / 6.0;
		phi_out[k + 1] = tmp;
		k += 3;
	}
	phi_out[k - 1] = tmp;
	phi_out[k]     = 2 * tmp;
}

void BuildLumenFunctionSpaceIntegration::RHS2(const VectorN& segX, const VectorN& C, const VectorN& R, VectorN& vec_out) {
	double tmp1;
	double tmp2 = (3.0 * C[0] + C[1]) * (segX[1] - segX[0]) / 12.0;
	double tmp3 = (      C[0] + C[1]) * (segX[1] - segX[0]) / 12.0;
	vec_out[0] = tmp2 * R[0] + tmp3 * R[1];
	for (int i = 1; i < Params::N - 1; i++) {
		tmp1 = (      C[i - 1] + C[i]) * (segX[i] - segX[i - 1]) / 12.0;
		tmp2 = (C[i - 1] * (segX[i] - segX[i - 1]) + 3.0 * C[i] * (segX[i + 1] - segX[i - 1]) + C[i + 1] * (segX[i + 1] - segX[i])) / 12.0;
		tmp3 = (C[i] + C[i + 1]) * (segX[i + 1] - segX[i]) / 12.0;
		vec_out[i] = tmp1 * R[i - 1] + tmp2 * R[i] + tmp3 * R[i + 1];
	}
	tmp1 = (C[Params::N - 2] +		 C[Params::N - 1]) * (segX[Params::N - 1] - segX[Params::N - 2]) / 12.0;
	tmp2 = (C[Params::N - 2] + 3.0 * C[Params::N - 1]) * (segX[Params::N - 1] - segX[Params::N - 2]) / 12.0;
	vec_out[Params::N - 1] = tmp1 * R[Params::N - 2] + tmp2 * R[Params::N - 1];
}
