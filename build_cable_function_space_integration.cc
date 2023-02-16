#include "build_cable_function_space_integration.h"

void BuildCableFunctionSpaceIntegration::LHS(const VectorN1& x, const VectorN1& p, TriMatN1& phi_out) {
	phi_out[0] = (x[1] - x[0]) * (3.0 * p[0] + p[1]) / 12.0;
	double tmp = (x[1] - x[0]) * (p[0] + p[1]) / 12.0;
	phi_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N; i++) {
		phi_out[k - 1] = tmp;
		tmp = ((x[i] - x[i - 1]) * (p[i - 1] + 3.0 * p[i]) + (x[i + 1] - x[i]) * (3.0 * p[i] + p[i + 1])) / 12.0;
		phi_out[k] = tmp;
		tmp = (x[i + 1] - x[i]) * (p[i] + p[i + 1]) / 12.0;
		phi_out[k + 1] = tmp;
		k += 3;
	}
	phi_out[k - 1] = tmp;
	phi_out[k] = (x[Params::N] - x[Params::N - 1]) * (p[Params::N - 1] + 3.0 * p[Params::N]) / 12.0;
}

void BuildCableFunctionSpaceIntegration::RHS(const VectorN1& x, const VectorN1& a, TriMatN1& eta_out) {
	double tmp = 0.5 * Params::GA * (a[0] + a[1]) / (x[1] - x[0]);
	eta_out[0] = -tmp;
	eta_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N; i++) {
		eta_out[k - 1] = tmp;
		eta_out[k] = -tmp;
		tmp = 0.5 * Params::GA * (a[i] + a[i + 1]) / (x[i + 1] - x[i]);
		eta_out[k] -= tmp;
		eta_out[k + 1] = tmp;
		k += 3;
	}
	eta_out[k - 1] = tmp;
	eta_out[k] = -tmp;
}
