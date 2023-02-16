#include "build_axial_current_function.h"

void BuildAxialCurrentFunction::buildXi(const VectorN1& x, TriMatN1& xi_out) {
	xi_out[0] = (x[1] - x[0]) / 3.0;
	double tmp = (x[1] - x[0]) / 6.0;
	xi_out[1] = tmp;
	int k = 3;
	for (int i = 1; i < Params::N; i++) {
		xi_out[k - 1] = tmp;
		tmp = (x[i + 1] - x[i - 1]) / 3.0;
		xi_out[k] = tmp;
		tmp = (x[i + 1] - x[i]) / 6.0;
		xi_out[k + 1] = tmp;
		k += 3;
	}
	xi_out[k - 1] = tmp;
	xi_out[k] = (x[Params::N] - x[Params::N - 1]) / 3.0;
}

void BuildAxialCurrentFunction::buildChi(const VectorN1& a, TriMatN1& chi_out) {
	double tmp = Params::GA * (2.0 * a[0] + a[1]) / 6.0;
	chi_out[0] = tmp;
	chi_out[1] = -tmp;
	int k = 3;
	for (int i = 1; i < Params::N; i++) {
		tmp = Params::GA * (a[i - 1] + 2.0 * a[i]) / 6.0;
		chi_out[k - 1] = tmp;
		chi_out[k] = -tmp;
		tmp = Params::GA * (2.0 * a[i] + a[i + 1]) / 6.0;

		chi_out[k] += tmp;
		chi_out[k + 1] = -tmp;
		k += 3;
	}
	tmp = Params::GA * (a[Params::N - 1] + 2.0 * a[Params::N]) / 6.0;
	chi_out[k - 1] = tmp;
	chi_out[k] = -tmp;
}
