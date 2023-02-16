#include "build_cable_function_time_integration.h"

void BuildCableFunctionTimeIntegration::integrate1(const TriMatN1& phi, TriMatN1& eta_in_out, TriMatN1& psi_out) {
	for (int k = 0; k <= 3 * Params::N; k++) {
		psi_out[k] = phi[k] * Params::CM;
		double tmp = 0.5 * Params::DT * eta_in_out[k];
		eta_in_out[k] = psi_out[k] + tmp;
		psi_out[k] -= tmp;
	}
}

void BuildCableFunctionTimeIntegration::integrate2(const TriMatN1& phi, const TriMatN1& psi, const TriMatN1& eta, const VectorN1& g, TriMatN1& al_out, TriMatN1& ar_out) {
	double tmp = 0.5 * Params::DT * phi[0] * g[0];
	al_out[0] = psi[0] + tmp;
	ar_out[0] = eta[0] - tmp;
	tmp = 0.5 * Params::DT * phi[1] * g[1];
	al_out[1] = psi[1] + tmp;
	ar_out[1] = eta[1] - tmp;
	int k = 3;
	for (int i = 1; i < Params::N; i++) {
		tmp = 0.5 * Params::DT * phi[k - 1] * g[i - 1];
		al_out[k - 1] = psi[k - 1] + tmp;
		ar_out[k - 1] = eta[k - 1] - tmp;
		tmp = 0.5 * Params::DT * phi[k] * g[i];
		al_out[k] = psi[k] + tmp;
		ar_out[k] = eta[k] - tmp;
		tmp = 0.5 * Params::DT * phi[k + 1] * g[i + 1];
		al_out[k + 1] = psi[k + 1] + tmp;
		ar_out[k + 1] = eta[k + 1] - tmp;
		k += 3;
	}
	tmp = 0.5 * Params::DT * phi[k - 1] * g[Params::N - 1];
	al_out[k - 1] = psi[k - 1] + tmp;
	ar_out[k - 1] = eta[k - 1] - tmp;
	tmp = 0.5 * Params::DT * phi[k] * g[Params::N];
	al_out[k] = psi[k] + tmp;
	ar_out[k] = eta[k] - tmp;
}

void BuildCableFunctionTimeIntegration::integrate3(const TriMatN1& ar, const VectorN1& v, const TriMatN1& phi, const VectorN1& d, const VectorN1& x, double t, double t0, double& a_now, VectorN1& b_out) {
	TriMatOperations::triMatVecMultN1(ar, v, b_out);
	VectorN1 g{};
	TriMatOperations::triMatVecMultN1(phi, d, g);

	for (int i = 0; i < Params::N + 1; i++) {
		b_out[i] += Params::DT * g[i];
	}
	double a_next = Params::EX_INPUT_AMP;			//Injected current at node x = EX_INPUT_NODE at (timeStep+1)	[mu amp]
	if (t + 1 > t0 + Params::EX_INPUT_DURATION) {
		a_next = 0.0;
	}
	double tmp = 0.5 * Params::DT * (a_now + a_next);
	double ex_input_node = Params::EX_INPUT_NODE;
	if (ex_input_node > 0){
		b_out[ex_input_node - 1] += tmp * (x[ex_input_node] - x[ex_input_node - 1]) / 6;
	}
	b_out[ex_input_node]	 += tmp * (x[ex_input_node + 1] - x[ex_input_node - 1]) / 3;
	b_out[ex_input_node + 1] += tmp * (x[ex_input_node + 1] - x[ex_input_node])		/ 6;
	a_now = a_next;
}
