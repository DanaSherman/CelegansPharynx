#include "build_lumen_function_time_integration.h"

void BuildLumenFunctionTimeIntegration::integrate1(const TriMatN& ksi, TriMatN& eta_in_out, TriMatN& psi_out) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	for (int k = 0; k <= 3 * (Params::N - 1); k++) {
		double tmp = 0.5 * DT_SEC * ksi[k];
		psi_out[k]    = - eta_in_out[k] - tmp;
		eta_in_out[k] = - eta_in_out[k] + tmp;
	}
}

void BuildLumenFunctionTimeIntegration::integrate2(const TriMatN& psi, const TriMatN& eta, TriMatN& al_out, TriMatN& ar_out) {
	for (int k = 0; k <= 3 * (Params::N - 1); k++) {
		al_out[k] = psi[k];
		ar_out[k] = eta[k];
	}
}

void BuildLumenFunctionTimeIntegration::integrate3(const TriMatN& ar, const VectorN& lCurr, const TriMatN& chi, const TriMatN& phi, const VectorN& F_CE05, const VectorN& F_OE05, const VectorN& rho, const VectorN& vCurr, const VectorN& vNext, VectorN& b_out) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	TriMatOperations::triMatVecMultN(ar, lCurr, b_out);
	VectorN v1_0{}, v_vec{}, F_CE_vec{}, F_OE_vec{};
	for (int i = 0; i < Params::N; i++) {
		v1_0[i] = vNext[i] - vCurr[i];														//v(timeStep + 1) - v(timeStep)
	}
	TriMatOperations::triMatVecMultN(chi, v1_0, v_vec);
	TriMatOperations::triMatVecMultN(phi, F_CE05, F_CE_vec);
	TriMatOperations::triMatVecMultN(phi, F_OE05, F_OE_vec);

	for (int i = 0; i < Params::N; i++) {
		b_out[i] += DT_SEC * (2 * F_CE_vec[i] + F_OE_vec[i] - rho[i]) - v_vec[i];		//[g*cm^2/s]
	}
}
