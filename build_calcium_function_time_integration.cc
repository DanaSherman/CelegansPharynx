#include "build_calcium_function_time_integration.h"

void BuildCalciumFunctionTimeIntegration::integrate1(const Params& params, const VectorN& tauCaSec, const TriMatN& phi, const TriMatN& ksi_shell, TriMatN& eta_out, TriMatN& psi_out, int shell_terms) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	constexpr double D_CA = Params::CA_DIFFUSE_COEF;
	double tmp;
	tmp = 0.5 * DT_SEC * (1 / tauCaSec[0] * phi[0] + shell_terms * D_CA * ksi_shell[0]);
	psi_out[0] = phi[0] + tmp;
	eta_out[0] = phi[0] - tmp;
	tmp = 0.5 * DT_SEC * (1 / tauCaSec[1] * phi[1] + shell_terms * D_CA * ksi_shell[1]);
	psi_out[1] = phi[1] + tmp;
	eta_out[1] = phi[1] - tmp;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp = 0.5 * DT_SEC * (1 / tauCaSec[i - 1] * phi[k - 1] + shell_terms * D_CA * ksi_shell[k - 1]);
		psi_out[k - 1] = phi[k - 1] + tmp;
		eta_out[k - 1] = phi[k - 1] - tmp;
		tmp = 0.5 * DT_SEC * (1 / tauCaSec[i] * phi[k] + shell_terms * D_CA * ksi_shell[k]);
		psi_out[k] = phi[k] + tmp;
		eta_out[k] = phi[k] - tmp;
		tmp = 0.5 * DT_SEC * (1 / tauCaSec[i + 1] * phi[k + 1] + shell_terms * D_CA * ksi_shell[k + 1]);
		psi_out[k + 1] = phi[k + 1] + tmp;
		eta_out[k + 1] = phi[k + 1] - tmp;
		k += 3;
	}
	tmp = 0.5 * DT_SEC * (1 / tauCaSec[Params::N - 2] * phi[k - 1] + shell_terms * D_CA * ksi_shell[k - 1]);
	psi_out[k - 1] = phi[k - 1] + tmp;
	eta_out[k - 1] = phi[k - 1] - tmp;
	tmp = 0.5 * DT_SEC * (1 / tauCaSec[Params::N - 1] * phi[k] + shell_terms * D_CA * ksi_shell[k]);
	psi_out[k] = phi[k] + tmp;
	eta_out[k] = phi[k] - tmp;
}

void BuildCalciumFunctionTimeIntegration::integrate2(const TriMatN& mu_shell, const VectorN& l_muscle05, const VectorN& v_cell05, const TriMatN& nu_left, const TriMatN& nu_right, const VectorN& a_x_cell05_left, const VectorN& a_x_cell05_right,
													 const TriMatN& mu_cell, const TriMatN& psi, const TriMatN& eta, TriMatN& al_out, TriMatN& ar_out, int shell_terms) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	constexpr double D_CA = Params::CA_DIFFUSE_COEF;

	double tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[0] / l_muscle05[0] + (1.0 - shell_terms) * (1 / v_cell05[0]) * (- nu_right[0] * a_x_cell05_right[0] + 2 * mu_cell[0] / l_muscle05[0]));
	al_out[0] = psi[0] + tmp;
	ar_out[0] = eta[0] - tmp;
	tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[1] / l_muscle05[1] + (1.0 - shell_terms) * (1 / v_cell05[1]) * (- nu_right[1] * a_x_cell05_right[1] + 2 * mu_cell[1] / l_muscle05[1]));
	al_out[1] = psi[1] + tmp;
	ar_out[1] = eta[1] - tmp;
	int k = 3;
	for (int i = 1; i < Params::N - 1; i++) {
		tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[k - 1] / l_muscle05[i - 1] + (1.0 - shell_terms) * (1 / v_cell05[i - 1]) * (nu_left[k - 1] * a_x_cell05_left[i - 1] + 2 * mu_cell[k - 1] / l_muscle05[i - 1]));
		al_out[k - 1] = psi[k - 1] + tmp;
		ar_out[k - 1] = eta[k - 1] - tmp;
		tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[k] / l_muscle05[i] + (1.0 - shell_terms) * (1 / v_cell05[i]) * (nu_left[k] * a_x_cell05_left[i] - nu_right[k] * a_x_cell05_right[i] + 2 * mu_cell[k] / l_muscle05[i]));
		al_out[k] = psi[k] + tmp;
		ar_out[k] = eta[k] - tmp;
		tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[k + 1] / l_muscle05[i + 1] + (1.0 - shell_terms) * (1 / v_cell05[i + 1]) * (- nu_right[k + 1] * a_x_cell05_right[i + 1] + 2 * mu_cell[k + 1] / l_muscle05[i + 1]));
		al_out[k + 1] = psi[k + 1] + tmp;
		ar_out[k + 1] = eta[k + 1] - tmp;
		k += 3;
	}
	tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[k - 1] / l_muscle05[Params::N - 2] + (1.0 - shell_terms) * (1 / v_cell05[Params::N - 2]) * (nu_left[k - 1] * a_x_cell05_left[Params::N - 2] + 2 * mu_cell[k - 1] / l_muscle05[Params::N - 2]));
	al_out[k - 1] = psi[k - 1] + tmp;
	ar_out[k - 1] = eta[k - 1] - tmp;
	tmp = 0.5 * DT_SEC * D_CA * (shell_terms * 2 * mu_shell[k] / l_muscle05[Params::N - 1] + (1.0 - shell_terms) * (1 / v_cell05[Params::N - 1]) * (nu_left[k] * a_x_cell05_left[Params::N - 1] + 2 * mu_cell[k] / l_muscle05[Params::N - 1]));
	al_out[k] = psi[k] + tmp;
	ar_out[k] = eta[k] - tmp;
}

void BuildCalciumFunctionTimeIntegration::integrate3(const Params& params, const VectorN& tauCaSec, const TriMatN& ar, const TriMatN& chi_left, const TriMatN& chi_right, const VectorN1& Ica05, const TriMatN& mu_shell, const TriMatN& mu_cell,
													 const VectorN& CaShellCurr, const VectorN& CaShellPrev, const VectorN& CaCellCurr, const VectorN& CaCellPrev, const VectorN& l_muscle05, const VectorN& v_cell05, const VectorN& rho, VectorN& b_out, int shell_terms) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	constexpr double D_CA = Params::CA_DIFFUSE_COEF;
	constexpr double CA_REST = Params::CA_IN_REST;
	VectorN CaCellToL05{}, CaShellToLVcell05{}, Ica05_left{}, Ica05_right{}, Ica_left_vec{}, Ica_right_vec{}, mu_shell05_vec{}, mu_cell05_vec{};
	for (int i = 0; i < Params::N; i++) {
		CaCellToL05[i] = 0.5 * (3 * CaCellCurr[i] - CaCellPrev[i]) / l_muscle05[i];
		CaShellToLVcell05[i]  = 0.5 * (3 * CaShellCurr[i] - CaShellPrev[i]) / (l_muscle05[i] * v_cell05[i]);
		Ica05_left[i]  = Ica05[i];
		Ica05_right[i] = Ica05[i + 1];
	}
	if (shell_terms == 1){
		TriMatOperations::triMatVecMultN(ar, CaShellCurr, b_out);
	} else {
		TriMatOperations::triMatVecMultN(ar, CaCellCurr, b_out);
	}
	TriMatOperations::triMatVecMultN(chi_left,  Ica05_left,  Ica_left_vec);
	TriMatOperations::triMatVecMultN(chi_right, Ica05_right, Ica_right_vec);
	TriMatOperations::triMatVecMultN(mu_shell, CaCellToL05,  mu_shell05_vec);
	TriMatOperations::triMatVecMultN(mu_cell, CaShellToLVcell05, mu_cell05_vec);

	constexpr double K = 1.0 / (Params::Z * Params::F);
	for (int i = 0; i < Params::N; i++) {
		b_out[i] += DT_SEC * (CA_REST / tauCaSec[i] * rho[i]
							  +		 shell_terms  * (D_CA * 2 * mu_shell05_vec[i] - K * (Ica_left_vec[i] + Ica_right_vec[i]))
							  + (1 - shell_terms) * (D_CA * 2 * mu_cell05_vec[i]));
	}
}
