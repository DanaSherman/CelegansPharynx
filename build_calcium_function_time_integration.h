#ifndef BUILD_CALCIUM_FUNCTION_TIME_INTEGRATION_H_
#define BUILD_CALCIUM_FUNCTION_TIME_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildCalciumFunctionTimeIntegration {
 public:
	static void integrate1(const Params& params, const VectorN& tauCaSec, const TriMatN& phi, const TriMatN& ksi_shell, TriMatN& eta_out, TriMatN& psi_out, int shell_terms);
	static void integrate2(const TriMatN& mu_shell, const VectorN& l_muscle05, const VectorN& v_cell05, const TriMatN& nu_left, const TriMatN& nu_right, const VectorN& a_x_cell05_left, const VectorN& a_x_cell05_right, const TriMatN& mu_cell, const TriMatN& psi, const TriMatN& eta, TriMatN& al_out, TriMatN& ar_out, int shell_terms);
	static void integrate3(const Params& params, const VectorN& tauCaSec, const TriMatN& ar, const TriMatN& chi_left, const TriMatN& chi_right, const VectorN1& Ica05, const TriMatN& muCa_shell, const TriMatN& muCa_cell,
						   const VectorN& ca_shell, const VectorN& ca_shell_prev, const VectorN& ca_cell, const VectorN& ca_cell_prev, const VectorN& l_muscle05, const VectorN& v_cell05, const VectorN& rhoCa, VectorN& b_out, int shell_terms);
};

#endif	// BUILD_CALCIUM_FUNCTION_TIME_INTEGRATION_H_
