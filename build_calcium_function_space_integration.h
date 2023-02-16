#ifndef BUILD_CALCIUM_FUNCTION_SPACE_INTEGRATION_H_
#define BUILD_CALCIUM_FUNCTION_SPACE_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildCalciumFunctionSpaceIntegration {
 public:
	static void LHS(const VectorN& segX, TriMatN& phi_out);
	static void RHS1(const VectorN& segX, const VectorN& a, const VectorN& v, TriMatN& chi_out);
	static void RHS2(const VectorN& segX, VectorN& rho_out);
	static void RHS3_left (TriMatN& nu_left_out);
	static void RHS3_right(TriMatN& nu_right_out);
	static void RHS4_left_minus_right(const VectorN& a_l, const VectorN& a_r, const VectorN& v, TriMatN& ksi_out);
	static void RHS5(const VectorN& segX, const VectorN& a, const VectorN& v, TriMatN& mu_shell_out);
	static void RHS6(const VectorN& segX, const VectorN& a, TriMatN& mu_cell_out);
};

#endif	// BUILD_CALCIUM_FUNCTION_SPACE_INTEGRATION_H_
