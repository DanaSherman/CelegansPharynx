#ifndef BUILD_LUMEN_FUNCTION_TIME_INTEGRATION_H_
#define BUILD_LUMEN_FUNCTION_TIME_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildLumenFunctionTimeIntegration {
 public:
	static void integrate1 (const TriMatN& ksi, TriMatN& eta_in_out, TriMatN& psi_out);
	static void integrate2 (const TriMatN& psi, const TriMatN& eta, TriMatN& al_out, TriMatN& ar_out);
	static void integrate3 (const TriMatN& ar, const VectorN& l_muscle, const TriMatN& chi, const TriMatN& phi, const VectorN& F_CE, const VectorN& F_OE, const VectorN& rho, const VectorN& v_muscle_curr, const VectorN& v_muscle_next, VectorN& b_out);
};

#endif /* BUILD_LUMEN_FUNCTION_TIME_INTEGRATION_H_ */
