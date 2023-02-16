#ifndef BUILD_CABLE_FUNCTION_TIME_INTEGRATION_H_
#define BUILD_CABLE_FUNCTION_TIME_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildCableFunctionTimeIntegration {
 public:
	static void integrate1(const TriMatN1& phi, TriMatN1& eta_in_out, TriMatN1& psi_out);
	static void integrate2(const TriMatN1& phi, const TriMatN1& psi, const TriMatN1& eta, const VectorN1& g, TriMatN1& al_out, TriMatN1& ar_out);
	static void integrate3(const TriMatN1& ar, const VectorN1& v, const TriMatN1& phi, const VectorN1& d, const VectorN1& x, double t, double t0, double& a_now, VectorN1& b_out);
};

#endif /* BUILD_CABLE_FUNCTION_TIME_INTEGRATION_H_ */
