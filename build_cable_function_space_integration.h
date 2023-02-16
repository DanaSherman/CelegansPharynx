#ifndef BUILD_CABLE_FUNCTION_SPACE_INTEGRATION_H_
#define BUILD_CABLE_FUNCTION_SPACE_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildCableFunctionSpaceIntegration {
 public:
	static void LHS(const VectorN1& x, const VectorN1& p, TriMatN1& out_phi);
	static void RHS(const VectorN1& x, const VectorN1& a, TriMatN1& out_eta);
};

#endif	// BUILD_CABLE_FUNCTION_SPACE_INTEGRATION_H_
