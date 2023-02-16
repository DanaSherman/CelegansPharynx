#ifndef BUILD_LUMEN_FUNCTION_SPACE_INTEGRATION_H_
#define BUILD_LUMEN_FUNCTION_SPACE_INTEGRATION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildLumenFunctionSpaceIntegration {
 public:
	static void LHS(const VectorN& segX, const VectorN& C, TriMatN& mat_out);
	static void RHS1(const VectorN& segX, TriMatN& phi_out);
	static void RHS2(const VectorN& segX, const VectorN& K, const VectorN& R, VectorN& vec_out);
};

#endif /* BUILD_LUMEN_FUNCTION_SPACE_INTEGRATION_H_ */
