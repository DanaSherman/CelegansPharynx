#ifndef BUILD_AXIAL_CURRENT_FUNCTION_H_
#define BUILD_AXIAL_CURRENT_FUNCTION_H_

#include "base.h"
#include "params.h"
#include "tri_mat_operations.h"

class BuildAxialCurrentFunction {
 public:
	static void buildXi(const VectorN1& x, TriMatN1& out_xi);
	static void buildChi(const VectorN1& a, TriMatN1& out_chi);
};

#endif /* BUILD_AXIAL_CURRENT_FUNCTION_H_ */
