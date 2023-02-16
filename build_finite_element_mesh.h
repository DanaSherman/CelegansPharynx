#ifndef BUILD_FINITE_ELEMENT_MESH_H_
#define BUILD_FINITE_ELEMENT_MESH_H_

#include "base.h"
#include "params.h"

class BuildFiniteElementMesh {
 public:
	static void buildMesh(VectorN1& x_out, VectorN1& r_out);
	static void setPerimeter(const VectorN1& r, VectorN1& p_out);
	static void setCrossSection(const VectorN1& r, VectorN1& a_out);
};

#endif	// BUILD_FINITE_ELEMENT_MESH_H_
