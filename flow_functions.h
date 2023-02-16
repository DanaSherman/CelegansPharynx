#ifndef FLOW_FUNCTIONS_H_
#define FLOW_FUNCTIONS_H_

#include "base.h"
#include "params.h"

class FlowFunctions {
 public:
	static void loadZMax (VectorN& z_max);
	static void calcLumenArea (const VectorN1& z_nodes_max, const VectorN1& z05_nodes, const VectorN& z05_max, const VectorN& z05, Vector2N1& area_out, boolVecN1* print_in_out);
	static void calcLumenVolume (const VectorN1& x, const VectorN& segL, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, const VectorN& z05_max, const VectorN& z05, Vector2N& volume_out, boolVecN1* print_in_out);
	static void calcFlow (const VectorN& segL, const Vector2N& volume, const Vector2N& volume_prev, Vector2N& flow_out);
	static void calcVelocity (const VectorN1& z_nodes_max, const VectorN& z_max, const Vector2N1& area, const Vector2N& flow, int mean, Vector2N1& velocity_out);
	static void calcParticlePosition (int timeStep, const Vector2N1& x_all, const Vector2N1& z05_all, double v_coef, const Vector2N1& velocity, int mean, Particle* in_out_part);
};

#endif /* FLOW_FUNCTIONS_H_ */
