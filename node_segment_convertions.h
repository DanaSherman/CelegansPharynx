#ifndef NODE_SEGMENT_CONVERTIONS_H_
#define NODE_SEGMENT_CONVERTIONS_H_

#include "base.h"
#include "params.h"

class NodeSegmentConvertions {
 public:
	static void calcLMuscle05Seg(const VectorN& l_muscle_next, const VectorN& l_muscle, VectorN& l_muscle05_out);
	static void calcZ(const VectorN& segR, const VectorN& l_muscle, VectorN& z_out);
	static void calcZNodes(const VectorN1& x, const VectorN& segX, const VectorN& z, VectorN1& z_nodes_out);
	static void setZAll(const VectorN1& z_nodes, const VectorN& z, Vector2N1& z_all_out);
	static void calcZNodesMax(const VectorN1& x, const VectorN& segX, const VectorN& z_max, VectorN1& z_nodes_max_out);
	static void calcAxCell05Left (const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& a_x_cell_left_out);
	static void calcAxCell05Right(const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& a_x_cell_right_out);
	static void calcCaShellNodes(const VectorN& ca_shell, const VectorN& v_shell_left,  const VectorN& v_shell_right,  VectorN1& ca_shell_nodes_out);
	static void calcCaShell05Nodes(const VectorN& ca_shell_next, const VectorN& ca_shell, const VectorN& v_shell_left,  const VectorN& v_shell_right,  VectorN1& ca_shell05_nodes_out);
	static void calcVCell05(const VectorN& segL, const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& v_cell05_out);
};

#endif	// NODE_SEGMENT_CONVERTIONS_H_
