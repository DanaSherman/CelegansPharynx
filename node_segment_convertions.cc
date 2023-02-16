#include "node_segment_convertions.h"

void NodeSegmentConvertions::calcLMuscle05Seg(const VectorN& l_muscle_next, const VectorN& l_muscle, VectorN& l_muscle05_out) {
	for (int i = 0; i < Params::N; i++) {
		l_muscle05_out[i] = 0.5 * (3 * l_muscle_next[i] - l_muscle[i]);					//l_muscle[i+0.5, timeStep+1.5]
	}
}

void NodeSegmentConvertions::calcZ(const VectorN& segR, const VectorN& l_muscle, VectorN& seg_z_out) {
	for (int i = 0; i < Params::N; i++) {
		seg_z_out[i] = segR[i] - l_muscle[i];
		if (seg_z_out[i] < 0) {seg_z_out[i] = 0;}
	}
}

void NodeSegmentConvertions::calcZNodes(const VectorN1& x, const VectorN& segX, const VectorN& z, VectorN1& z_nodes_out) {
	for (int i = 1; i < Params::N; i++) {												//z_node[i]
		z_nodes_out[i] = 0.5 * (z[i - 1] + z[i]);
	}
	z_nodes_out[0] = 2 * z[0] - z_nodes_out[1];											//z_node[0]
	if (z_nodes_out[0] < 0) {z_nodes_out[0] = 0;}
	z_nodes_out[Params::N] = 2 * z[Params::N - 1] - z_nodes_out[Params::N - 1];			//z05_node[N]
	if (z_nodes_out[Params::N] < 0) {z_nodes_out[Params::N] = 0;}
}

void NodeSegmentConvertions::setZAll(const VectorN1& z_nodes, const VectorN& z, Vector2N1& z_all_out) {
	for (int i = 0; i < Params::N; i++) {
		z_all_out[2 * i]	 = z_nodes[i];
		z_all_out[2 * i + 1] = z[i];
	}
	z_all_out[2 * Params::N] = z_nodes[Params::N];
}

void NodeSegmentConvertions::calcZNodesMax(const VectorN1& x, const VectorN& segX, const VectorN& z_max, VectorN1& z_nodes_max_out) {
	double tmp;
	for (int i = 1; i < Params::N; i++) {
		z_nodes_max_out[i] = 0.5 * (z_max[i] + z_max[i + 1]);								//z_node_max[i, timeStep+0.5]
	}
	z_nodes_max_out[0] = 2 * z_max[0] - z_nodes_max_out[1];									//z_node_max[0, timeStep+0.5]
	z_nodes_max_out[Params::N] = 2 * z_max[Params::N - 1] - z_nodes_max_out[Params::N - 1];	//z_node_max[N, timeStep+0.5]
}

void NodeSegmentConvertions::calcAxCell05Left(const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& a_x_cell05_left_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;																	//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		a_x_cell05_left_out[i] = M_PI * square(CELL_R * r[i]) - 3 * sqrt(3) * z_nodes_max[i] * z05_nodes[i];				//a_x_cell05_left_out[i, timeStep+1.5]
	}
}

void NodeSegmentConvertions::calcAxCell05Right(const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& a_x_cell05_right_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;																	//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		a_x_cell05_right_out[i] = M_PI * square(CELL_R * r[i + 1]) - 3 * sqrt(3) * z_nodes_max[i + 1] * z05_nodes[i + 1];	//a_x_cell05_right_out[i, timeStep+1.5]
	}
}

void NodeSegmentConvertions::calcCaShellNodes(const VectorN& CaShellNext, const VectorN& VShellLeft, const VectorN& VShellRight, VectorN1& ca_shell_nodes_out) {
	double caShellLeft;							//Ca[i-0.5, timeStep+1]
	double CaShellRight = CaShellNext[0];		//Ca[i+0.5, timeStep+1]
	ca_shell_nodes_out[0] = CaShellRight;
	for (int i = 1; i < Params::N; i++) {
		caShellLeft = CaShellRight;
		CaShellRight = CaShellNext[i];
		ca_shell_nodes_out[i] = (caShellLeft * VShellLeft[i] + CaShellRight * VShellRight[i]) / (VShellLeft[i] + VShellRight[i]);	//Ca[i, timeStep+1]
																																	//Notes: (1): VShell(Lef,Righ)t are VectorN, with VShellLeft[0] and VShellRight[0] being irrelevant for calculations
																																	//		 (2): VShell(Lef,Righ)t are in [m^3]
	}
	ca_shell_nodes_out[Params::N] = CaShellRight;
}

void NodeSegmentConvertions::calcCaShell05Nodes(const VectorN& CaShellNext, const VectorN& CaShellCurr, const VectorN& VShellLeft, const VectorN& VShellRight, VectorN1& ca_shell05_nodes_out) {
	double caShellLeft05;														//Ca[i-0.5, timeStep+1.5]
	double CaShellRight05 = 0.5 * (3 * CaShellNext[0] - CaShellCurr[0]);		//Ca[i+0.5, timeStep+1.5]
	ca_shell05_nodes_out[0] = CaShellRight05;
	for (int i = 1; i < Params::N; i++) {
		caShellLeft05 = CaShellRight05;
		CaShellRight05 = 0.5 * (3 * CaShellNext[i] - CaShellCurr[i]);
		ca_shell05_nodes_out[i] = (caShellLeft05 * VShellLeft[i] + CaShellRight05 * VShellRight[i]) / (VShellLeft[i] + VShellRight[i]);		//Ca[i, timeStep+1.5]
																																			//Notes: (1) VShell(Lef,Righ)t05 are VectorN, with VShellLeft[0] and VShellRight[0] being irrelevant for calculations
																																			//		 (2) VShell(Lef,Righ)t are in [m^3]
	}
	ca_shell05_nodes_out[Params::N] = CaShellRight05;
}

void NodeSegmentConvertions::calcVCell05(const VectorN& segL, const VectorN1& r, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, VectorN& v_cell05_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;								//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		v_cell05_out[i] = 1 / 3.0 * M_PI * segL[i] * square(CELL_R) * (square(r[i]) + square(r[i + 1]) + r[i] * r[i + 1])
						  - sqrt(3) * segL[i] * (z05_nodes[i] * z_nodes_max[i] + z05_nodes[i + 1] * z_nodes_max[i + 1] + sqrt(z05_nodes[i] * z_nodes_max[i] * z05_nodes[i + 1] * z_nodes_max[i + 1]));	//v_cell05_out[i, timeStep+1.5]
	}
}
