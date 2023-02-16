#include "flow_functions.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <sstream>
using namespace std;

void FlowFunctions::loadZMax (VectorN& z_max) {
	std::ifstream maximal_inner_r_file("Maximal_inner_radius_um.csv");
	std::cout << "\nDon't forget to update Maximal_inner_radius_um.csv for flow calculations\n" << std::endl;
	if(!maximal_inner_r_file.is_open()) throw std::runtime_error("Could not open file\n");
	std::string line;
	if (std::getline(maximal_inner_r_file, line)) {
		std::stringstream ss(line);
		int col = 0;
		double val;
		while (ss >> val && col < Params::N){													// Extract each double in the current line
			z_max[col] = val * 1E-4;															//Convert from [um] to [cm]
			// If the next token is a comma, ignore it and move on
			if(ss.peek() == ',') ss.ignore();
			col++;
		}
	}
	maximal_inner_r_file.close();
}

void FlowFunctions::calcLumenArea(const VectorN1& z_nodes_max, const VectorN1& z05_nodes, const VectorN& z_max, const VectorN& z05, Vector2N1& a_out, boolVecN1* print_in_out) {
	for (int i = 0; i < Params::N; i++) {
		a_out[2 * i]     = 3 * std::sqrt(3) * z05_nodes[i] * z_nodes_max[i];
		a_out[2 * i + 1] = 3 * std::sqrt(3) * z05[i] 	   * z_max[i];
	}
	a_out[2 * Params::N] = 3 * std::sqrt(3) * z05_nodes[Params::N] * z_nodes_max[Params::N];
}

void FlowFunctions::calcLumenVolume(const VectorN1& x, const VectorN& segL, const VectorN1& z_nodes_max, const VectorN1& z05_nodes, const VectorN& z_max, const VectorN& z05, Vector2N& v_out, boolVecN1* print_in_out) {
	for (int i = 0; i < Params::N; i++) {
		v_out[2 * i]	 = sqrt(3) * 0.5 * segL[i] * (z05_nodes[i] * z_nodes_max[i] + z05[i]		   * z_max[i]			+ sqrt(z05_nodes[i] * z_nodes_max[i] * z05[i]			* z_max[i]));
		v_out[2 * i + 1] = sqrt(3) * 0.5 * segL[i] * (z05[i]	   * z_max[i]		+ z05_nodes[i + 1] * z_nodes_max[i + 1] + sqrt(z05[i]		* z_max[i]		 * z05_nodes[i + 1] * z_nodes_max[i + 1]));
	}
}

void FlowFunctions::calcFlow(const VectorN& segL, const Vector2N& V, const Vector2N& V_prev, Vector2N& flow_out) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	int most_posterior_open_segment = Params::POST_ISTH_START_NODE;
	for (int i = 2 * Params::N - 1; i >= 2 * most_posterior_open_segment; i--) {
		flow_out[i] = 0.0;
	}
	flow_out[2 * most_posterior_open_segment - 1] = (V[2 * most_posterior_open_segment - 1] - V_prev[2 * most_posterior_open_segment - 1]) / DT_SEC;
	for (int i = 2 * most_posterior_open_segment - 2; i >= 0; i--) {
		flow_out[i] = (V[i] - V_prev[i]) / DT_SEC + flow_out[i + 1];
	}
}

void FlowFunctions::calcVelocity(const VectorN1& z_nodes_max, const VectorN& z_max, const Vector2N1& A, const Vector2N& flow, int mean, Vector2N1& velocity_out) {
	double A_max_node, A_max_seg, O_node, O_seg;
	double v_coef_node, v_coef_seg;						//Coefficient of the mean velocity when particle is positioned at the center of the lumen [unit-less]
	if (mean == 1) {
		v_coef_node = 1.0;
		v_coef_seg  = 1.0;
		for (int i = 0; i < Params::POST_ISTH_START_NODE; i++) {
			if (A[2 * i] > 0)	  {velocity_out[2 * i]	   = v_coef_node * flow[2 * i] / A[2 * i];}
			else {velocity_out[2 * i] = 0.0;}
			if (A[2 * i + 1] > 0) {velocity_out[2 * i + 1] = v_coef_seg  * flow[2 * i + 1] / A[2 * i + 1];}
			else {velocity_out[2 * i + 1] = 0.0;}
		}
	} else {
		for (int i = 0; i < Params::POST_ISTH_START_NODE; i++) {
			A_max_node = 3 * std::sqrt(3) * square(z_nodes_max[i]);
			A_max_seg  = 3 * std::sqrt(3) * square(z_max[i]);
			O_node = A[2 * i]	  / A_max_node;
			O_seg  = A[2 * i + 1] / A_max_seg;
			v_coef_node = 1.0097 * square(O_node) - 2.366 * O_node + 3.5587;			//Eq. taken from Avery and Shtonda, 2003
			v_coef_seg  = 1.0097 * square(O_seg)  - 2.366 * O_seg  + 3.5587;
			if (A[2 * i] > 0)	  {velocity_out[2 * i]	   = v_coef_node * flow[2 * i] / A[2 * i];}
			else {velocity_out[2 * i] = 0.0;}
			if (A[2 * i + 1] > 0) {velocity_out[2 * i + 1] = v_coef_seg  * flow[2 * i + 1] / A[2 * i + 1];}
			else {velocity_out[2 * i + 1] = 0.0;}
		}
	}
	for (int i = 2 * Params::POST_ISTH_START_NODE; i < 2 * Params::N + 1; i++) {
		velocity_out[i] = 0.0;
	}
}

void FlowFunctions::calcParticlePosition(int timeStep, const Vector2N1& x_all, const Vector2N1& z_all, double v_coef, const Vector2N1& velocity, int mean, Particle* particle) {
	constexpr double DT_SEC = 1E-3 * Params::DT;
	double most_post_x = x_all[2 * Params::POST_ISTH_START_NODE];
	Particle& part = *particle;
	double curr_idx = part.seg, curr_pos = part.position, delta_t, delta_x;
	Vector2N1 vel;
	double v_mean_coef = 1.0;
	if (mean == 1) {v_mean_coef = v_coef;}
	for (int i = 0; i < 2 * Params::N + 1; i++) {
		vel[i] = v_mean_coef * velocity[i];
	}
	double left_boundary_r, right_boundary_r, curr_r;
	left_boundary_r = z_all[curr_idx];
	right_boundary_r = z_all[curr_idx + 1];
	curr_r = left_boundary_r  + (right_boundary_r - left_boundary_r) * (curr_pos - x_all[curr_idx]) / (x_all[curr_idx + 1] - x_all[curr_idx]);

	bool stuck = false, reach_end = false, time_over = false;
	if (2 * curr_r <= part.diameter) {stuck = true;}
	if ((curr_pos == most_post_x && vel[curr_idx] > 0) || (curr_pos == x_all[0] && vel[curr_idx] < 0)) {reach_end = true;}

	double time_left = DT_SEC;
	while (!stuck && !reach_end && !time_over) {

		// Try to reach left/right end of current segment
		if (vel[curr_idx] == 0) {
			time_over = true;
			delta_x = 0;
			delta_t = 0;
		} else {
			if (vel[curr_idx] > 0) {
				delta_x = x_all[curr_idx + 1] - curr_pos;
			} else {
				delta_x = curr_pos - x_all[curr_idx];
			}
		}
		delta_t = delta_x / std::abs(vel[curr_idx]);

		// Check if have enough time to reach left/right end of current segment
		if (delta_t > time_left) {
			time_over = true;
			delta_t = time_left;
			delta_x = time_left * std::abs(vel[curr_idx]);
		}

		// Check if get stuck before reaching left/right end of current segment
		if (vel[curr_idx] == 0) {
			// no need to update curr_pos, curr_idx
		} else {
			if (vel[curr_idx] > 0) {
				if (2 * right_boundary_r <= part.diameter) {
					stuck = true;
					curr_pos = x_all[curr_idx] + (0.5 * part.diameter - left_boundary_r) * (x_all[curr_idx + 1] - x_all[curr_idx]) / (right_boundary_r - left_boundary_r);
					// curr_idx doesn't change
				} else {
					curr_pos += delta_x;
					if (!time_over) {curr_idx++;}																	//otherwise stay at the current node
					time_left -= delta_t;
					if (curr_pos >= most_post_x) {
						reach_end = true;
						curr_pos = most_post_x;
						curr_idx = 2 * Params::POST_ISTH_START_NODE - 1;
					}
				}
			} else {
				if (2 * left_boundary_r <= part.diameter) {
					stuck = true;
					curr_pos = x_all[curr_idx] + (0.5 * part.diameter - left_boundary_r) * (x_all[curr_idx + 1] - x_all[curr_idx]) / (right_boundary_r - left_boundary_r);
					// curr_idx doesn't change
				} else {
					curr_pos -= delta_x;
					if (!time_over) {curr_idx--;}																	//otherwise stay at the current node
					time_left -= delta_t;

					if (curr_pos <= x_all[0]) {
						reach_end = true;
						curr_pos = x_all[0];
						curr_idx = 0;
					}
				}
			}
		}

		// Update left- and right boundary positions
		left_boundary_r = z_all[curr_idx];
		right_boundary_r = z_all[curr_idx + 1];
		curr_r = left_boundary_r  + (right_boundary_r - left_boundary_r) * (curr_pos - x_all[curr_idx]) / (x_all[curr_idx + 1] - x_all[curr_idx]);
	}

	// Update particle data
	part.position = curr_pos;
	part.seg = curr_idx;
}
