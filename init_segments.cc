#include "init_segments.h"
#include <cmath>

/* Calculate the position at the middle of each segment along the pharynx (x-axis) */
void InitSegments::setSegX(const VectorN1& x, VectorN& segX_out) {
	for (int i = 0; i < Params::N; i++) {
		segX_out[i] = 0.5 * (x[i] + x[i + 1]);
	}
}

void InitSegments::setXAll(const VectorN1& x, const VectorN& segX, Vector2N1& x_all_out) {
	for (int i = 0; i < Params::N; i++) {
		x_all_out[2 * i]	 = x[i];
		x_all_out[2 * i + 1] = segX[i];
	}
	x_all_out[2 * Params::N] = x[Params::N];
}

void InitSegments::setSegLength(const VectorN1& x, VectorN& segL_out) {
	for (int i = 0; i < Params::N; i++) {
		segL_out[i] = x[i + 1] - x[i];
	}
}

/* Calculate the radius of each segment at rest along the pharynx (z-axis) */
void InitSegments::setSegRadius(const VectorN1& x, const VectorN& segX, const VectorN1& r, VectorN& segR_out) {
	for (int i = 0; i < Params::N; i++) {
		segR_out[i] = 0.5 * (r[i] + r[i + 1]);
	}
}

void InitSegments::setAIca(const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& a_Ica_left_out, VectorN& a_Ica_right_out) {
	double l_left, l_right;
	for (int i = 0; i < Params::N; i++) {
		l_left  = sqrt(square(0.5 * segL[i]) + square(r[i]    - segR[i]));
		l_right = sqrt(square(0.5 * segL[i]) + square(segR[i] - r[i + 1]));
		a_Ica_left_out[i]  = M_PI * l_left  * (r[i]    + segR[i]);
		a_Ica_right_out[i] = M_PI * l_right * (segR[i] + r[i + 1]);
	}
}

void InitSegments::setAxShellLeft(const VectorN1& r, VectorN& a_x_shell_left_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		a_x_shell_left_out[i]  = M_PI * (1 - square(CELL_R)) * square(r[i]);
	}
}

void InitSegments::setAxShellRight(const VectorN1& r, VectorN& a_x_shell_right_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		a_x_shell_right_out[i]  = M_PI * (1 - square(CELL_R)) * square(r[i + 1]);
	}
}

void InitSegments::setAz(const VectorN& segL, const VectorN1& r, VectorN& a_z_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	double l;
	for (int i = 0; i < Params::N; i++) {
		l = sqrt(square(segL[i]) + square(CELL_R * (r[i] - r[i + 1])));
		a_z_out[i]  = M_PI * l * CELL_R * (r[i] + r[i + 1]);
	}
}

void InitSegments::setVShell(const VectorN& segL, const VectorN1& r, VectorN& v_shell_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		v_shell_out[i]  = 1 / 3.0 * M_PI * segL[i] * (1 - square(CELL_R)) * (square(r[i]) + square(r[i + 1]) + r[i] * r[i + 1]);
	}
}

void InitSegments::setVShellLeft(const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& v_shell_left_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	v_shell_left_out[0] = 0;
	v_shell_left_out[0] *= 1E-6;								//Convert from [cm^3] to [m^3]
	for (int i = 1; i < Params::N; i++) {
		v_shell_left_out[i] = 1 / 3.0 * M_PI * 0.5 * segL[i - 1] * (1 - square(CELL_R)) * (square(segR[i - 1]) + square(r[i]) + segR[i - 1] * r[i]);
		v_shell_left_out[i] *= 1E-6;							//Convert from [cm^3] to [m^3]
	}
}

void InitSegments::setVShellRight(const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& v_shell_right_out) {
	constexpr double CELL_R = 1 - Params::SHELL_THICKNESS;		//in % of segment's radius
	for (int i = 0; i < Params::N; i++) {
		v_shell_right_out[i] = 1 / 3.0 * M_PI * 0.5 * segL[i] * (1 - square(CELL_R)) * (square(r[i]) + square(segR[i]) + r[i] * segR[i]);
		v_shell_right_out[i] *= 1E-6;							//Convert from [cm^3] to [m^3]
	}
}

void InitSegments::setVSegmentInit(const VectorN& segL, const VectorN1& r, VectorN& v_seg_init_out) {
		for (int i = 0; i < Params::N; i++) {
			v_seg_init_out[i]  = 1 / 3.0 * M_PI * segL[i] * (square(r[i]) + square(r[i + 1]) + r[i] * r[i + 1]);
		}
}
