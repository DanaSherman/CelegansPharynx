#ifndef INIT_SEGMENTS_H_
#define INIT_SEGMENTS_H_

#include "base.h"
#include "params.h"

class InitSegments {
 public:
	static void setSegX(const VectorN1& x, VectorN& segX_out);
	static void setXAll(const VectorN1& x, const VectorN& segX, Vector2N1& x_all_out);
	static void setSegLength(const VectorN1& x, VectorN& segL_out);
	static void setSegRadius(const VectorN1& x, const VectorN& segX, const VectorN1& r, VectorN& segR_out);
	static void setAIca(const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& a_Ica_left_out, VectorN& a_Ica_right_out);
	static void setAxShellLeft(const VectorN1& r, VectorN& a_x_shell_left_out);
	static void setAxShellRight(const VectorN1& r, VectorN& a_x_shell_right_out);
	static void setAz(const VectorN& segL, const VectorN1& r, VectorN& a_z_out);
	static void setVShell(const VectorN& segL, const VectorN1& r, VectorN& v_shell_out);
	static void setVShellLeft (const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& v_shell_left_out);
	static void setVShellRight(const VectorN& segL, const VectorN1& r, const VectorN& segR, VectorN& v_shell_right_out);
	static void setVSegmentInit(const VectorN& segL, const VectorN1& r, VectorN& v_seg_init_out);
};

#endif	// INIT_SEGMENTS_H_
