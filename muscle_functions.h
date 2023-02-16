#ifndef MUSCLE_FUNCTIONS_H_
#define MUSCLE_FUNCTIONS_H_

#include "base.h"
#include "params.h"

class MuscleFunctions {
 public:
	static void calc_m_PED_CEK_PEK (const VectorN& v_seg_init, VectorN& m, VectorN& PE_D, VectorN& CE_K, VectorN& PE_K, boolVecN1* print_in_out);
	static double calcK1 (int i, double ca_in, boolVecN1* print);
	static double calcMuscleContractionFactor(double t, double AP_t0, int i, boolVecN1* print, double x);
	static double calcMLCP_P (int i, double ca_in, double CF, double MLCP_P_prev, boolVecN1* print);
	static double calcMuscleActivation (double L0, double l_muscle);
	static void   buildDistributionsFromNdistMoments (int seg, const MomentsVec& Q_prev, Distribution& n_in_out, Distribution& n0l_in_out);
	static double calcBeta (int moment, double M, double Mp);
	static double calcMusclePhi1 (int moment, double Q0_prev, Distribution n0l, double M, double Mp);
	static void	  calcMusclePhi2 (int moment, double Q0_prev, Distribution n, double AM, double AMp, double* phi2s);
	static double calcF_OE (int i, const VectorN& l_muscle, const VectorN& segR);
	static double G (int moment, double z, Distribution dist);
	static double F (int m, double z);
	static int    factorial (int n);
};

#endif /* MUSCLE_FUNCTIONS_H_ */
