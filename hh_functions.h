
#ifndef HH_FUNCTIONS_H_
#define HH_FUNCTIONS_H_

#include "base.h"
#include "params.h"

class HHFunctions {
 public:
	static double calcGateOpeningProb  (const Params& params, GatePart0 gate_part0, GatePart1 gate_part1, double v, double ca_shell05_nodes, double HH);
	static double calcGateOpeningProbHP(const Params& params, GatePart0 gate_part0, GatePart1 gate_part1, double v, double HH);
	static AlphaBeta alphaBeta_m_cca1  (const Params& params, double V);
	static AlphaBeta alphaBeta_h_cca1  (const Params& params, double V);
	static AlphaBeta alphaBeta_m_egl19 (const Params& params, double V);
	static AlphaBeta alphaBeta_f_egl19 (const Params& params, double V, double ca_shell05_nodes);
	static AlphaBeta alphaBeta_m_exp2  (const Params& params, double V);
	static AlphaBeta alphaBeta_h_exp2  (const Params& params, double V);
	static AlphaBeta alphaBeta_m_x     (const Params& params, double V);
	static AlphaBeta alphaBeta_h_x     (const Params& params, double V);
	static double calcC05(const Params& params, double VPrev, double VCurr, double hm205_cca1, double fm205_egl19);
	static double calcE05(double VPrev, double VCurr);
	static double calcGca(double VPrev, double VCurr, double CaIn05);
 private:
	static AlphaBeta alphaBeta_m_cca1_HP (const Params& params, double V);
	static AlphaBeta alphaBeta_m_exp2_HP (const Params& params, double V);
	static AlphaBeta alphaBeta_h_exp2_HP (const Params& params, double V);
};

#endif	// HH_FUNCTIONS_H_
