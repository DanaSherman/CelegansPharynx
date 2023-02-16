#include "hh_functions.h"

#include <cmath>

double HHFunctions::calcGateOpeningProb(const Params& params, GatePart0 gate_part0, GatePart1 gate_part1, double v, double ca_shell05_nodes, double HH) {
	double aval;
	double bval;
	AlphaBeta abvals{};
	if (gate_part1 == GatePart1::cca1) {
		if (gate_part0 == GatePart0::m) {
			abvals = alphaBeta_m_cca1(params, v);
		} else {
			abvals = alphaBeta_h_cca1(params, v);
		}
	} else if (gate_part1 == GatePart1::egl19) {
		if (gate_part0 == GatePart0::m) {
			abvals = alphaBeta_m_egl19(params, v);
		} else {
			abvals = alphaBeta_f_egl19(params, v, ca_shell05_nodes);
		}
	} else {	//(gate_part1 == GatePart1::exp2)
		if (gate_part0 == GatePart0::m) {
			abvals = alphaBeta_m_exp2(params, v);
		} else {
			abvals = alphaBeta_h_exp2(params, v);
		}
	}
	aval = Params::DT * abvals.a_prob;
	bval = Params::DT * abvals.b_prob;

	double tmp = 0.5 * (aval + bval);
	return ((aval + (1.0 - tmp) * HH) / (1.0 + tmp));
}

double HHFunctions::calcGateOpeningProbHP(const Params& params, GatePart0 gate_part0, GatePart1 gate_part1, double v, double HH) {
	double aval;
	double bval;
	AlphaBeta abvals{};
	if (gate_part1 == GatePart1::cca1) {
		abvals = alphaBeta_m_cca1_HP(params, v);
	} else {
		if (gate_part0 == GatePart0::m) {
			abvals = alphaBeta_m_exp2_HP(params, v);
		} else {
			abvals = alphaBeta_h_exp2_HP(params, v);
		}
	}
	aval = Params::DT * abvals.a_prob;
	bval = Params::DT * abvals.b_prob;

	double tmp = 0.5 * (aval + bval);
	return ((aval + (1.0 - tmp) * HH) / (1.0 + tmp));
}

AlphaBeta HHFunctions::alphaBeta_m_cca1(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = params.get(GeneType::cca1_mInf_A2);
	double v0 = params.get(GeneType::cca1_mInf_v0);
	double dv = params.get(GeneType::cca1_mInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double m_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::cca1_mTau_act_A1);
	A2 = params.get(GeneType::cca1_mTau_act_A2);
	v0 = params.get(GeneType::cca1_mTau_act_v0);
	dv = params.get(GeneType::cca1_mTau_act_dv);
	tmp = std::exp(-(v - v0) / dv);
	double tau_m = A1 + (A2 - A1) / (1 + tmp);

	double aProb = m_inf / tau_m;
	double bProb = (1 - m_inf) / tau_m;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_h_cca1(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::cca1_hInf_v0);
	double dv = params.get(GeneType::cca1_hInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double h_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::cca1_hTau_A1);
	A2 = params.get(GeneType::cca1_hTau_A2);
	v0 = params.get(GeneType::cca1_hTau_v0);
	dv = params.get(GeneType::cca1_hTau_dv);
	tmp = std::exp(-(v - v0) / dv);
	double tau_h = A1 + (A2 - A1) / (1 + tmp);

	double aProb = h_inf / tau_h;
	double bProb = (1 - h_inf) / tau_h;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_m_egl19(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::egl19_mInf_v0);
	double dv = params.get(GeneType::egl19_mInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double m_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::egl19_mTau_A1);
	A2 = params.get(GeneType::egl19_mTau_A2);
	v0 = params.get(GeneType::egl19_mTau_v0);
	dv = params.get(GeneType::egl19_mTau_dv);
	tmp = std::exp(-(v - v0) / dv);
	double tau_m = A1 + (A2 - A1) / (1 + tmp);

	double aProb = m_inf / tau_m;
	double bProb = (1 - m_inf) / tau_m;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_f_egl19(const Params& params, double V, double ca_shell05_node) {
	double A1 = params.get(GeneType::egl19_fInf_A1);
	double A2 = 1.0;
	double ca05 = params.get(GeneType::egl19_fInf_ca05);
	double dca = params.get(GeneType::egl19_fInf_dca);
	double tmp = std::exp(-(ca_shell05_node - ca05) / dca);
	double f_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::egl19_fTau_A1);
	A2 = params.get(GeneType::egl19_fTau_A2);
	ca05 = params.get(GeneType::egl19_fTau_ca05);
	dca = params.get(GeneType::egl19_fTau_dca);
	tmp = std::exp(-(ca_shell05_node - ca05) / dca);
	double tau_f = A1 + (A2 - A1) / (1 + tmp);

	double aProb = f_inf / tau_f;
	double bProb = (1 - f_inf) / tau_f;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_m_exp2(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::exp2_mInf_v0);
	double dv = params.get(GeneType::exp2_mInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double m_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::exp2_mTau_act_A1);
	A2 = params.get(GeneType::exp2_mTau_act_A2);
	v0 = params.get(GeneType::exp2_mTau_act_v0);
	dv = params.get(GeneType::exp2_mTau_act_dv);
	double pow = params.get(GeneType::exp2_mTau_act_pow);
	double tau_m = A1 + A2 * std::exp(dv * std::pow(std::abs(v - v0),pow));

	double aProb = m_inf / tau_m;
	double bProb = (1 - m_inf) / tau_m;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_h_exp2(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::exp2_hInf_v0);
	double dv = params.get(GeneType::exp2_hInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double h_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::exp2_hTau_inact_A1);
	A2 = params.get(GeneType::exp2_hTau_inact_A2);
	v0 = params.get(GeneType::exp2_hTau_inact_v0);
	dv = params.get(GeneType::exp2_hTau_inact_dv);
	tmp = std::exp(-(v - v0) / dv);
	double tau_h = A1 + (A2 - A1) / (1 + tmp);

	double aProb = h_inf / tau_h;
	double bProb = (1 - h_inf) / tau_h;
	return {aProb, bProb};
}

double HHFunctions::calcC05(const Params& params, double VPrev, double VCurr, double hm205_cca1, double fm205_egl19) {
	double Pcca1 = params.get(GeneType::p_cca1_ca);		//[cm/s]
	double Pegl19 = params.get(GeneType::p_egl19_ca);	//[cm/s]
	double vPrev = 1E-3 * (VPrev + Params::VM);
	double vCurr = 1E-3 * (VCurr + Params::VM);
	double v = 0.5 * (3 * vCurr - vPrev);
	double tmp1 = std::exp(-Params::Z * Params::F * v / (Params::R * Params::T));
	double tmp2 = square(Params::Z) * square(Params::F) / (Params::R * Params::T);

	// Use first/second line below when only CCA/EGL CALCIUM current is used (in main function)
	//return ((Pcca1 * hm205_cca1) * tmp2 * v / (1 - tmp1));
	//return ((Pegl19 * fm205_egl19) * tmp2 * v / (1 - tmp1));
	return ((Pcca1 * hm205_cca1 + Pegl19 * fm205_egl19) * tmp2 * v / (1 - tmp1));
}

double HHFunctions::calcE05(double VPrev, double VCurr) {
	double vPrev = 1E-3 * (VPrev + Params::VM);
	double vCurr = 1E-3 * (VCurr + Params::VM);
	double v = 0.5 * (3 * vCurr - vPrev);

	return (std::exp(-Params::Z * Params::F * v / (Params::R * Params::T)));
}

double HHFunctions::calcGca(double VPrev, double VCurr, double CaIn05) {
	double vPrev = 1E-3 * (VPrev + Params::VM);		//v[i][timeStep-1]
	double vCurr = 1E-3 * (VCurr + Params::VM);		//v[i][timeStep]
	double v = 0.5 * (3 * vCurr - vPrev);			//v[i][timeStep+0.5]
	double tmp = std::exp(-Params::Z * Params::F * v / (Params::R * Params::T));

	return (square(Params::Z) * square(Params::F) * v / (Params::R * Params::T) * (CaIn05 - Params::CA_OUT * tmp) / (1 - tmp));		//use ca_shell[i][timeStep+0.5]
}

AlphaBeta HHFunctions::alphaBeta_m_cca1_HP(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = params.get(GeneType::cca1_mInf_A2);
	double v0 = params.get(GeneType::cca1_mInf_v0);
	double dv = params.get(GeneType::cca1_mInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double m_inf = A1 + (A2 - A1) / (1 + tmp);
	double tau_m = params.get(GeneType::cca1_mTau_deact);	//1.0;

	double aProb = m_inf / tau_m;
	double bProb = (1 - m_inf) / tau_m;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_m_exp2_HP(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::exp2_mInf_v0);
	double dv = params.get(GeneType::exp2_mInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double m_inf = A1 + (A2 - A1) / (1 + tmp);

	A1 = params.get(GeneType::exp2_mTau_deact_A1);
	A2 = params.get(GeneType::exp2_mTau_deact_A2);
	v0 = params.get(GeneType::exp2_mTau_deact_v0);
	dv = params.get(GeneType::exp2_mTau_deact_dv);
	tmp = std::exp(-(v - v0) / dv);
	double tau_m = A1 + (A2 - A1) / (1 + tmp);

	double aProb = m_inf / tau_m;
	double bProb = (1 - m_inf) / tau_m;
	return {aProb, bProb};
}

AlphaBeta HHFunctions::alphaBeta_h_exp2_HP(const Params& params, double V) {
	double A1 = 0.0;
	double A2 = 1.0;
	double v0 = params.get(GeneType::exp2_hInf_v0);
	double dv = params.get(GeneType::exp2_hInf_dv);
	double v = V + Params::VM;
	double tmp = std::exp(-(v - v0) / dv);
	double h_inf = A1 + (A2 - A1) / (1 + tmp);

	double tau_h = params.get(GeneType::exp2_hTau_deinact);

	double aProb = h_inf / tau_h;
	double bProb = (1 - h_inf) / tau_h;
	return {aProb, bProb};
}
