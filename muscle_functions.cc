#include "muscle_functions.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>
using namespace std;

void MuscleFunctions::calc_m_PED_CEK_PEK (const VectorN& v_seg_init, VectorN& m, VectorN& PE_D, VectorN& CE_K, VectorN& PE_K, boolVecN1* print_in_out) {
	//boolVecN1& print = *print_in_out;
	constexpr double MASS = Params::M;
	constexpr double PEK = Params::PE_K;
	constexpr double CEK = Params::CE_K;
	double PED_pm12_factor = 0.6, PED_pm3_factor = 1.0, PED_co_factor = 1.0, PED_isth_factor = 1.5, PED_TB_factor = 1.0;
	double CEK_pm12_factor = 1.0, CEK_pm3_factor = 1.0, CEK_co_factor = 1.0, CEK_isth_factor = 1.0, CEK_TB_factor = 0.08;
	double PED_factor_0, PED_factor_1, PED_factor_i, idx_50, slope, tmp, CEK_factor_0, CEK_factor_1, CEK_factor_i;

	for (int i = 0; i < Params::N; i++) {
		m[i] = MASS * v_seg_init[i] / 3;														//Divide by the number of muscle cells (or sarcomere groups) per segment, due to threefold radial symmetry
		PE_K[i] = PEK;
		if (i >= Params::PM6_START_NODE) {														// at the TB
			PE_D[i] = PED_TB_factor * std::sqrt(4 * m[i] * PE_K[i]);							// Critical damping ratio for no oscillations when returning to resting muscle length
			CE_K[i] = CEK_TB_factor * CEK;
		} else {
			if (i >= Params::ANT_ISTH_START_NODE + 3 * Params::NODE_NUM_TRANSITION_CO) {		// at the isthmus
				PE_D[i] = PED_isth_factor * std::sqrt(4 * m[i] * PE_K[i]);
				CE_K[i] = CEK_isth_factor * CEK;
			} else {
				if (i >= Params::ANT_ISTH_START_NODE) {											// at the pm4-5 border (PE_D)
					PED_factor_0 = PED_co_factor;
					PED_factor_1 = PED_isth_factor;
					idx_50 = Params::ANT_ISTH_START_NODE + 1.5 * Params::NODE_NUM_TRANSITION_CO;
					slope = 15;//20;//15;
					tmp = std::exp(- (i - idx_50) / slope);
					PED_factor_i = PED_factor_0 + (PED_factor_1 - PED_factor_0) / (1 + tmp);
					PE_D[i] = PED_factor_i * std::sqrt(4 * m[i] * PE_K[i]);
					CE_K[i] = CEK_isth_factor * CEK;
				} else {																		// at the corpus
					if (i >= Params::PM4_START_NODE) {											// at pm4
						if (i >= Params::PM4_PEAK) {											// at the pm4-5 border (CE_K)
							PE_D[i] = PED_co_factor * std::sqrt(4 * m[i] * PE_K[i]);
							CEK_factor_0 = CEK_co_factor;
							CEK_factor_1 = CEK_isth_factor;
							idx_50 = Params::ANT_ISTH_START_NODE - 6 * Params::NODE_NUM_TRANSITION_CO;
							slope = 40;
							tmp = std::exp(- (i - idx_50) / slope);
							CEK_factor_i = CEK_factor_0 + (CEK_factor_1 - CEK_factor_0) / (1 + tmp);
							CE_K[i] = CEK_factor_i * CEK;
						} else {																// at the pm3-4 border
							PED_factor_0 = PED_pm3_factor;
							PED_factor_1 = PED_co_factor;
							CEK_factor_0 = CEK_pm3_factor;
							CEK_factor_1 = CEK_co_factor;
							idx_50 = Params::PM4_START_NODE + 6 * Params::NODE_NUM_TRANSITION_CO;
							slope = 40;
							tmp = std::exp(- (i - idx_50) / slope);
							PED_factor_i = PED_factor_0 + (PED_factor_1 - PED_factor_0) / (1 + tmp);
							PE_D[i] = PED_factor_i * std::sqrt(4 * m[i] * PE_K[i]);
							CEK_factor_i = CEK_factor_0 + (CEK_factor_1 - CEK_factor_0) / (1 + tmp);
							CE_K[i] = CEK_factor_i * CEK;
						}
					} else {
						if (i < 9 * Params::PM3_START_NODE) {									// at the pm12-3 border
							PED_factor_0 = PED_pm12_factor;
							PED_factor_1 = PED_pm3_factor;
							CEK_factor_0 = CEK_pm12_factor;
							CEK_factor_1 = CEK_pm3_factor;
							idx_50 = 4 * Params::PM3_START_NODE;
							slope = 60;
							tmp = std::exp(- (i - idx_50) / slope);
							PED_factor_i = PED_factor_0 + (PED_factor_1 - PED_factor_0) / (1 + tmp);
							PE_D[i] = PED_factor_i * std::sqrt(4 * m[i] * PE_K[i]);
							CEK_factor_i = CEK_factor_0 + (CEK_factor_1 - CEK_factor_0) / (1 + tmp);
							CE_K[i] = CEK_factor_i * CEK;
						} else {																// at pm3
							PE_D[i] = PED_pm3_factor * std::sqrt(4 * m[i] * PE_K[i]);
							CE_K[i] = CEK_pm3_factor * CEK;
						}
					}
				}
			}

		}
		/*if (print[i]){
			cout << i << "\t" << CEK_factor_i << "\t" << PED_factor_i << std::endl;
			print[i] = false;
		}*/
	}
}

double MuscleFunctions::calcK1(int i, double ca_in, boolVecN1* print_in_out) {

	//boolVecN1& print = *print_in_out;
	double pm12_delta = 3.5E-5;
	double pm3_delta  = 0.3E-5;
	double pm4_delta  =	(Params::FLAT_META_CO == 0) ? - 1.4E-5 : pm3_delta;
	double isth_delta = 5E-5;
	double TB_delta   = - 4.5E-5;

	double delta_K1_B_0, delta_K1_B_1, idx_50, slope, tmp, delta_K1_B;

	if (i >= Params::PM6_START_NODE) {																		//at the TB
		return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + TB_delta) + quad(ca_in)));
	} else {
		if (i >= Params::ANT_ISTH_START_NODE) {																//at the isthmus
			return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + isth_delta) + quad(ca_in)));
		} else {
			if (i >= Params::PM4_START_NODE) {																//at the meta-corpus
				if (i >= Params::PM4_PEAK) {																//at the pm4-5 border
					delta_K1_B_0 = pm4_delta;
					delta_K1_B_1 = isth_delta;
					idx_50 = Params::ANT_ISTH_START_NODE - 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K1_B = delta_K1_B_0 + (delta_K1_B_1 - delta_K1_B_0) / (1 + tmp);
					return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + delta_K1_B) + quad(ca_in)));
				} else {																					//at the pm3-4 border
					delta_K1_B_0 = pm3_delta;
					delta_K1_B_1 = pm4_delta;
					idx_50 = Params::PM4_START_NODE + 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K1_B = delta_K1_B_0 + (delta_K1_B_1 - delta_K1_B_0) / (1 + tmp);
					return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + delta_K1_B) + quad(ca_in)));
				}
			} else {																						//at the pro-corpus
				if (i < 9 * Params::PM3_START_NODE) {														//at the pm12-3 border
					delta_K1_B_0 = pm12_delta;
					delta_K1_B_1 = pm3_delta;
					idx_50 = 4 * Params::PM3_START_NODE;
					slope = 60;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K1_B = delta_K1_B_0 + (delta_K1_B_1 - delta_K1_B_0) / (1 + tmp);
					return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + delta_K1_B) + quad(ca_in)));
				} else {																					//at pm3
					return (Params::K1_A * quad(ca_in) / (quad(Params::K1_B + pm3_delta) + quad(ca_in)));
				}
			}
		}
	}
}

double MuscleFunctions::calcMuscleContractionFactor(double t, double AP_t0, int i, boolVecN1* print_in_out, double x) {
	double t_CF_start = 0;
	double t_CF_start_isthmus = 0;
	double t_CF_start_TB = 0;
	double t_CF_start_shift_corpus = 0.01;
	double t_CF_start_shift_isthmus = 0.04;
	double t_CF_start_shift_TB = 0.01;
	double t_CF_stop = 200;
	double t_CF_stop_pm12 = -20;
	double t_CF_stop_pm3 = -6;
	double t_CF_stop_pm4 = (Params::FLAT_META_CO == 0) ? 0 : t_CF_stop_pm3;
	double t_CF_stop_shift_isthmus = 0;
	double t_CF_stop_isthmus = (Params::FLAT_META_CO == 0) ? -14 : -5;
	double t_CF_stop_TB = -96;
	double KCF_activation_level = 0.0;
	double KCF_activate_t05_pm12 = 107;
	double KCF_activate_slope_pm12 = 0.2;
	double KCF_activate_t05_pm3 = 105;
	double KCF_activate_slope_pm3 = 0.2;
	double KCF_activate_t05_pm4 = (Params::FLAT_META_CO == 0) ? 101 : KCF_activate_t05_pm3;
	double KCF_activate_slope_pm4 = 0.2;
	double KCF_activate_t05_isth = 123;
	double KCF_activate_slope_isth = 0.2;
	double KCF_activate_t05_TB = 99;
	double KCF_activate_slope_TB = 0.2;
	double KCF_activate_extend_dur_pm12 = 12;
	double KCF_activate_extend_slope_pm12 = 1;
	double KCF_activate_extend_t05_pm12 = 0.5;
	double KCF_activate_extend_dur_pm3 = (Params::FLAT_META_CO == 0) ? 29 : 22;
	double KCF_activate_extend_slope_pm3 = 1;
	double KCF_activate_extend_t05_pm3 = 0.6;
	double KCF_activate_extend_dur_pm4 = (Params::FLAT_META_CO == 0) ? 5 : KCF_activate_extend_dur_pm3;
	double KCF_activate_extend_slope_pm4 = 0.5;
	double KCF_activate_extend_t05_pm4 = 0.6;
	double KCF_activate_extend_dur_isth = 48;
	double KCF_activate_extend_slope_isth = 4;
	double KCF_activate_extend_t05_isth = 0.5;
	double KCF_activate_extend_dur_TB = 14;
	double KCF_activate_extend_slope_TB = 1;
	double KCF_activate_extend_t05_TB = 0.5;
	double t_CF_start_by_idx;

	//boolVecN1& print = *print_in_out;
	if (Params::ANT_ISTH_START_NODE <= i && i < Params::PM6_START_NODE) {									//at the isthmus
		double isthmus_stop = t_CF_stop + t_CF_stop_isthmus + t_CF_stop_shift_isthmus * (i - Params::PM6_START_NODE);
		double KCF_activate_extend_dur = KCF_activate_extend_dur_isth;
		double KCF_activate_t05 = KCF_activate_t05_isth;
		double KCF_activate_slope = KCF_activate_slope_isth;
		double KCF_activate_extend_t05 = KCF_activate_extend_t05_isth;
		double KCF_activate_extend_slope = KCF_activate_extend_slope_isth;
		double t_CF_stop_0, t_CF_stop_1, idx_50, slope, tmp,
			   KCF_activate_t05_0, KCF_activate_t05_1, KCF_activate_slope_0, KCF_activate_slope_1,
			   KCF_activate_extend_t05_0, KCF_activate_extend_t05_1,
			   KCF_activate_extend_dur_0, KCF_activate_extend_dur_1,
			   KCF_activate_extend_slope_0, KCF_activate_extend_slope_1;
		if (i < Params::ANT_ISTH_START_NODE + 3 * Params::NODE_NUM_TRANSITION_CO) {							//at the pm4-5 border
			t_CF_stop_0 = t_CF_stop + t_CF_stop_pm4;
			t_CF_stop_1 = t_CF_stop + t_CF_stop_isthmus + t_CF_stop_shift_isthmus * (Params::ANT_ISTH_START_NODE + 3 * Params::NODE_NUM_TRANSITION_CO - Params::PM6_START_NODE);
			KCF_activate_extend_dur_0 = KCF_activate_extend_dur_pm4;
			KCF_activate_extend_dur_1 = KCF_activate_extend_dur_isth;
			KCF_activate_t05_0 = KCF_activate_t05_pm4;
			KCF_activate_t05_1 = KCF_activate_t05_isth;
			KCF_activate_slope_0 = KCF_activate_slope_pm4;
			KCF_activate_slope_1 = KCF_activate_slope_isth;
			KCF_activate_extend_t05_0 = KCF_activate_extend_t05_pm4;
			KCF_activate_extend_t05_1 = KCF_activate_extend_t05_isth;
			KCF_activate_extend_slope_0 = KCF_activate_extend_slope_pm4;
			KCF_activate_extend_slope_1 = KCF_activate_extend_slope_isth;
			idx_50 = Params::ANT_ISTH_START_NODE + 1.5 * Params::NODE_NUM_TRANSITION_CO;
			slope = 15;
			tmp = std::exp(- (i - idx_50) / slope);
			isthmus_stop = t_CF_stop_0 + (t_CF_stop_1 - t_CF_stop_0) / (1 + tmp);
			KCF_activate_extend_dur = KCF_activate_extend_dur_0 + (KCF_activate_extend_dur_1 - KCF_activate_extend_dur_0) / (1 + tmp);
			KCF_activate_t05 = KCF_activate_t05_0 + (KCF_activate_t05_1 - KCF_activate_t05_0) / (1 + tmp);
			KCF_activate_slope = KCF_activate_slope_0 + (KCF_activate_slope_1 - KCF_activate_slope_0) / (1 + tmp);
			KCF_activate_extend_t05 = KCF_activate_extend_t05_0 + (KCF_activate_extend_t05_1 - KCF_activate_extend_t05_0) / (1 + tmp);
			KCF_activate_extend_slope = KCF_activate_extend_slope_0 + (KCF_activate_extend_slope_1 - KCF_activate_extend_slope_0) / (1 + tmp);
		}

		if (t < AP_t0 + isthmus_stop + KCF_activate_extend_dur) {
			t_CF_start_by_idx = t_CF_start + t_CF_start_isthmus + t_CF_start_shift_isthmus * (i - Params::ANT_ISTH_START_NODE);
			if (t < AP_t0 + isthmus_stop) {
				KCF_activation_level = 1 / (1 + std::exp(- KCF_activate_slope * (t - (AP_t0 + t_CF_start_by_idx + KCF_activate_t05))));
			} else {
				double KCF_start = 1;
				double KCF_end   = 0;
				double t_50 = AP_t0 + isthmus_stop + KCF_activate_extend_t05 * KCF_activate_extend_dur;
				double slope = KCF_activate_extend_slope;
				double tmp = std::exp(- (t - t_50) / slope);
				KCF_activation_level = KCF_start + (KCF_end - KCF_start) / (1 + tmp);
			}

			/*if (print[i]){
				//cout << i << "\t" << t_RF_start_by_idx << "\t" << corpus_stop << "\t" << corpus_stop + KRF_generate_extend_dur_pm12 << std::endl;
				cout << x * 1E4 << "\t" << t_CF_start_by_idx << "\t" << isthmus_stop << "\t" << isthmus_stop + KCF_activate_extend_t05 * KCF_activate_extend_dur << std::endl;
				print[i] = false;
			}*/

		}
	} else {
		if (i < Params::ANT_ISTH_START_NODE) {																//at the corpus
			double t_CF_stop_0, t_CF_stop_1, idx_50, slope, tmp, KCF_activate_t05_0, KCF_activate_t05_1,
				   KCF_activate_slope_0, KCF_activate_slope_1, KCF_activate_extend_t05_0, KCF_activate_extend_t05_1,
				   KCF_activate_extend_dur_0, KCF_activate_extend_dur_1, KCF_activate_extend_slope_0, KCF_activate_extend_slope_1;
			if (i >= Params::PM4_START_NODE) {
				if (i >= Params::PM4_PEAK) {																//at the pm4-5 border
					t_CF_stop_0 = t_CF_stop + t_CF_stop_pm4;
					t_CF_stop_1 = t_CF_stop + t_CF_stop_pm4;
					KCF_activate_t05_0 = KCF_activate_t05_pm4;
					KCF_activate_t05_1 = KCF_activate_t05_pm4;
					KCF_activate_slope_0 = KCF_activate_slope_pm4;
					KCF_activate_slope_1 = KCF_activate_slope_pm4;
					KCF_activate_extend_dur_0 = KCF_activate_extend_dur_pm4;
					KCF_activate_extend_dur_1 = KCF_activate_extend_dur_pm4;
					KCF_activate_extend_t05_0 = KCF_activate_extend_t05_pm4;
					KCF_activate_extend_t05_1 = KCF_activate_extend_t05_pm4;
					KCF_activate_extend_slope_0 = KCF_activate_extend_slope_pm4;
					KCF_activate_extend_slope_1 = KCF_activate_extend_slope_pm4;
					idx_50 = Params::ANT_ISTH_START_NODE - 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
				} else {																					//at the pm3-4 border
					t_CF_stop_0 = t_CF_stop + t_CF_stop_pm3;
					t_CF_stop_1 = t_CF_stop + t_CF_stop_pm4;
					KCF_activate_t05_0 = KCF_activate_t05_pm3;
					KCF_activate_t05_1 = KCF_activate_t05_pm4;
					KCF_activate_slope_0 = KCF_activate_slope_pm3;
					KCF_activate_slope_1 = KCF_activate_slope_pm4;
					KCF_activate_extend_dur_0 = KCF_activate_extend_dur_pm3;
					KCF_activate_extend_dur_1 = KCF_activate_extend_dur_pm4;
					KCF_activate_extend_t05_0 = KCF_activate_extend_t05_pm3;
					KCF_activate_extend_t05_1 = KCF_activate_extend_t05_pm4;
					KCF_activate_extend_slope_0 = KCF_activate_extend_slope_pm3;
					KCF_activate_extend_slope_1 = KCF_activate_extend_slope_pm4;
					idx_50 = Params::PM4_START_NODE + 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
				}
			} else {
				if (i < 9 * Params::PM3_START_NODE) {														//at the pm12-3 border
					t_CF_stop_0 = t_CF_stop + t_CF_stop_pm12;
					t_CF_stop_1 = t_CF_stop + t_CF_stop_pm3;
					KCF_activate_t05_0 = KCF_activate_t05_pm12;
					KCF_activate_t05_1 = KCF_activate_t05_pm3;
					KCF_activate_slope_0 = KCF_activate_slope_pm12;
					KCF_activate_slope_1 = KCF_activate_slope_pm3;
					KCF_activate_extend_dur_0 = KCF_activate_extend_dur_pm12;
					KCF_activate_extend_dur_1 = KCF_activate_extend_dur_pm3;
					KCF_activate_extend_t05_0 = KCF_activate_extend_t05_pm12;
					KCF_activate_extend_t05_1 = KCF_activate_extend_t05_pm3;
					KCF_activate_extend_slope_0 = KCF_activate_extend_slope_pm12;
					KCF_activate_extend_slope_1 = KCF_activate_extend_slope_pm3;
					idx_50 = 4 * Params::PM3_START_NODE;
					slope = 60;
				} else {																					//at pm3
					t_CF_stop_0 = t_CF_stop + t_CF_stop_pm3;
					t_CF_stop_1 = t_CF_stop + t_CF_stop_pm3;
					KCF_activate_t05_0 = KCF_activate_t05_pm3;
					KCF_activate_t05_1 = KCF_activate_t05_pm3;
					KCF_activate_slope_0 = KCF_activate_slope_pm3;
					KCF_activate_slope_1 = KCF_activate_slope_pm3;
					KCF_activate_extend_dur_0 = KCF_activate_extend_dur_pm3;
					KCF_activate_extend_dur_1 = KCF_activate_extend_dur_pm3;
					KCF_activate_extend_t05_0 = KCF_activate_extend_t05_pm3;
					KCF_activate_extend_t05_1 = KCF_activate_extend_t05_pm3;
					KCF_activate_extend_slope_0 = KCF_activate_extend_slope_pm3;
					KCF_activate_extend_slope_1 = KCF_activate_extend_slope_pm3;
					idx_50 = Params::PM3_START_NODE + 0.5 * (Params::PM4_START_NODE - Params::PM3_START_NODE);
					slope = 10;
				}
			}
			tmp = std::exp(- (i - idx_50) / slope);
			double corpus_stop = t_CF_stop_0 + (t_CF_stop_1 - t_CF_stop_0) / (1 + tmp);
			double KCF_activate_t05 = KCF_activate_t05_0 + (KCF_activate_t05_1 - KCF_activate_t05_0) / (1 + tmp);
			double KCF_activate_slope = KCF_activate_slope_0 + (KCF_activate_slope_1 - KCF_activate_slope_0) / (1 + tmp);
			double KCF_activate_extend_dur = KCF_activate_extend_dur_0 + (KCF_activate_extend_dur_1 - KCF_activate_extend_dur_0) / (1 + tmp);
			double KCF_activate_extend_t05 = KCF_activate_extend_t05_0 + (KCF_activate_extend_t05_1 - KCF_activate_extend_t05_0) / (1 + tmp);
			double KCF_activate_extend_slope = KCF_activate_extend_slope_0 + (KCF_activate_extend_slope_1 - KCF_activate_extend_slope_0) / (1 + tmp);

			if (t < AP_t0 + corpus_stop + KCF_activate_extend_dur) {
				t_CF_start_by_idx = t_CF_start + t_CF_start_shift_corpus * (Params::ANT_ISTH_START_NODE - i);
				if (t < AP_t0 + corpus_stop) {
					KCF_activation_level = 1 / (1 + std::exp(- KCF_activate_slope * (t - (AP_t0 + t_CF_start_by_idx + KCF_activate_t05))));
				} else {
					double KCF_start = 1;
					double KCF_end   = 0;
					double t_50 = AP_t0 + corpus_stop + KCF_activate_extend_t05 * KCF_activate_extend_dur;
					double slope = KCF_activate_extend_slope;
					double tmp = std::exp(- (t - t_50) / slope);
					KCF_activation_level = KCF_start + (KCF_end - KCF_start) / (1 + tmp);
				}

				/*if (print[i]){
					////cout << i << "\t" << t_RF_start_by_idx << "\t" << corpus_stop << "\t" << corpus_stop + KRF_generate_extend_dur_pm12 << std::endl;
					cout << x * 1E4 << "\t" << t_CF_start_by_idx << "\t" << corpus_stop << "\t" << corpus_stop + KCF_activate_extend_t05 * KCF_activate_extend_dur << std::endl;
					print[i] = false;
				}*/

			}
		} else {		// (Params::PM6_START_NODE <= i)													//at the TB
			double TB_stop = t_CF_stop + t_CF_stop_TB;
			if (t < AP_t0 + TB_stop + KCF_activate_extend_dur_TB) {
				t_CF_start_by_idx = t_CF_start + t_CF_start_shift_isthmus * (Params::PM6_START_NODE - Params::ANT_ISTH_START_NODE) + t_CF_start_TB + t_CF_start_shift_TB * (i - Params::PM6_START_NODE);
				if (t < AP_t0 + TB_stop) {
					KCF_activation_level = 1 / (1 + std::exp(- KCF_activate_slope_TB * (t - (AP_t0 + t_CF_start_by_idx + KCF_activate_t05_TB))));
				} else {
					double KCF_start = 1;
					double KCF_end   = 0;
					double t_50 = AP_t0 + TB_stop + KCF_activate_extend_t05_TB * KCF_activate_extend_dur_TB;
					double slope = KCF_activate_extend_slope_TB;
					double tmp = std::exp(- (t - t_50) / slope);
					KCF_activation_level = KCF_start + (KCF_end - KCF_start) / (1 + tmp);
				}

				/*if (print[i]){
					//cout << i << "\t" << t_RF_start_by_idx << "\t" << TB_stop << "\t" << TB_stop + KRF_generate_extend_dur_TB << std::endl;
					cout << x * 1E4 << "\t" << t_CF_start_by_idx << "\t" << TB_stop << "\t" << TB_stop + KCF_activate_extend_t05_TB * KCF_activate_extend_dur_TB << std::endl;
					print[i] = false;
				}*/
			}
		}
	}
	return KCF_activation_level;
}

double MuscleFunctions::calcMLCP_P(int i, double ca_in, double CF, double MLCP_P_prev, boolVecN1* print_in_out) {

	//boolVecN1& print = *print_in_out;
	double pm12_delta = 3.5E-5;
	double pm3_delta  = 0.3E-5;
	double pm4_delta  =	(Params::FLAT_META_CO == 0) ? - 1.4E-5 : pm3_delta;
	double isth_delta = 5E-5;
	double TB_delta   = - 4.5E-5;

	double delta_K2_ON3_0, delta_K2_ON3_1, idx_50, slope, tmp, delta_K2_ON3;

	if (i >= Params::PM6_START_NODE) {																	//at the TB
		delta_K2_ON3 = TB_delta;
	} else {
		if (i >= Params::ANT_ISTH_START_NODE) {															//at the isthmus
			delta_K2_ON3 = isth_delta;
		} else {
			if (i >= Params::PM4_START_NODE) {
				if (i >= Params::PM4_PEAK) {															//at the pm4-5 border
					delta_K2_ON3_0 = pm4_delta;
					delta_K2_ON3_1 = isth_delta;
					idx_50 = Params::ANT_ISTH_START_NODE - 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K2_ON3 = delta_K2_ON3_0 + (delta_K2_ON3_1 - delta_K2_ON3_0) / (1 + tmp);
				} else {																				//at the pm3-4 border
					delta_K2_ON3_0 = pm3_delta;
					delta_K2_ON3_1 = pm4_delta;
					idx_50 = Params::PM4_START_NODE + 6 * Params::NODE_NUM_TRANSITION_CO;
					slope = 40;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K2_ON3 = delta_K2_ON3_0 + (delta_K2_ON3_1 - delta_K2_ON3_0) / (1 + tmp);
				}
			} else {
				if (i < 9 * Params::PM3_START_NODE) {													//at the pm12-3 border
					delta_K2_ON3_0 = pm12_delta;
					delta_K2_ON3_1 = pm3_delta;
					idx_50 = 4 * Params::PM3_START_NODE;
					slope = 60;
					tmp = std::exp(- (i - idx_50) / slope);
					delta_K2_ON3 = delta_K2_ON3_0 + (delta_K2_ON3_1 - delta_K2_ON3_0) / (1 + tmp);
				} else {																				//at pm3
					delta_K2_ON3 = pm3_delta;
				}
			}
		}
	}
	double K2_on = Params::K2_ON1 + Params::K2_ON2 * square(ca_in) / (square(Params::K2_ON3 + delta_K2_ON3) + square(ca_in));
	double K2_off = Params::K2_OFF1 + CF / (Params::K2_OFF3 + CF);
	double DT_SEC = 1E-3 * Params::DT;
	return (MLCP_P_prev + DT_SEC / Params::MLCP_TAU * (K2_on * (1 - MLCP_P_prev) - K2_off * MLCP_P_prev));
}

double MuscleFunctions::calcMuscleActivation(double L0, double l) {

	double L0_um = L0 * 1E4;		//Convert L0 from [cm] to [um]
	double l_um = l * 1E4;			//Convert l  from [cm] to [um]
	double y2_val = 0.4;
	double y34_val = 1.0;
	double y15_val = 0.0;
	double x12_distance =  0.433;
	double y23_slope	=  0.31;
	double x34_distance =  0.17;
	double y45_slope	= -0.7;
	double m, n, x1, x2, y1, y2;

	x2 = L0_um + (y2_val - y34_val - y23_slope * x34_distance) / y23_slope - x12_distance; y2 = y15_val;
	if (l_um <= x2) {
		return y15_val;
	}
	x1 = x2; x2 += x12_distance;
	y1 = y2; y2 = y2_val;
	if (x1 < l_um && l_um <= x2) {
		m = (y2 - y1) / (x2 - x1); n = y2 - m * x2;
		return (Params::CE_ACT_MAX * (m * l_um + n));
	}
	x1 = x2; x2 -= (y2_val - y34_val) / y23_slope;
	y1 = y2; y2 = y34_val;
	if (x1 < l_um && l_um <= x2) {
		m = y23_slope; n = y2 - m * x2;
		return (Params::CE_ACT_MAX * (m * l_um + n));
	}
	x1 = x2; x2 += x34_distance;
	y1 = y2; y2 = y34_val;
	if (x1 < l_um && l_um <= x2) {
		return  Params::CE_ACT_MAX;
	}
	x1 = x2; x2 += (y15_val - y34_val) / y45_slope;
	y1 = y2; y2 = y15_val;
	if (x1 < l_um && l_um <= x2) {
		m = y45_slope; n = y2 - m * x2;
		return (Params::CE_ACT_MAX * (m * l_um + n));
	}
	return y15_val;	//if (x2 < l_um)
}

void MuscleFunctions::buildDistributionsFromNdistMoments(int seg, const MomentsVec& Q, Distribution& n, Distribution& n0l) {
	//Assuming n is between [-0.5C*H, 0.5C*H]
	constexpr double LtoH = Params::L_MYOSIN_BINDING_SITES / Params::H;
	constexpr int n_LtoH_length = 2 * Params::MAXIMAL_NORMALIZED_BOND_LENGTH;
	constexpr int n_points_num = Params::CE_POINTS;
	int n0l_points_num = n_points_num / n_LtoH_length;
	if (Q[0] == 0) {
		double xi_min = 0.0;
		double xi_max = 1.0;
		n.mean = 0.5 * (xi_min + xi_max);							// assume n is a uniform distribution between x[0,H] (<=> xi[0,1]) at rest
		n.standard_deviation = std::sqrt(1.0 / 12.0) * (xi_max - xi_min);
		n0l.mean = n.mean;
		n0l.standard_deviation = n.standard_deviation;
	} else {
		n.mean = Q[1] / Q[0];
		n.standard_deviation = std::sqrt(Q[2] / Q[0] - square(Q[1] / Q[0]));

		if (std::isnan(Q[1] / Q[0])) {cout << "Q[1] / Q[0] is nan at segment " << seg << std::endl; exit(EXIT_FAILURE);}
		if (Q[2] / Q[0] < square(Q[1] / Q[0])) {cout << "Q[2] / Q[0] < square(Q[1] / Q[0]) at segment " << seg << std::endl; exit(EXIT_FAILURE);}	//"Cannot calculate n0l_dist.SD. Check if M decreases enough (in/decrease K1/2)" << std::endl;}

		double n_vals[n_points_num];
		double n0l_vals[n0l_points_num];
		double xi;													// x along the sarcomere, normalized by H
		for (int i = 0; i < n_points_num; i++) {
			xi = n_LtoH_length * ((double)i / n_points_num - 0.5) * LtoH;
			n_vals[i] = Q[0] / (std::sqrt(2 * M_PI) * n.standard_deviation) * std::exp(- 0.5 * square((xi - n.mean) / n.standard_deviation));
		}
		double n0l_total = 0.0;
		n0l.mean = 0.0;
		int idx;
		for (int i = 0; i < n0l_points_num; i++) {
			n0l_vals[i] = 0.0;
			for (int LtoH_seg = 0; LtoH_seg < n_LtoH_length; LtoH_seg++) {
				idx = i + LtoH_seg * n0l_points_num;
				n0l_vals[i] += n_vals[idx];
			}
			xi = (double)i / n0l_points_num * LtoH;
			n0l.mean += n0l_vals[i] * xi;
			n0l_total += n0l_vals[i];
		}
		if (n0l_total == 0.0) {cout << "Cannot calculate n0l_dist.mean (n0l_total = 0). Increase DISCRETE_SEGMENTS_NUM" << std::endl; exit(EXIT_FAILURE);}
		n0l.mean = n0l.mean / n0l_total;
		n0l.standard_deviation = 0.0;
		for (int i = 0; i < n0l_points_num; i++) {
			xi = (double)i / n0l_points_num * LtoH;
			n0l.standard_deviation += square(xi - n0l.mean);
		}
		n0l.standard_deviation = std::sqrt(n0l.standard_deviation / n0l_points_num);
	}
}

double MuscleFunctions::calcBeta(int moment, double M, double Mp) {
	if (M + Mp == 0){return 0.0;}
	else {return (Mp / (M + Mp) * Params::K3_F1 / (moment + 2));}
}

double MuscleFunctions::calcMusclePhi1(int moment, double Q0, Distribution n0l, double M, double Mp) {
	if (M + Mp == 0){return 0.0;}
	else {
		double G1 = G(moment + 1, ((1 - n0l.mean) / n0l.standard_deviation), n0l);
		double G2 = G(moment + 1, (- n0l.mean / n0l.standard_deviation), n0l);
		return (Mp / (M + Mp) * Q0 * Params::K3_F1 * (G1 - G2));
	}
}

void MuscleFunctions::calcMusclePhi2(int moment, double Q0, Distribution n, double AM, double AMp, double* phi2s) {
	double phi2_1 = 0.0, phi2_2 = 0.0;
	if (AM + AMp > 0) {
		double G1 = G(moment, (- n.mean / n.standard_deviation), n);
		double G2 = G(moment + 1, ((1 - n.mean) / n.standard_deviation), n);
		double G3 = G(moment + 1, (- n.mean / n.standard_deviation), n);
		double G4 = G(moment, INFINITY, n);
		double G5 = G(moment, ((1 - n.mean) / n.standard_deviation), n);
		double G6 = G(moment + 1, INFINITY, n);
		double G7 = G(moment + 1, ((1 - n.mean) / n.standard_deviation), n);
		phi2_1 = AMp / (AM + AMp) * Q0 * Params::K4_G1 * (Params::K_C1 * G1 + (G2 - G3) + (Params::K_D1 - Params::K_D2) * (G4 - G5) + Params::K_D2 * (G6 - G7));
		phi2_2 = AM  / (AM + AMp) * Q0 * Params::K7_G2 * (Params::K_C2 * G1 + (G2 - G3) + (Params::K_D3 - Params::K_D4) * (G4 - G5) + Params::K_D4 * (G6 - G7));
	}
	phi2s[0] = phi2_1;
	phi2s[1] = phi2_2;
}

double MuscleFunctions::calcF_OE(int i, const VectorN& l, const VectorN& segR) {
	double F_OE_left = 0.0, F_OE_right = 0.0;
	double left_z_muscle  = segR[i - 1] - l[i - 1];				//z_muscle[i-1] at (timeStep)
	double z_muscle		  = segR[i]     - l[i];					//z_muscle[i]   at (timeStep)
	double right_z_muscle = segR[i + 1] - l[i + 1];				//z_muscle[i+1] at (timeStep)
	if (i > 0)			   {F_OE_left  = Params::OE_K * (left_z_muscle  - z_muscle);}
	if (i < Params::N - 1) {F_OE_right = Params::OE_K * (right_z_muscle - z_muscle);}
	return (F_OE_left + F_OE_right);
}

double MuscleFunctions::G(int k, double z, Distribution dist) {
	double p = dist.mean, q = dist.standard_deviation, result = 0;
	for (int m = 0; m <= k; m++) {
		result += factorial(k) / (factorial(m) * factorial(k - m)) * std::pow(p, k - m) * std::pow(q, m) * F(m, z);
	}
	return result;
}

double MuscleFunctions::F(int m, double z) {
	if (m == 0) {
		return (0.5 * (1 + std::erf(z / std::sqrt(2))));
	}
	if (m == 1) {
		return (- 1 / std::sqrt(2 * M_PI) * std::exp(- 0.5 * square(z)));
	}
	// if (m > 1):
	if (std::isinf(z)) {return ((m - 1) * F(m - 2, z));}		// Need to handle z = INFINITY differently, since in c++ INFINITY * 0 = NaN instead of 0
	else {return ((m - 1) * F(m - 2, z) - std::pow(z, m - 1) / std::sqrt(2 * M_PI) * std::exp(- 0.5 * square(z)));}
}

int MuscleFunctions::factorial(int n) {
	return ((n == 1 || n == 0) ? 1 : factorial(n - 1) * n);
}
