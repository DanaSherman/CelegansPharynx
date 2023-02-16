#include "gene.h"

const std::array<GeneRange, Params::NUM_GENES> GENE_RANGES = {{
		{   0.1,	20.0,	100.0,	  1.0},												// cca1_mInf_A2
		{-100.0,	20.0,	100.0,	-30.0},												// cca1_mInf_v0
		{   0.1,	 0.1,	 50.0,	  5.0},												// cca1_mInf_dv
		{   0.1,	20.0,	100.0,	  2.7},												// cca1_mTau_act_A1
		{   0.1,	20.0,	100.0,	 22.3},												// cca1_mTau_act_A2 > cca1_mTau_A1
		{-100.0,	20.0,	100.0,	-52.0},												// cca1_mTau_act_v0
		{ -50.0,	 0.1,	 -0.1,	 -8.0},												// cca1_mTau_act_dv
		{   0.1,	20.0,	100.0,	  1.0},												// cca1_mTau_deact
		{-100.0,	20.0,	100.0,	-76.0},												// cca1_hInf_v0
		{ -50.0,	 0.1,	 -0.1,	 -4.0},												// cca1_hInf_dv
		{   0.1,	20.0,	100.0,	  1.0},												// cca1_hTau_A1
		{   0.1,	20.0,	100.0,  250.0},												// cca1_hTau_A2 > cca1_hTau_A1
		{-100.0,	20.0,	100.0,	  0.0},												// cca1_hTau_v0
		{ -50.0,	 0.1,	 -0.1,	 -4.5},												// cca1_hTau_dv
		{-100.0,	 1.0,	100.0,	  5.6},												// egl19_mInf_v0
		{   0.1,	 0.1,	 50.0,	  7.5},												// egl19_mInf_dv
		{   0.1,	 1.0,	100.0,	  0.5455},											// egl19_mTau_A1
		{   0.1,	 1.0,	100.0,	  0.9635},											// egl19_mTau_A2 > egl19_mTau_A1
		{-100.0,	 1.0,	100.0,	  8.5},												// egl19_mTau_v0
		{ -50.0,	 0.1,	 -0.1,	  3},												// egl19_mTau_dv
		{   0.1,	20.0,	100.0,	  0.8},												// egl19_fInf_A1
		{   1E-6,	 1E-6,	  1E-3,	  1.3E-4},											// egl19_fInf_ca05
		{  -1E-7,	 1E-7,	 -1E-4,	 -5E-6},											// egl19_fInf_dca
		{   0.1,	 1.0,	100.0,	  5.0},												// egl19_fTau_A1
		{   0.1,	 1.0,	100.0,	 50.0},												// egl19_fTau_A2 > egl19_fTau_A1
		{   1E-6,	 1E-6,	  1E-3,	  1E-4},											// egl19_fTau_ca05
		{   1E-7,	 1E-7,	  1E-4,	  1E-5},											// egl19_fTau_dca
		{-100.0,	 1.0,	100.0,	-15.5},												// exp2_mInf_v0
		{   0.1,	 0.1,	 50.0,	  7.0},												// exp2_mInf_dv
		{   0.1,	 1.0,	100.0,	 15.0},												// exp2_mTau_act_A1
		{   0.1,	 1.0,	100.0,	730.0},												// exp2_mTau_act_A2 > exp2_mTau_act_A1
		{-100.0,	 1.0,	100.0,	-20.0},												// exp2_mTau_act_v0
		{   0.1,	 0.1,	 50.0,	 -0.15},											// exp2_mTau_act_dv
		{   1.65,	 0.05,	  2.0,	  0.8},												// exp2_mTau_act_pow
		{-100.0,	 1.0,	100.0,	 -5.0},												// exp2_hInf_v0
		{ -50.0,	 0.1,	 -0.1,	 -1.0},												// exp2_hInf_dv
		{   0.1,	 1.0,	100.0,	  0.3},												// exp2_hTau_inact_A1
		{   0.1,	 1.0,	100.0,	 18.0},												// exp2_hTau_inact_A2 > exp2_hTau_inact_A1
		{-100.0,	 1.0,	100.0,	-60.0},												// exp2_hTau_inact_v0
		{   0.1,	 0.1,	 50.0,	-20.0},												// exp2_hTau_inact_dv
		{   0.1,	 1.0,	100.0,	 65.0},												// exp2_mTau_deact_A1
		{   0.1,	 1.0,	100.0,	 75.2},												// exp2_mTau_deact_A2 > exp2_mTau_deact_A1
		{-100.0,	 1.0,	100.0,	-83.9},												// exp2_mTau_deact_v0
		{   0.1,	 0.1,	 50.0,	 -3.4},												// exp2_mTau_deact_dv
		{   0.1,	 1.0,	100.0,	  1.0},												// exp2_hTau_deinact
		{   0.0001,	 0.01, 	 10.0,	  0.009},											// p_cca1_ca
		{  10.0,	 0.01,	200.0,	  2.5},												// g_cca1_na
		{   0.0001,	 0.01,	 10.0,	  0.001},											// p_egl19_ca
		{  10.0,	0.05,	500.0,	 (Params::FLAT_META_CO == 0) ? 0.1523 : 0.1509},	// g_egl19_na
		{   0.01,	0.01,	200.0,	  0.8},												// g_exp2
		{  10.0,    1.0,	300.0,	 50.0},												// ca_in_tau [ms]
}};
