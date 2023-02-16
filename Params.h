#ifndef PARAMS_H_
#define PARAMS_H_

#include "base.h"

class Params {
 public:
	static constexpr int	EQUAL_VOLUME_ELEMENTS = 0;									//0 = equal-length elements; 1 = equal-volume elements
	static constexpr int	FLAT_META_CO = 1;//0;											//0 = spherical; 1 = flat
	static constexpr int	XY_INPUT_LENGTH =  3661;
	static constexpr int	XY_OUTPUT_LENGTH = (EQUAL_VOLUME_ELEMENTS == 0) ? 3660 : 1587;	//change the 1587
	static constexpr int	N = XY_OUTPUT_LENGTH;									 	//Number of nodes along the x-axis (original value: 3000)
	static constexpr int	PM3_START_NODE = 85;
	static constexpr int	PM4_START_NODE = 1282;
	static constexpr int	PM4_PEAK = 1616;
	static constexpr int	ANT_ISTH_START_NODE = 2000;
	static constexpr int	POST_ISTH_START_NODE = 2621; 								//(x=122.67) when setting the isthmus to lie between x=93.6-140, and assuming anterior isthmus length is ~2/3 of entire isthmus
	static constexpr int	PM6_START_NODE = 2991;
	static constexpr int	TB_PEAK = 3400;
	static constexpr int	EX_INPUT_NODE  = 2000;										//Node's index at which exogenous input current is injected
																						//Set to be at the metacorpus-isthmus border
	// For anterior-to-posterior AP along the pharynx:
	//static constexpr int	EX_INPUT_NODE  = 0;											//0 = the anterior-most point
	static constexpr int	NODE_NUM_TRANSITION_CO = 30;								//Number of segments for gradual opening/closing at the border between different pharyngeal areas
	static constexpr int	NUM_GENES = 51;												//For genetic algorithm
	static constexpr double CM = 1.0;													//Specific membrane capacitance [mu F/cm^2]
	static constexpr double GA = 0.5;													//Specific intracellular conductance [m Siemens/cm]
	static constexpr double R = 8.68815102;												//Gas constant [J/K*mol]
	static constexpr double T = 295.65;													//Temperature [K] in Jospin et al. 2002 experiments (i.e., most suitable for the EGL19 channel's parameters)
	static constexpr double Z = 2;														//Valency of calcium ions []
	static constexpr double F = 96485;													//Faraday constant [C/mol]
	static constexpr double CA_IN_REST = 5E-8;											//Intracellular calcium concentration at rest [mol/L]
	static constexpr double CA_OUT = 0.003;												//Extracellular calcium concentration [mol/L] in Shtonda & Avery 2005 experiments
	static constexpr double CA_DIFFUSE_COEF = 3.5E-6;									//Calcium diffusion coefficient in muscle's cytoplasm [cm^2/sec]
	static constexpr double SHELL_THICKNESS = 0.02;										//Thickness of the segment's outer-most shell (2% of segment's radius)
	static constexpr double VM = -73.0;													//Equilibrium transmembrane voltage [mV]
	static constexpr double V_NA = 69.0;												//Reversal potential for sodium [mV]; calculated for Na out = 150, Na in = 10 mM
	static constexpr double V_K = -82.0;												//Reversal potential for potassium [mV]; calculated for K out = 6, K in = 150 mM
	static constexpr double V_L = -73.0;												//Leakage reversal potential [mV] (original value: -49.402)
	static constexpr double EX_INPUT_AMP = 62;											//Injected current [mu amp/cm]				//Injected current [mu amp]	(original value: -5.0)	//-8000.0 for Na alone with this geometry
	static constexpr double EX_INPUT_DURATION = 40.0;									//[ms], based on Steger et al. 2005 Fig. 4A (down)
	static constexpr double H = 1.56E-6;												//The range of bond lengths over which myosin has a nonzero bonding rate [cm]
	static constexpr double CE_ACT_MAX = 0.7;											//Maximal fraction of myosin heads available for cross-bridges in the CE
	static constexpr double L_MYOSIN_BINDING_SITES = H / CE_ACT_MAX;					//Distance between successive binding sites for myosin along the actin filaments [cm]
	static constexpr int	MAXIMAL_NORMALIZED_BOND_LENGTH = 10;						//Maximal myosin-actin bond length, expressed by multiples of normalized L_ACTIN_BINDING_SITES, i.e., of L_ACTIN_BINDING_SITES / H.
	static constexpr int	CE_POINTS = MAXIMAL_NORMALIZED_BOND_LENGTH * 2 * 10;		//Number of discrete points of the n (AM + AMp) distribution along the z-axis, from which the n' distribution is calculated.
	static constexpr double K1_A = 1000.0;												//Maximal activity (phosphorylation) rate of MLCK [1/s]
	static constexpr double K1_B = 5.2E-5;												//Dissociation constant of MLCK from calcium [mol/L]
	static constexpr double MLCP_TAU = 1E-4;											//MLCP time constant [s]
	static constexpr double K2_ON1 = 0.0;												//Basal   activation rate of MLCP (unit-less)
	static constexpr double K2_ON2 = 0.1;												//Maximal activation rate of MLCP (unit-less)
	static constexpr double K2_ON3 = 5.2E-5;											//Dissociation constant of MLCP from calcium [mol/L]
	static constexpr double K2_OFF1 = 0.0;												//Basal   inhibition rate of MLCP (unit-less)
	static constexpr double K2_OFF3 = 1.0;												//Maximal inhibition rate of MLCP (unit-less)
	static constexpr double K2_MAX = 10000.0;											//Maximal activity (de-phosphorylation) rate of MLCP [1/s]
	static constexpr double K4_AVG = 1000.0;											//Average unbonding rate of AMp [1/s]
	static constexpr double K4_G1 = 1.5 * K4_AVG;										//Maximal unbonding rate of AMp [1/s]
	static constexpr double K7_AVG = 0.3 * K4_AVG;										//Average unbonding rate of AM [1/s]
	static constexpr double K7_G2 = 1.5 * K7_AVG;										//Maximal unbonding rate of AM [1/s]
	static constexpr double K3_AVG = 0.4 * K4_AVG;										//Average bonding rate of Mp to A [1/s]
	static constexpr double K3_F1 = 2 * K3_AVG;											//Maximal bounding rate of Mp to A [1/s]
	static constexpr double K_D1 = 2.0;													//Assuming d1=d2=d3=d4 (unit-less)
	static constexpr double K_D2 = K_D1;
	static constexpr double K_D3 = K_D1;
	static constexpr double K_D4 = K_D1;
	static constexpr double K_C1 = 10 * K_D1;											//Assuming c1=c2 (unit-less)
	static constexpr double K_C2 = K_C1;
	static constexpr double M = 1.074;													//Muscle's mass density [g/cm3]
	static constexpr int	CE_NUM_AM = 3E3;											//Number of myosin heads per HALF sarcomere group (unit-less)
	static constexpr double CE_K = 2.2E-7;												//Active element's spring constant [g/s^2]
	static constexpr double PE_K = 1.1E-7;												//Passive element's spring constant [g/s^2]
	static constexpr double OE_K = 3E-6;												//Connective tissue's spring constant [g/s^2]
	static constexpr double DT = 1E-1; //1E-2;											//Time integration step [ms]
	static constexpr int    TEND = 300;	//500;											//Total duration of simulation [ms]
	static constexpr int	TIME_POINTS_SAMPLING_RATE = (int) (1 / DT);
	static constexpr int	SIMULATED_TIME_POINTS = (int) (TEND * TIME_POINTS_SAMPLING_RATE);
	static constexpr int	EX_INPUT_INTERVAL = (int) (1450 * TIME_POINTS_SAMPLING_RATE);	//based on Ca+2inDynamics_Kerr_Shimonozo.xlsx [ms]
																							//350 ms during high-rate pumping, and 1450 ms during low-rate pumping.

	inline double get(GeneType gt) const { return genes[static_cast<int>(gt)]; }

	Vector<NUM_GENES> genes;
};

using VectorN   = Vector<Params::N>;
using VectorN1  = Vector<Params::N + 1>;
using Vector2N  = Vector<2 * Params::N>;
using Vector2N1 = Vector<2 * Params::N + 1>;
using boolVecN1 = boolVec<Params::N + 1>;
using TriMatN   = Vector<3 * Params::N - 2>;
using TriMatN1  = Vector<3 * (Params::N + 1) - 2>;
using MomentsN  = std::array<MomentsVec, Params::N>;
using DistN     = std::array<Distribution, Params::N>;

#endif	// PARAMS_H_
