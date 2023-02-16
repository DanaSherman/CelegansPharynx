#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
using namespace std;

#include "base.h"
#include "build_axial_current_function.h"
#include "build_cable_function_space_integration.h"
#include "build_cable_function_time_integration.h"
#include "build_calcium_function_space_integration.h"
#include "build_calcium_function_time_integration.h"
#include "build_finite_element_mesh.h"
#include "build_lumen_function_space_integration.h"
#include "build_lumen_function_time_integration.h"
#include "flow_functions.h"
#include "gene.h"
#include "hh_functions.h"
#include "init_segments.h"
#include "muscle_functions.h"
#include "node_segment_convertions.h"
#include "params.h"

/* Propagation of an action potential along a cylinder using a 1-D
   Finite Element representation of the cable equation.
   Motion is initiated by injected current at x = EX_INPUT_NODE.
   Action potential code is based on Altenberger et al. 2001 code */

void runPharynxSimulator(const Params& params) {

	/*****************
     * Initializations
     *****************/

	/* Set basic values */

	double v_l	= Params::V_L  - Params::VM;		//Leak		reversal potential, relative to transmembrane equilibrium potential [mV]
	double v_k	= Params::V_K  - Params::VM;		//Potassium reversal potential, relative to transmembrane equilibrium potential [mV]
	double v_na = Params::V_NA - Params::VM;		//Sodium    reversal potential, relative to transmembrane equilibrium potential [mV]
	double g_l	= 0.057651;							//Leak current conductance [m Siemens/cm^2]
	double aNow = Params::EX_INPUT_AMP;				//Injected current at node x = EX_INPUT_NODE at (timeStep)	[mu amp]

	boolVecN1 print;
	for (int i = 0; i < Params::N + 1; i++) {
		print[i] = true;
	}

	/* Initialize geometry
	 *
	 * Parameters:
	 * x						Points along x-axis at the edges of the segments, defining segments position/length [cm]
	 * 							Segment i is bounded at its left (right) end by x[i] (x[i+1]). Have N+1 x's (x[0]..x[N])
	 * r						Pharyngeal outer radius at x [cm]
	 * p						Pharyngeal outer perimeter at x [cm]
	 * a						Pharyngeal cross-section area at x [cm^2]
	 * segX						Points along x-axis at the centers of the segments [cm]. Have N segX's (segX[0]..segX[N-1])
	 * segL						Segment's length [cm]. Have N segments (segL[0]..segL[N-1])
	 * segR						Muscle's length at rest	[cm]
	 * a_Ica_left				Sub-membrane shell's area through which trans-membrane calcium current flows at left  half of segment [cm^2]. Have N a_left's  (a_left[0]..a_left[N-1])
	 * a_Ica_right				Sub-membrane shell's area through which trans-membrane calcium current flows at right half of segment [cm^2]. Have N a_right's (a_right[0]..a_right[N-1])
	 * a_x_shell				Sub-membrane shell's cross-section area of x-axis calcium diffusion between adjacent x-axis segments [cm^2]. Have N-1 a_x's (a_x[0]..a_x[N-2])
	 * a_x_cell05				Inner cell shell's   cross-section area of x-axis calcium diffusion between adjacent x-axis segments [cm^2], at (timeStep+0.5). Have N-1 a_x's (a_x[0]..a_x[N-2])
	 * a_z						Shell's area of z-axis calcium diffusion between adjacent shells of a segment, i.e., sub-membrane- and inner-cell shells [cm^2]. Have N a_z's  (a_z[0]..a_z[N-1])
	 * v_shell					Sub-membrane shell's volume [cm^3]. Have N v's (v[0]..v[N-1])
	 * v_shell_{lef,righ}t		Volume of a half sub-membrane shell found to the {lef,righ}t of node i [m^3]. Have N v's (v[0]..v[N-1])
	 * 							Note: need these volumes to be in [m^3] in order to properly calculate calcium concentration, which is in mol/L
	 * v_cell					Inner cell shell's volume [cm^3]. Have N v's (v[0]..v[N-1])
	 */

	auto x_p				= make_unique<VectorN1>();	auto& x					= *x_p;					x = {};
	auto r_p				= make_unique<VectorN1>();	auto& r					= *r_p;					r = {};
	auto p_p				= make_unique<VectorN1>();	auto& p					= *p_p;					p = {};
	auto a_p				= make_unique<VectorN1>();	auto& a					= *a_p;					a = {};
	auto segX_p				= make_unique<VectorN>();	auto& segX				= *segX_p;				segX = {};
	auto x_all_p			= make_unique<Vector2N1>(); auto& x_all				= *x_all_p;				x_all = {};
	auto segL_p				= make_unique<VectorN>();	auto& segL				= *segL_p;				segL = {};
	auto segR_p				= make_unique<VectorN>();	auto& segR				= *segR_p;				segR = {};
	auto a_Ica_left_p	    = make_unique<VectorN>();	auto& a_Ica_left	  	= *a_Ica_left_p;		a_Ica_left = {};
	auto a_Ica_right_p	    = make_unique<VectorN>();	auto& a_Ica_right	  	= *a_Ica_right_p;		a_Ica_right = {};
	auto a_x_shell_left_p   = make_unique<VectorN>();	auto& a_x_shell_left  	= *a_x_shell_left_p;	a_x_shell_left = {};
	auto a_x_shell_right_p  = make_unique<VectorN>();	auto& a_x_shell_right 	= *a_x_shell_right_p;	a_x_shell_right = {};
	auto a_x_cell05_left_p  = make_unique<VectorN>();	auto& a_x_cell05_left 	= *a_x_cell05_left_p;	a_x_cell05_left = {};
	auto a_x_cell05_right_p = make_unique<VectorN>();	auto& a_x_cell05_right	= *a_x_cell05_right_p;	a_x_cell05_right = {};
	auto a_z_p			    = make_unique<VectorN>();	auto& a_z			  	= *a_z_p;				a_z = {};
	auto v_seg_init_p	    = make_unique<VectorN>();	auto& v_seg_init	  	= *v_seg_init_p;		v_seg_init = {};
	auto v_shell_p		    = make_unique<VectorN>();	auto& v_shell		  	= *v_shell_p;			v_shell = {};
	auto v_shell_left_p     = make_unique<VectorN>();	auto& v_shell_left	  	= *v_shell_left_p;		v_shell_left = {};
	auto v_shell_right_p    = make_unique<VectorN>();	auto& v_shell_right	  	= *v_shell_right_p;		v_shell_right = {};
	auto v_cell05_p		    = make_unique<VectorN>();	auto& v_cell05		  	= *v_cell05_p;			v_cell05 = {};

	std::cout << "Try opening file" << std::endl;
	BuildFiniteElementMesh::buildMesh(x, r);
	std::cout << "Opened file successfully" << std::endl;
	BuildFiniteElementMesh::setPerimeter(r, p);
	BuildFiniteElementMesh::setCrossSection(r, a);
	InitSegments::setSegX(x, segX);
	InitSegments::setXAll(x, segX, x_all);
	InitSegments::setSegLength(x, segL);
	InitSegments::setSegRadius(x, segX, r, segR);
	InitSegments::setAIca(segL, r, segR, a_Ica_left, a_Ica_right);
	InitSegments::setAxShellLeft(r, a_x_shell_left);
	InitSegments::setAxShellRight(r, a_x_shell_right);
	InitSegments::setAz(segL, r, a_z);
	InitSegments::setVShell(segL, r, v_shell);
	InitSegments::setVShellLeft(segL, r, segR, v_shell_left);
	InitSegments::setVShellRight(segL, r, segR, v_shell_right);
	InitSegments::setVSegmentInit(segL, r, v_seg_init);

	/* Initialize electrophysiology
	 *
	 * Parameters:
	 * v					Membrane potential at each x along the axon [mV], at (timeStep). Have N+1 v's (v[0]..v[N])
	 * v_prev				Membrane potential at each x along the axon [mV], at (timeStep-1)
	 * h_cca1				Opening probability of H&H voltage-dependent inactivation gate of cca-1  channels, at (timeStep+0.5)
	 * m_cca1				Opening probability of H&H voltage-dependent activation   gate of cca-1  channels, at (timeStep+0.5)
	 * f_egl19				Opening probability of H&H calcium-dependent inactivation gate of egl-19 channels, at (timeStep+0.5)
	 * m_egl19				Opening probability of H&H voltage-dependent activation   gate of egl-19 channels, at (timeStep+0.5)
	 * h_exp2				Opening probability of H&H voltage-dependent inactivation gate of exp-2  channels, at (timeStep+0.5)
	 * m_exp2				Opening probability of H&H voltage-dependent activation   gate of exp-2  channels, at (timeStep+0.5)
	 * Gca					Gca(V,CaIn(V,t),[Ca+2]out) [C/L]
	 * segCa_(sh,c)ell		[Ca+2]in at the sub-membrane- or inner-cell shell, at the middle of each segment, at the current  time step [mol/L]. Have N segCa's (ca[0]..ca[N-1])
	 * segCa_(sh,c)ell_prev	[Ca+2]in at the sub-membrane- or inner-cell shell, at the middle of each segment, at the previous time step [mol/L]
	 * ca_(sh,c)ell			[Ca+2]in at the sub-membrane- or inner-cell shell, at each x, at the current  time step [mol/L]. Have N+1 ca's (ca[0]..ca[N])
	 * ca_(sh,c)ell_prev	[Ca+2]in at the sub-membrane- or inner-cell shell, at each x, at the previous time step [mol/L]
	 * C05					[amp*cm/mol]
	 * E05					[unit-less]
	 */

	auto v_p				= make_unique<VectorN1>();	auto& v				   = *v_p;					v = {};
	auto v_prev_p			= make_unique<VectorN1>();	auto& v_prev		   = *v_prev_p;				v_prev = {};
	auto h_cca1_p			= make_unique<VectorN1>();	auto& h_cca1		   = *h_cca1_p;				h_cca1 = {};
	auto m_cca1_p			= make_unique<VectorN1>();	auto& m_cca1		   = *m_cca1_p;				m_cca1 = {};
	auto f_egl19_p			= make_unique<VectorN1>();	auto& f_egl19		   = *f_egl19_p;			f_egl19 = {};
	auto m_egl19_p			= make_unique<VectorN1>();	auto& m_egl19		   = *m_egl19_p;			m_egl19 = {};
	auto gca05_p			= make_unique<VectorN1>();	auto& gca05			   = *gca05_p;				gca05 = {};
	auto c05_p				= make_unique<VectorN1>();	auto& c05			   = *c05_p;				c05 = {};
	auto e05_p				= make_unique<VectorN1>();	auto& e05			   = *e05_p;				e05 = {};
	auto h_exp2_p			= make_unique<VectorN1>();	auto& h_exp2		   = *h_exp2_p;				h_exp2 = {};
	auto m_exp2_p			= make_unique<VectorN1>();	auto& m_exp2		   = *m_exp2_p;				m_exp2 = {};
	auto ca_shell_p			= make_unique<VectorN>();	auto& ca_shell		   = *ca_shell_p;			ca_shell = {};
	auto ca_shell_prev_p	= make_unique<VectorN>();	auto& ca_shell_prev	   = *ca_shell_prev_p;		ca_shell_prev = {};
	auto ca_cell_p			= make_unique<VectorN>();	auto& ca_cell		   = *ca_cell_p;			ca_cell = {};
	auto ca_cell_prev_p		= make_unique<VectorN>();	auto& ca_cell_prev	   = *ca_cell_prev_p;		ca_cell_prev = {};
	auto ca_shell05_nodes_p	= make_unique<VectorN1>();	auto& ca_shell05_nodes = *ca_shell05_nodes_p;	ca_shell05_nodes = {};
	auto ca_shell_nodes_p	= make_unique<VectorN1>();	auto& ca_shell_nodes   = *ca_shell_nodes_p;		ca_shell_nodes = {};
	auto tau_ca_sec_p		= make_unique<VectorN>();	auto& tau_ca_sec	   = *tau_ca_sec_p;			tau_ca_sec = {};

	for (int seg = 0; seg < Params::N; seg++) {
		ca_shell[seg]	   = Params::CA_IN_REST;
		ca_shell_prev[seg] = Params::CA_IN_REST;
		ca_cell[seg]	   = Params::CA_IN_REST;
		ca_cell_prev[seg]  = Params::CA_IN_REST;
		tau_ca_sec[seg]	= 1E-3 * params.get(GeneType::ca_in_tau);
	}
	for (int seg = Params::POST_ISTH_START_NODE; seg < Params::PM6_START_NODE; seg++) {
		tau_ca_sec[seg]	= 1E-3 * 1.0;
	}

	for (int i = 0; i < Params::N + 1; i++) {
		ca_shell05_nodes[i] = Params::CA_IN_REST;
		ca_shell_nodes[i] = Params::CA_IN_REST;
		v[i] = 1.4;
		v_prev[i] = 1.4;
		AlphaBeta h_cca1_ab = HHFunctions::alphaBeta_h_cca1(params, v[i]);
		h_cca1[i] = 1.0 / (1.0 + h_cca1_ab.b_prob / h_cca1_ab.a_prob);
		AlphaBeta m_cca1_ab = HHFunctions::alphaBeta_m_cca1(params, v[i]);
		m_cca1[i] = 1.0 / (1.0 + m_cca1_ab.b_prob / m_cca1_ab.a_prob);
		AlphaBeta f_egl19_ab = HHFunctions::alphaBeta_f_egl19(params, v[i], ca_shell_nodes[i]);
		f_egl19[i] = 1.0 / (1.0 + f_egl19_ab.b_prob / f_egl19_ab.a_prob);
		AlphaBeta m_egl19_ab = HHFunctions::alphaBeta_m_egl19(params, v[i]);
		m_egl19[i] = 1.0 / (1.0 + m_egl19_ab.b_prob / m_egl19_ab.a_prob);
		gca05[i] = HHFunctions::calcGca(v_prev[i], v[i], ca_shell05_nodes[i]);
		AlphaBeta h_exp2_ab = HHFunctions::alphaBeta_h_exp2(params, v[i]);
		h_exp2[i] = 1.0 / (1.0 + h_exp2_ab.b_prob / h_exp2_ab.a_prob);
		AlphaBeta m_exp2_ab = HHFunctions::alphaBeta_m_exp2(params, v[i]);
		m_exp2[i] = 1.0 / (1.0 + m_exp2_ab.b_prob / m_exp2_ab.a_prob);
	}

	/* Initialize muscles
	 *
	 * Parameters:
	 * M, Mp, AM, AMp	Fraction of bounded (AM, AMp) and unbounded (M, Mp) myosin heads to actin, relative to the total number of myosin heads,
	 * 					L_ACTIN_BINDING_SITES/H, at the current time step [unit-less]
	 * MLCP_P			Activity-level of MLCP (fraction of activated MLCP) [0..1]
	 * K1, K6			Rate of myosin (M), actomyosin (AM) phosphorylation by MLCK [1/s]
	 * RF				Concentration of muscle's relaxation-inducing factor [mol/L]
	 * l_muscle_prev	Muscle's length [cm] of a (WHOLE) sarcomere, at each x, at the previous time step
	 * l_muscle			Muscle's length [cm] of a (WHOLE) sarcomere, at each x, at the current  time step.  Have N+1 l_muscle's (l_muscle[0]..l_muscle[N])
	 * v_muscle_prev	Muscle's speed of shortening [cm/s] (> 0), at the previous time step
	 * v_muscle			Muscle's speed of shortening [cm/s] (> 0), at the current  time step
	 * u_sarc			Normalized speed of shortening [1/s] (> 0) of a HALF sarcomere, at the current time step
	 * z05				Lumen's radius [cm], at the middle of the segment, at (timeStep+0.5). Have N   z05's (z05[0]..z05[N-1])
	 * z				Lumen's radius [cm], at the middle of the segment, at (timeStep).	  Have N   z's (z[0]..z[N-1])
	 * z05_nodes		Lumen's radius [cm], at the edge   of the segment, at (timeStep+0.5). Have N+1 z05_nodes (z05_nodes[0]..z05_nodes[N])
	 * z_nodes			Lumen's radius [cm], at the edge   of the segment, at (timeStep).	  Have N+1 z_nodes (z_nodes[0]..z_nodes[N])
	 * m				1/3 of total mass of muscles found in each segment [g] (1/3 since each segment is composed of 3 muscle cells, and for each muscle cell calculate force of the WHOLE muscle <=> sarcomere)
	 * F_CE				Muscle's contractile element's force, of a HALF sarcomere, at (timeStep) [gr·cm/s^2]
	 * F_PE_K			Pharyngeal passive element's spring force, of the (WHOLE) sarcomere, at (timeStep) [gr·cm/s^2]
	 * F_PE_D			Pharyngeal passive element's damper force, of the (WHOLE) sarcomere, at (timeStep) [gr·cm/s^2]
	 * PE_K				Passive element's spring constant [g/s^2]
	 * PE_D				Passive element's damper constant [g/s]
	 * F_OE				Connective tissue's force that acts on a (WHOLE) sarcomere, applied by neighboring segments at (timeStep) [gr·cm/s^2]
	 * Q, Q_estim		Distributions' moments during the current and next time steps, respectively [unit-less]
	 * n_dist, n0l_dist	dist.mean, dist.standard_deviation [unit-less]
	 */

	auto M_p				 = make_unique<VectorN>();		auto& M					= *M_p;					M = {};
	auto Mp_p				 = make_unique<VectorN>();		auto& Mp				= *Mp_p;				Mp = {};
	auto AM_p				 = make_unique<VectorN>();		auto& AM				= *AM_p;				AM = {};
	auto AMp_p			 	 = make_unique<VectorN>();		auto& AMp				= *AMp_p;				AMp = {};
	auto CF_p			 	 = make_unique<VectorN>();		auto& CF				= *CF_p;				CF = {};
	auto MLCP_P_p			 = make_unique<VectorN>();		auto& MLCP_P			= *MLCP_P_p;			MLCP_P = {};
	auto K2vec_p			 = make_unique<VectorN>();		auto& K2vec				= *K2vec_p;				K2vec = {};
	auto K1_p				 = make_unique<VectorN>();		auto& K1				= *K1_p;				K1 = {};
	auto l_muscle_p			 = make_unique<VectorN>();		auto& l_muscle			= *l_muscle_p;			l_muscle = {};
	auto l_muscle_prev_p	 = make_unique<VectorN>();		auto& l_muscle_prev		= *l_muscle_prev_p;		l_muscle_prev = {};
	auto l_muscle05_p		 = make_unique<VectorN>();		auto& l_muscle05		= *l_muscle05_p;		l_muscle05 = {};
	auto v_muscle_p			 = make_unique<VectorN>();		auto& v_muscle			= *v_muscle_p;			v_muscle = {};
	auto v_muscle_prev_p	 = make_unique<VectorN>();		auto& v_muscle_prev		= *v_muscle_prev_p; 	v_muscle_prev = {};
	auto u_sarc_p			 = make_unique<VectorN>();		auto& u_sarc			= *u_sarc_p;			u_sarc = {};
	auto z05_p				 = make_unique<VectorN>();		auto& z05				= *z05_p;				z05 = {};
	auto z_p				 = make_unique<VectorN>();		auto& z					= *z_p;					z = {};
	auto m_p				 = make_unique<VectorN>();		auto& m					= *m_p;					m = {};
	auto CE_K_p				 = make_unique<VectorN>();		auto& CE_K				= *CE_K_p;				CE_K = {};
	auto PE_K_p				 = make_unique<VectorN>();		auto& PE_K				= *PE_K_p;				PE_K = {};
	auto PE_D_p				 = make_unique<VectorN>();		auto& PE_D				= *PE_D_p;				PE_D = {};
	auto F_CE_p				 = make_unique<VectorN>();		auto& F_CE				= *F_CE_p;				F_CE = {};
	auto F_CE_prev_p		 = make_unique<VectorN>();		auto& F_CE_prev			= *F_CE_prev_p;			F_CE_prev = {};
	auto F_CE05_p			 = make_unique<VectorN>();		auto& F_CE05			= *F_CE05_p;			F_CE05 = {};
	auto F_OE_p				 = make_unique<VectorN>();		auto& F_OE				= *F_OE_p;				F_OE = {};
	auto F_OE_prev_p		 = make_unique<VectorN>();		auto& F_OE_prev			= *F_OE_prev_p;			F_OE_prev = {};
	auto F_OE05_p			 = make_unique<VectorN>();		auto& F_OE05			= *F_OE05_p;			F_OE05 = {};
	auto F_PE_K_p			 = make_unique<VectorN>();		auto& F_PE_K			= *F_PE_K_p;			F_PE_K = {};
	auto F_PE_D_p			 = make_unique<VectorN>();		auto& F_PE_D			= *F_PE_D_p;			F_PE_D = {};
	auto muscle_activation_p = make_unique<VectorN>();		auto& muscle_activation = *muscle_activation_p;	muscle_activation = {};
	auto z05_nodes_p		 = make_unique<VectorN1>();		auto& z05_nodes			= *z05_nodes_p;			z05_nodes = {};
	auto z_nodes_p			 = make_unique<VectorN1>();		auto& z_nodes			= *z_nodes_p;			z_nodes = {};
	auto z_all_p			 = make_unique<Vector2N1>();	auto& z_all				= *z_all_p;				z_all = {};
	auto Q_p				 = make_unique<MomentsN>(); 	auto& Q					= *Q_p;					Q = {};
	auto Q_estim_p			 = make_unique<MomentsN>();		auto& Q_estim			= *Q_estim_p;			Q_estim = {};
	auto n_dist_p			 = make_unique<DistN>();		auto& n_dist			= *n_dist_p;			n_dist = {};
	auto n0l_dist_p			 = make_unique<DistN>();		auto& n0l_dist			= *n0l_dist_p;			n0l_dist = {};
	auto CE_num_AM_p		 = make_unique<VectorN>();		auto& CE_num_AM			= *CE_num_AM_p;			CE_num_AM = {};

	double minimal_volume = v_seg_init[0];
	MuscleFunctions::calc_m_PED_CEK_PEK(v_seg_init, m, PE_D, CE_K, PE_K, &print);
	for (int i = 0; i < Params::N; i++) {
		K1[i] = MuscleFunctions::calcK1(i, ca_cell[i], &print);
		CF[i] = 0.0;
		MLCP_P[i] = 1.0;
		K2vec[i] = Params::K2_MAX * MLCP_P[i];
		l_muscle_prev[i] = segR[i];
		l_muscle[i]   = l_muscle_prev[i];
		l_muscle05[i] = l_muscle_prev[i];
		v_muscle_prev[i] = 0.0;
		v_muscle[i] = 0.0;
		u_sarc[i] = 0.5 * v_muscle[i] / Params::H;
		z05[i] = 0.0;
		z[i] = 0.0;
		F_OE[i] = 0.0;
		CE_num_AM[i] = v_seg_init[i] / minimal_volume * Params::CE_NUM_AM;
		AMp[i] = 0.0;
		AM[i]  = 0.0;
		Mp[i]  = 0.0;
		muscle_activation[i] = MuscleFunctions::calcMuscleActivation(segR[i], l_muscle[i]);
		Q[i][0] = AM[i] + AMp[i];
		M[i] = muscle_activation[i] * Params::L_MYOSIN_BINDING_SITES / Params::H - Q[i][0] - Mp[i];
		if (Q[i][0] != 0){
			double x_min = 0.0;
			double x_max = Params::L_MYOSIN_BINDING_SITES / Params::H;
			double n_mean = 0.5 * (x_min + x_max);									// assume n is a uniform distribution between [0, l/h] at rest
			double n_standard_deviation = std::sqrt(1.0 / 12.0) * (x_max - x_min);
			Q[i][1] = n_mean * Q[i][0];
			Q[i][2] = Q[i][0] * (square(n_standard_deviation) + square(Q[i][1] / Q[i][0]));
		} else {
			Q[i][1] = 0.0;
			Q[i][2] = 0.0;
		}
		F_CE[i] = muscle_activation[i] * CE_num_AM[i] * CE_K[i] * square(Params::H) / (2 * Params::L_MYOSIN_BINDING_SITES) * Q[i][1];
	}
	for (int i = 0; i < Params::N + 1; i++) {
		z05_nodes[i] = 0.0;
		z_nodes[i] = 0.0;
	}

	/* Initialize particle motion
	 *
	 * Parameters:
	 * inner_r			Radius of the lumen, at each x, at the current time step [cm]. Have N+1 inner_r's (inner_r[0]..inner_r[N])
	 * inner_r_max		Maximal inner radius [cm]. For all segments, assume that the lumen opens to a triangle at 2z_max
	 * a_lumen			Cross-sectional area of the lumen  [cm^2], at (timeStep)
	 * v_lumen			Volume of the lumen of a segment   [cm^2], at (timeStep)
	 * v_lumen_prev		Volume of the lumen of a segment   [cm^2], at (timeStep-1)
	 * flow				Flow through the anterior boundary [cm^3/s], at (timeStep)
	 * v_center			Velocity at the center of a triradiate-shaped lumen [cm/s], at (timeStep)
	 */

	auto z_max_p		= make_unique<VectorN>();	auto& z_max		   = *z_max_p;			z_max = {};
	auto z_nodes_max_p	= make_unique<VectorN1>();	auto& z_nodes_max  = *z_nodes_max_p;	z_nodes_max = {};
	auto a_lumen_p		= make_unique<Vector2N1>(); auto& a_lumen	   = *a_lumen_p;		a_lumen = {};
	auto v_lumen_p		= make_unique<Vector2N>(); 	auto& v_lumen	   = *v_lumen_p;		v_lumen = {};
	auto v_lumen_prev_p	= make_unique<Vector2N>();	auto& v_lumen_prev = *v_lumen_prev_p;	v_lumen_prev = {};
	auto flow_p			= make_unique<Vector2N>();	auto& flow		   = *flow_p;			flow = {};
	auto vel_mean_p		= make_unique<Vector2N1>();	auto& vel_mean	   = *vel_mean_p;		vel_mean = {};
	auto vel_center_p	= make_unique<Vector2N1>();	auto& vel_max	   = *vel_center_p;		vel_max = {};
	Particle part1, part2, part3, part4, part5;

	FlowFunctions::loadZMax(z_max);
	NodeSegmentConvertions::calcZNodesMax(x, segX, z_max, z_nodes_max);
	NodeSegmentConvertions::calcAxCell05Left (r, z_nodes_max, z05_nodes, a_x_cell05_left);
	NodeSegmentConvertions::calcAxCell05Right(r, z_nodes_max, z05_nodes, a_x_cell05_right);
	NodeSegmentConvertions::calcVCell05(segL, r, z_nodes_max, z05_nodes, v_cell05);
	FlowFunctions::calcLumenArea(z_nodes_max, z_nodes, z_max, z, a_lumen, &print);
	FlowFunctions::calcLumenVolume(x, segL, z_nodes_max, z_nodes, z_max, z, v_lumen, &print);
	part1.diameter = 0.52 * 1E-4;					// Mean diameter of E. coli, OP50 strain [cm] (Fang-Yen et al., 2009, Table S3)
	part1.seg = 0;
	part1.position = x[part1.seg];					// <=> at 0 um; Particle position at the current time step [cm]
	part2.diameter = 0.52 * 1E-4;
	part2.seg = 0;
	part2.position = x[part2.seg];					// <=> at 0 um; Particle position at the current time step [cm]
	part3.diameter = 0.52 * 1E-4;
	part3.seg = 2 * 71;
	part3.position = x[int(0.5 * part3.seg)];		// <=> at 3.3 um; Particle position at the current time step [cm]
	part4.diameter = 0.52 * 1E-4;
	part4.seg = 2 * 529;
	part4.position = x[int(0.5 * part4.seg)];		// <=> at 24.8 um; Particle position at the current time step [cm]
	part5.diameter = 0.52 * 1E-4;
	part5.seg = 2 * 2205;
	part5.position = x[int(0.5 * part5.seg)];		// <=> at 103.2 um; Particle position at the current time step [cm]

	/* Build the expanded form of the modeled functions */

	auto xi_p			= make_unique<TriMatN1>();	auto& xi		   = *xi_p;				xi = {};				//[cm]
	auto chi_p			= make_unique<TriMatN1>();	auto& chi		   = *chi_p;			chi = {};				//[m Siemens·cm]
	auto phiV_p			= make_unique<TriMatN1>();	auto& phiV		   = *phiV_p;			phiV = {};				//[cm^2]
	auto etaV_p			= make_unique<TriMatN1>();	auto& etaV		   = *etaV_p;			etaV = {};				//[m Siemens] --> [mu F]
	auto psiV_p			= make_unique<TriMatN1>();	auto& psiV		   = *psiV_p;			psiV = {};				//[mu F]
	auto phiCa_p		= make_unique<TriMatN>();	auto& phiCa		   = *phiCa_p;			phiCa = {};				//[cm]
	auto chiCa_left_p	= make_unique<TriMatN>();	auto& chiCa_left   = *chiCa_left_p;		chiCa_left = {};		//[unit-less]
	auto chiCa_right_p	= make_unique<TriMatN>();	auto& chiCa_right  = *chiCa_right_p;	chiCa_right = {};		//[unit-less]
	auto rhoCa_p		= make_unique<VectorN>();	auto& rhoCa		   = *rhoCa_p;			rhoCa = {};				//[cm]
	auto ksiCa_shell_p	= make_unique<TriMatN>();	auto& ksiCa_shell  = *ksiCa_shell_p;	ksiCa_shell = {};		//[1/cm]
	auto ksiCa_cell05_p = make_unique<TriMatN>();	auto& ksiCa_cell05 = *ksiCa_cell05_p;	ksiCa_cell05 = {};		//[1/cm]
	auto muCa_shell_p	= make_unique<TriMatN>();	auto& muCa_shell   = *muCa_shell_p;		muCa_shell = {};		//[1/cm]
	auto muCa_cell_p	= make_unique<TriMatN>();	auto& muCa_cell	   = *muCa_cell_p;		muCa_cell = {};			//[1/cm]
	auto nuCa_left_p	= make_unique<TriMatN>();	auto& nuCa_left	   = *nuCa_left_p;		nuCa_left = {};			//[cm/s]
	auto nuCa_right_p	= make_unique<TriMatN>();	auto& nuCa_right   = *nuCa_right_p;		nuCa_right = {};		//[cm/s]
	auto etaCa_shell_p	= make_unique<TriMatN>();	auto& etaCa_shell  = *etaCa_shell_p;	etaCa_shell = {};		//[cm]
	auto etaCa_cell_p	= make_unique<TriMatN>();	auto& etaCa_cell   = *etaCa_cell_p;		etaCa_cell = {};		//[cm]
	auto psiCa_shell_p	= make_unique<TriMatN>();	auto& psiCa_shell  = *psiCa_shell_p;	psiCa_shell = {};		//[cm]
	auto psiCa_cell_p	= make_unique<TriMatN>();	auto& psiCa_cell   = *psiCa_cell_p;		psiCa_cell = {};		//[cm]
	auto chiL_p			= make_unique<TriMatN>();	auto& chiL		   = *chiL_p;			chiL = {};				//[gr·cm]
	auto phiL_p			= make_unique<TriMatN>();	auto& phiL		   = *phiL_p;			phiL = {};				//[cm]
	auto ksiL_p			= make_unique<TriMatN>();	auto& ksiL		   = *ksiL_p;			ksiL = {};				//[gr·cm/s^2]
	auto rhoL_p			= make_unique<VectorN>();	auto& rhoL		   = *rhoL_p;			rhoL = {};				//[gr·cm^2/s^2]
	auto etaL_p			= make_unique<TriMatN>();	auto& etaL		   = *etaL_p;			etaL = {};				//[gr·cm/s]
	auto psiL_p			= make_unique<TriMatN>();	auto& psiL		   = *psiL_p;			psiL = {};				//[gr·cm/s]

	BuildAxialCurrentFunction::buildXi(x, xi);
	BuildAxialCurrentFunction::buildChi(a, chi);
	BuildCableFunctionSpaceIntegration::LHS(x, p, phiV);
	BuildCableFunctionSpaceIntegration::RHS(x, a, etaV);
	BuildCalciumFunctionSpaceIntegration::LHS(segX, phiCa);
	BuildCalciumFunctionSpaceIntegration::RHS1(segX, a_Ica_left,  v_shell, chiCa_left);
	BuildCalciumFunctionSpaceIntegration::RHS1(segX, a_Ica_right, v_shell, chiCa_right);
	BuildCalciumFunctionSpaceIntegration::RHS2(segX, rhoCa);
	BuildCalciumFunctionSpaceIntegration::RHS3_left(nuCa_left);
	BuildCalciumFunctionSpaceIntegration::RHS3_right(nuCa_right);
	BuildCalciumFunctionSpaceIntegration::RHS4_left_minus_right(a_x_shell_left, a_x_shell_right, v_shell, ksiCa_shell);
	BuildCalciumFunctionSpaceIntegration::RHS5(segX, a_z, v_shell, muCa_shell);
	BuildCalciumFunctionSpaceIntegration::RHS6(segX, a_z, muCa_cell);
	BuildLumenFunctionSpaceIntegration::LHS(segX, m, chiL);
	BuildLumenFunctionSpaceIntegration::RHS1(segX, phiL);
	BuildLumenFunctionSpaceIntegration::LHS(segX, PE_K, ksiL);
	BuildLumenFunctionSpaceIntegration::RHS2(segX, PE_K, segR, rhoL);
	BuildLumenFunctionSpaceIntegration::LHS(segX, PE_D, etaL);
	BuildCableFunctionTimeIntegration::integrate1(phiV, etaV, psiV);
	BuildCalciumFunctionTimeIntegration::integrate1(params, tau_ca_sec, phiCa, ksiCa_shell, etaCa_shell, psiCa_shell, 1);
	BuildCalciumFunctionTimeIntegration::integrate1(params, tau_ca_sec, phiCa, ksiCa_shell, etaCa_cell,  psiCa_cell,  0);
	BuildLumenFunctionTimeIntegration::integrate1(ksiL, etaL, psiL);

	/* Initialize time-dependent variables
	 *
	 * Parameters:
	 * g																				[m Siemens/cm^2]
	 * d																				[mu amp/cm^2]
	 * amp_: l, cca1_ca, cca1_na, egl19_ca, egl19_na, exp2, c, m 						[mu amp/cm^2]
	 * amp_ca			Trans-membrane calcium current at (timeStep + 0.5)				[mu amp/cm^2]
	 * b
	 * BV
	 * BCaShell, BCaCell
	 * BBuffShell, BBuffCell
	 * BL																				[gr·cm^2/s]
	 * alV, arV
	 * alCaShell, arCaShell, alCaCell, arCaCell
	 * alBuffShell, arBuffShell, alBuffCell, arBuffCell
	 * alL, arL																			[gr·cm/s]
	 * K2vec, K5		Rate of myosin (Mp), actomyosin (AMp) dephosphorylation by MLCP [1/s]
	 * betaVec																			[1/s]
	 */

	auto g_p				 = make_unique<VectorN1>();	auto& g					= *g_p;					g = {};
	auto d_p				 = make_unique<VectorN1>(); auto& d					= *d_p;					d = {};
	auto amp_l_p			 = make_unique<VectorN1>(); auto& amp_l				= *amp_l_p;				amp_l = {};
	auto amp_cca1_ca_p		 = make_unique<VectorN1>(); auto& amp_cca1_ca		= *amp_cca1_ca_p;		amp_cca1_ca = {};
	auto amp_cca1_na_p		 = make_unique<VectorN1>(); auto& amp_cca1_na		= *amp_cca1_na_p;		amp_cca1_na = {};
	auto amp_egl19_ca_p		 = make_unique<VectorN1>(); auto& amp_egl19_ca		= *amp_egl19_ca_p;		amp_egl19_ca = {};
	auto amp_egl19_na_p		 = make_unique<VectorN1>(); auto& amp_egl19_na		= *amp_egl19_na_p;		amp_egl19_na = {};
	auto amp_ca05_p			 = make_unique<VectorN1>(); auto& amp_ca05			= *amp_ca05_p;			amp_ca05 = {};
	auto amp_exp2_p			 = make_unique<VectorN1>(); auto& amp_exp2			= *amp_exp2_p;			amp_exp2 = {};
	auto amp_c_p			 = make_unique<VectorN1>(); auto& amp_c				= *amp_c_p;				amp_c = {};
	auto amp_m_p			 = make_unique<VectorN1>(); auto& amp_m				= *amp_m_p;				amp_m = {};
	auto b_p				 = make_unique<VectorN1>(); auto& b					= *b_p;					b = {};
	auto BV_p				 = make_unique<VectorN1>(); auto& BV				= *BV_p;				BV = {};
	auto BCa_shell_p		 = make_unique<VectorN>();  auto& BCa_shell			= *BCa_shell_p;			BCa_shell = {};
	auto BCa_cell_p			 = make_unique<VectorN>();  auto& BCa_cell			= *BCa_cell_p;			BCa_cell = {};
	auto BL_p				 = make_unique<VectorN>();  auto& BL				= *BL_p;				BL = {};
	auto alV_p				 = make_unique<TriMatN1>(); auto& alV				= *alV_p;				alV = {};
	auto arV_p				 = make_unique<TriMatN1>(); auto& arV				= *arV_p;				arV = {};
	auto alCaShell_p		 = make_unique<TriMatN>();  auto& alCaShell			= *alCaShell_p;			alCaShell = {};
	auto arCaShell_p		 = make_unique<TriMatN>();  auto& arCaShell			= *arCaShell_p;			arCaShell = {};
	auto alCaCell_p			 = make_unique<TriMatN>();  auto& alCaCell			= *alCaCell_p;			alCaCell = {};
	auto arCaCell_p			 = make_unique<TriMatN>();  auto& arCaCell			= *arCaCell_p;			arCaCell = {};
	auto alL_p				 = make_unique<TriMatN>();  auto& alL				= *alL_p;				alL = {};
	auto arL_p				 = make_unique<TriMatN>();  auto& arL				= *arL_p;				arL = {};
	auto betaVec_p			 = make_unique<VectorN>();  auto& betaVec			= *betaVec_p;			betaVec = {};
	auto muscle_phi10vec_p	 = make_unique<VectorN>();	auto& muscle_phi10vec   = *muscle_phi10vec_p;   muscle_phi10vec = {};
	auto muscle_phi2s01vec_p = make_unique<VectorN>();	auto& muscle_phi2s01vec = *muscle_phi2s01vec_p; muscle_phi2s01vec = {};
	auto muscle_phi2s02vec_p = make_unique<VectorN>();	auto& muscle_phi2s02vec = *muscle_phi2s02vec_p; muscle_phi2s02vec = {};

	std::ofstream outputFile;
	std::ostringstream oss;
	oss << "Output.csv";
	std::string out_file = oss.str();
	outputFile.open(out_file);
	//outputFile << "t, vm, ca_shell, ca_cell, l, z, v, u, a, F_PE_K, F_PE_D, 2*F_CE, F_neigh, Q0, Q1, Q2, M, Mp, AM, AMp, K1, K2, beta, phi10 ,phi20_0 ,phi20_1, muscle_act, CF, MLCP, n_dist_mean, n_dist_SD, n0l_dist_mean, n0l_dist_SD" << endl;//, KRF_generate, KRF_disassemble, K2_on, K2_off, K1_on, K1_off" << endl;
	//outputFile << "t, vm, cca1_ca, cca1_na, egl19_ca, egl19_na, exp2, ca_shell, ca_cell, ca_shell05_nodes, m0_cca1, h0_cca1, m0_egl19, m0_egl19^2, f0_egl19, m0_exp2^4, h0_exp2^4, gca05, segR-2l_muscle, 2l_muscle, Eca" << endl;//, KRF_generate, KRF_disassemble, K2_on, K2_off, K1_on, K1_off" << endl;
	outputFile << "t, part1_pos, part2_pos, part3_pos, part4_pos, part5_pos" << endl;

	/*for (int seg = 0; seg < Params::N; seg++) {
		outputFile << seg << ",";
	}
	outputFile << endl;
	for (int seg = 0; seg < Params::N; seg++) {
		outputFile << std::setprecision(5) << segX[seg] * 1E4 << ",";
	}
	outputFile << endl;
	for (int seg = 0; seg < Params::N; seg++) {
		outputFile << std::setprecision(5) << segR[seg] * 1E4 << ",";
	}
	outputFile << endl;*/

	double AP_t0 = 0.0;
	for (int timeStep = 1; timeStep <= Params::SIMULATED_TIME_POINTS + 1; timeStep++) {															//add 1 timeStep to simulation in order to have the 1st timeStep at resting potential
		double t = timeStep * Params::DT;
		if (timeStep % Params::EX_INPUT_INTERVAL == 0){AP_t0 = t; std::cout << AP_t0 << std::endl;}
		double h0_cca1, m0_cca1, f0_egl19, m0_egl19, Gca0, h0_exp2, m0_exp2;
		for (int i = 0; i < Params::N + 1; i++) {
			g[i] = g_l;																															//g[i][timeStep+0.5]
			d[i] = g_l * v_l;																													//d[i][timeStep+0.5]
			m0_cca1 = m_cca1[i];																												//m[i][timeStep-0.5]
			if (v[i] >= v_prev[i]) {
				m_cca1[i] = HHFunctions::calcGateOpeningProb(params, GatePart0::m, GatePart1::cca1, v[i], ca_shell_nodes[i], m_cca1[i]);		//m[i][timeStep+0.5]
			} else {
				m_cca1[i] = HHFunctions::calcGateOpeningProbHP(params, GatePart0::m, GatePart1::cca1, v[i], m_cca1[i]);
			}
			m0_cca1 = 0.5 * (m_cca1[i] + m0_cca1);																								//m[i][timeStep]
			h0_cca1 = h_cca1[i];																												//h[i][timeStep-0.5]	(at timeStep=1, h[i] contains h[i][0])
			h_cca1[i] = HHFunctions::calcGateOpeningProb(params, GatePart0::h, GatePart1::cca1, v[i], ca_shell_nodes[i], h_cca1[i]);			//h[i][timeStep+0.5]
			h0_cca1 = 0.5 * (h_cca1[i] + h0_cca1);																								//h[i][timeStep]
			m0_egl19 = m_egl19[i];																												//m[i][timeStep-0.5]
			m_egl19[i] = HHFunctions::calcGateOpeningProb(params, GatePart0::m, GatePart1::egl19, v[i], ca_shell_nodes[i], m_egl19[i]);			//m[i][timeStep+0.5]
			m0_egl19 = 0.5 * (m_egl19[i] + m0_egl19);																							//m[i][timeStep]
			f0_egl19 = f_egl19[i];																												//f[i][timeStep-0.5]
			f_egl19[i] = HHFunctions::calcGateOpeningProb(params, GatePart0::f, GatePart1::egl19, v[i], ca_shell_nodes[i], f_egl19[i]);			//f[i][timeStep+0.5];	use v[i][timeStep] and ca_shell_nodes[i][timeStep]
			f0_egl19 = 0.5 * (f_egl19[i] + f0_egl19);																							//f[i][timeStep]
			Gca0 = gca05[i];																													//Gca[i][timeStep-0.5]
			gca05[i] = HHFunctions::calcGca(v_prev[i], v[i], ca_shell05_nodes[i]);																//Gca[i][timeStep+0.5]
			Gca0 = 0.5 * (gca05[i] + Gca0);																										//Gca[i][timeStep]
			double hm205_cca1  = h_cca1[i]  * square(m_cca1[i]);
			double fm205_egl19 = f_egl19[i] * square(m_egl19[i]);
			c05[i] = HHFunctions::calcC05(params, v_prev[i], v[i], hm205_cca1, fm205_egl19);													//C05[i][timeStep+0.5]
			e05[i] = HHFunctions::calcE05(v_prev[i], v[i]);																						//E05[i][timeStep+0.5]
			g[i] += params.get(GeneType::g_cca1_na)	 * h_cca1[i]  * square(m_cca1[i]);
			d[i] += params.get(GeneType::g_cca1_na)	 * h_cca1[i]  * square(m_cca1[i])  * v_na;
			d[i] -= params.get(GeneType::p_cca1_ca)	 * h_cca1[i]  * square(m_cca1[i])  * gca05[i];
			g[i] += params.get(GeneType::g_egl19_na) * f_egl19[i] * square(m_egl19[i]);
			d[i] += params.get(GeneType::g_egl19_na) * f_egl19[i] * square(m_egl19[i]) * v_na;
			d[i] -= params.get(GeneType::p_egl19_ca) * f_egl19[i] * square(m_egl19[i]) * gca05[i];
			h0_exp2 = h_exp2[i];																												//h[i][timeStep-0.5]
			if (v[i] >= v_prev[i]) {
				h_exp2[i] = HHFunctions::calcGateOpeningProb  (params, GatePart0::h, GatePart1::exp2, v[i], ca_shell_nodes[i], h_exp2[i]);		//h[i][timeStep+0.5]
			} else {
				h_exp2[i] = HHFunctions::calcGateOpeningProbHP(params, GatePart0::h, GatePart1::exp2, v[i], h_exp2[i]);
			}
			h0_exp2 = 0.5 * (h_exp2[i] + h0_exp2);																								//h[i][timeStep]
			m0_exp2 = m_exp2[i];																												//m[i][timeStep-0.5]
			if (v[i] >= v_prev[i]) {
				m_exp2[i] = HHFunctions::calcGateOpeningProb  (params, GatePart0::m, GatePart1::exp2, v[i], ca_shell_nodes[i], m_exp2[i]);		//m[i][timeStep+0.5]
			} else {
				m_exp2[i] = HHFunctions::calcGateOpeningProbHP(params, GatePart0::m, GatePart1::exp2, v[i], m_exp2[i]);
			}
			m0_exp2 = 0.5 * (m_exp2[i] + m0_exp2);																								//m[i][timeStep]
			g[i] += params.get(GeneType::g_exp2) * quad(h_exp2[i]) * quad(m_exp2[i]);															//g[i][timeStep+0.5]
			d[i] += params.get(GeneType::g_exp2) * quad(h_exp2[i]) * quad(m_exp2[i]) * v_k;														//d[i][timeStep+0.5]
			amp_l[i]		= g_l * (v[i] - v_l);																								//amp_l[i][timeStep]
			amp_ca05[i]		= params.get(GeneType::p_cca1_ca)  * h_cca1[i]  * square(m_cca1[i])  * gca05[i] +
							  params.get(GeneType::p_egl19_ca) * f_egl19[i] * square(m_egl19[i]) * gca05[i];									//amp_ca[i][timeStep+0.5]
			amp_cca1_ca[i]	= params.get(GeneType::p_cca1_ca)  * h0_cca1  * square(m0_cca1)  * Gca0;											//amp_cca1_ca[i][timeStep]
			amp_cca1_na[i] 	= params.get(GeneType::g_cca1_na)  * h0_cca1  * square(m0_cca1)  * (v[i] - v_na);									//amp_cca1_na[i][timeStep]
			amp_egl19_ca[i] = params.get(GeneType::p_egl19_ca) * f0_egl19 * square(m0_egl19) * Gca0;											//amp_egl19_ca[i][timeStep]
			amp_egl19_na[i] = params.get(GeneType::g_egl19_na) * f0_egl19 * square(m0_egl19) * (v[i] - v_na);
			amp_exp2[i]		= params.get(GeneType::g_exp2) * quad(h0_exp2) * quad(m0_exp2) * (v[i] - v_k);										//amp_exp2[i][timeStep]
			amp_m[i]		= amp_l[i] + amp_cca1_ca[i] + amp_cca1_na[i] + amp_egl19_ca[i] + amp_egl19_na[i] + amp_exp2[i];						//amp_m[i][timeStep]
			// No Na+:
			//amp_m[i]		= amp_l[i] + amp_cca1_ca[i] + amp_egl19_ca[i] + amp_exp2[i];
			// EGL-19 alone:
			//amp_m[i]		= amp_l[i] + amp_egl19_ca[i] + amp_egl19_na[i];
		}

		/* Calculate axial current */

		TriMatOperations::triMatVecMultN1(chi, v, b);			//b[i][timeStep]
		TriMatOperations::triMatSolve2N1(xi, b, amp_c);			//amp_c[i][timeStep]	[mu A]
		/* Invert sign of amp_c left of EX_INPUT_NODE, where axial current flows in opposite direction */
		for (int i = 0; i < Params::EX_INPUT_NODE; i++) {
			amp_c[i] = -amp_c[i];
		}
		for (int i = 0; i < Params::N + 1; i++) {
			amp_c[i] = amp_c[i] / (p[i] * segL[i]);			//[mu A] -> [mu A/cm^2]
		}

		if (timeStep % Params::TIME_POINTS_SAMPLING_RATE == 0) {
			double seg = Params::PM4_PEAK;

			/*double a_muscle = (F_PE_K[seg] + F_PE_D[seg] + 2 * F_CE[seg] + F_OE[seg]) / m[seg];
			outputFile << t << "," << v[seg] + Params::VM << "," << ca_shell[seg] << "," << ca_cell[seg] << "," << l_muscle[seg] * 1E4 << "," << z05[seg] * 1E4 << "," << v_muscle[seg] << ","
					   //<< 2 * u_muscle[seg] << "," << 0.5 * a_muscle << "," << F_PE_K05[seg] << "," << F_PE_D05[seg] << "," << F_CE05[seg] << "," << F_OE05[seg] << "," << Q[seg][0] << "," << Q[seg][1] << "," << Q[seg][2] << ","
					   << u_sarc[seg] << "," << a_muscle << "," << F_PE_K[seg] << "," << F_PE_D[seg] << "," << 2 * F_CE[seg] << "," << F_OE[seg] << "," << Q[seg][0] << "," << Q[seg][1] << "," << Q[seg][2] << ","
					   << M[seg] << "," << Mp[seg] << "," << AM[seg] << "," << AMp[seg] << "," << K1[seg] << "," << K2vec[seg] << ","
					   << betaVec[seg] << "," << muscle_phi10vec[seg] << "," << muscle_phi2s01vec[seg] << "," << muscle_phi2s02vec[seg] << "," << muscle_activation[seg] << ","
					   << CF[seg] << "," << MLCP_P[seg] << "," << n_dist[seg].mean << "," << n_dist[seg].standard_deviation << "," << n0l_dist[seg].mean << "," << n0l_dist[seg].standard_deviation << endl;*/

			/*double Eca = 1000 * Params::R * Params::T / (Params::Z * Params::F) * log(Params::CA_OUT / ca_shell05_nodes[seg]);
			outputFile << t << "," << v[seg] + Params::VM << "," << amp_cca1_ca[seg] << "," << amp_cca1_na[seg] << "," << amp_egl19_ca[seg] << "," << amp_egl19_na[seg] << "," << amp_exp2[seg] << "," << ca_shell[seg] << "," << ca_cell[seg] << "," << ca_shell05_nodes[seg] << ","
					   << m0_cca1 << "," << h0_cca1 << "," << m0_egl19 << "," << square(m0_egl19) << "," << f0_egl19 << "," << quad(m0_exp2) << "," << quad(h0_exp2) << "," << gca05[seg] << "," << (segR[seg] - 2 * l_muscle[seg]) * 1E4 << "," << 2 * l_muscle[seg] * 1E4 << "," << Eca << endl;*/

			outputFile << t - 1 << "," << part1.position * 1E4 << "," << part2.position * 1E4 << "," << part3.position * 1E4 << "," << part4.position * 1E4 << "," << part5.position * 1E4 << endl;

			/*for (int seg = 0; seg < Params::N; seg++) {
				//outputFile << std::setprecision(10) << (r[seg] - 2 * l_muscle[seg]) / r[seg] << ",";
				//outputFile << std::setprecision(10) << (r[seg] - 2 * l_muscle[seg]) * 1E4 << ",";
				outputFile << std::setprecision(10) << z[seg] * 1E4 << ",";
				//outputFile << std::setprecision(10) << l_muscle05[seg] * 1E4 << ",";
				//outputFile << std::setprecision(10) << n0l_dist[seg].mean << ",";
				//outputFile << std::setprecision(10) << v[seg] + Params::VM << ",";
				//outputFile << std::setprecision(10) << ca_shell05_nodes[seg] << ",";
				//outputFile << std::setprecision(10) << ca_shell[seg] << ",";
				//outputFile << std::setprecision(10) << ca_cell[seg] << ",";
				//outputFile << std::setprecision(10) << CF[seg] << ",";
				//outputFile << std::setprecision(10) << F_CE05[seg] << ",";
				//outputFile << std::setprecision(10) << AMp[seg] << ",";
				//outputFile << std::setprecision(10) << Mp[seg] << ",";
				//outputFile << std::setprecision(10) << v_center[seg] * 1E4 << ",";
				//outputFile << std::setprecision(10) << v_shell_right[seg] * 1E18 << ",";
				//outputFile << std::setprecision(10) << flow[seg] << ",";
				//outputFile << std::setprecision(10) << amp_ca05[seg] << ",";
			}
			outputFile << endl;*/
		}

		/* Calculate V(timeStep+1) */

		BuildCableFunctionTimeIntegration::integrate2(phiV, psiV, etaV, g, alV, arV);
		BuildCableFunctionTimeIntegration::integrate3(arV, v, phiV, d, x, t, AP_t0, aNow, BV);
		//if (t > AP_t0 + Params::EX_INPUT_DURATION) {
		/*if (t + 1 > AP_t0 + Params::EX_INPUT_DURATION) {
			aOld = 0.0;
		}*/
		v_prev = v;
		TriMatOperations::triMatSolve2N1(alV, BV, v);
		//For voltage clamp, use the following loop instead of the line above:
		/*for (int seg = 0; seg < Params::N + 1; seg++) {
			v[seg] = 123.0;
		}*/

		/* Calculate Ca(timeStep+1) */

		BuildCalciumFunctionTimeIntegration::integrate2(muCa_shell, l_muscle05, v_cell05, nuCa_left, nuCa_right, a_x_cell05_left, a_x_cell05_right, muCa_cell, psiCa_shell, etaCa_shell, alCaShell, arCaShell, 1);
		BuildCalciumFunctionTimeIntegration::integrate2(muCa_shell, l_muscle05, v_cell05, nuCa_left, nuCa_right, a_x_cell05_left, a_x_cell05_right, muCa_cell, psiCa_cell,  etaCa_cell,  alCaCell,  arCaCell,  0);
		BuildCalciumFunctionTimeIntegration::integrate3(params, tau_ca_sec, arCaShell, chiCa_left, chiCa_right, amp_ca05, muCa_shell, muCa_cell, ca_shell, ca_shell_prev, ca_cell, ca_cell_prev, l_muscle05, v_cell05, rhoCa, BCa_shell, 1);
		BuildCalciumFunctionTimeIntegration::integrate3(params, tau_ca_sec, arCaCell,  chiCa_left, chiCa_right, amp_ca05, muCa_shell, muCa_cell, ca_shell, ca_shell_prev, ca_cell, ca_cell_prev, l_muscle05, v_cell05, rhoCa, BCa_cell,  0);

		ca_shell_prev = ca_shell;
		TriMatOperations::triMatSolve2N(alCaShell, BCa_shell, ca_shell);
		ca_cell_prev = ca_cell;
		TriMatOperations::triMatSolve2N(alCaCell, BCa_cell, ca_cell);

		/* Calculate pharyngeal forces at timeStep+1 and update l_muscle and z accordingly */

		/* Time-dependent variables:
		 * muscle_activation				Fraction of myosin heads available for cross-bridges in the CE [0..1], at (timeStep)	[unit-less]
		 * etaL, psiL						[gr·cm/s]
		 * ksiL								[gr·cm/s^2]
		 * chiL								[gr·cm]
		 * phiL								[cm]
		 * muscle_phi2s, muscle_phi2s0		[1/s]
		 * muscle_phi1, muscle_phi10		[1/s]
		 */

		double DT_SEC = 1E-3 * Params::DT;
		double K6, K5, muscle_phi1, muscle_phi10;
		double muscle_phi2s[2], muscle_phi2s0[2];
		F_CE_prev = F_CE;
		F_OE_prev = F_OE;
		for (int i = 0; i < Params::N; i++) {
			K6 = K1[i];
			K5 = K2vec[i];
			MuscleFunctions::buildDistributionsFromNdistMoments(i, Q[i], n_dist[i], n0l_dist[i]);
			for (int moment = 0; moment < 3; moment++) {
				betaVec[i] = MuscleFunctions::calcBeta(moment, M[i], Mp[i]);
				muscle_phi1 = MuscleFunctions::calcMusclePhi1(moment, Q[i][0], n0l_dist[i], M[i], Mp[i]);
				MuscleFunctions::calcMusclePhi2(moment, Q[i][0], n_dist[i], AM[i], AMp[i], muscle_phi2s);
				Q_estim[i][moment] = Q[i][moment] + DT_SEC * (muscle_activation[i] * betaVec[i] - muscle_phi1 - (muscle_phi2s[0] + muscle_phi2s[1]));
				if (moment > 0) {
					Q_estim[i][moment] -= DT_SEC * moment * u_sarc[i] * Q[i][moment - 1];
				} else {	//(moment == 0)
					muscle_phi10 = muscle_phi1;
					muscle_phi2s0[0] = muscle_phi2s[0];
					muscle_phi2s0[1] = muscle_phi2s[1];
				}
			}
			F_PE_K[i] = PE_K[i] * (l_muscle[i] - segR[i]);
			F_PE_D[i] = - PE_D[i] * v_muscle[i];
			F_CE[i] = muscle_activation[i] * CE_num_AM[i] * CE_K[i] * square(Params::H) / (2 * Params::L_MYOSIN_BINDING_SITES) * Q[i][1];
			F_OE[i] = MuscleFunctions::calcF_OE(i, l_muscle, segR);

			M[i] = M[i] + DT_SEC * (- K1[i] * M[i] + K2vec[i] * Mp[i] + muscle_phi2s0[1]);
			AMp[i] = AMp[i] + DT_SEC * (muscle_phi10 - muscle_phi2s0[0] - K5 * AMp[i] + K6 * AM[i]);
			muscle_phi10vec[i] = muscle_phi10;
			muscle_phi2s01vec[i] = muscle_phi2s0[0];
			muscle_phi2s02vec[i] = muscle_phi2s0[1];
		}

		// Calculate relevant variables at (timeStep+1), which becomes (timeStep) at the next time step
		Q = Q_estim;
		v_muscle_prev = v_muscle;
		for (int i = 0; i < Params::N; i++) {
			v_muscle[i] = v_muscle[i] + DT_SEC * (F_PE_K[i] + F_PE_D[i] + 2 * F_CE[i] + F_OE[i]) / m[i];								//LHS: estimated v at (timeStep+1); RHS: v at the current (timeStep)
																																		//All eq.'s terms regard a whole muscle <=> sarcomere
			u_sarc[i] = 0.5 * v_muscle[i] / Params::H;																					//LHS and RHS; estimated velocities at (timeStep+1)

			K1[i] = MuscleFunctions::calcK1(i, ca_cell[i], &print);
			MLCP_P[i] = MuscleFunctions::calcMLCP_P(i, ca_cell_prev[i], CF[i], MLCP_P[i], &print);										//NOTE: use CF at (timeStep) in order to calculate MLCP_P at (timeStep+1)
			K2vec[i] = Params::K2_MAX * MLCP_P[i];
			CF[i] = MuscleFunctions::calcMuscleContractionFactor(t, AP_t0, i, &print, x[i]);											//now calculate CF at (timeStep+1), for the next time step

			// Estimate F_CE05 and F_OE05 at (timeStep+0.5)
			F_CE05[i] = 0.5 * (3 * F_CE[i] - F_CE_prev[i]);
			F_OE05[i] = 0.5 * (3 * F_OE[i] - F_OE_prev[i]);
		}
		BuildLumenFunctionTimeIntegration::integrate2(psiL, etaL, alL, arL);
		BuildLumenFunctionTimeIntegration::integrate3(arL, l_muscle, chiL, phiL, F_CE05, F_OE05, rhoL, v_muscle_prev, v_muscle, BL);	//NOTE: l_muscle and v_muscle_prev are at (timeStep), and v_muscle is at (timeStep+1)

		l_muscle_prev = l_muscle;
		TriMatOperations::triMatSolve2N(alL, BL, l_muscle);

		for (int i = 0; i < Params::N; i++) {
			AM[i] = Q[i][0] - AMp[i];
			muscle_activation[i] = MuscleFunctions::calcMuscleActivation(segR[i], l_muscle[i]);											//can calculate muscle_activation => Mp, only after estimating l_muscle at (timeStep+1)
			Mp[i] = muscle_activation[i] * Params::L_MYOSIN_BINDING_SITES / Params::H - Q[i][0] - M[i];
			if (Mp[i] < 0) {Mp[i] = 0.0;}																								// drops below 0 at the beginning of simualtion
		}

		NodeSegmentConvertions::calcCaShellNodes(ca_shell, v_shell_left, v_shell_right, ca_shell_nodes);								//ca_shell_nodes at (timStep+1)
		NodeSegmentConvertions::calcZ(segR, l_muscle_prev, z);																			//z at (timeStep)
		NodeSegmentConvertions::calcZNodes(x, segX, z, z_nodes);																		//z_nodes at (timeStep)
		// NOTE: the following variables are estimated for (timeStep+1.5), by using their values that were calculated at (timeStep) and at (timeStep+1) (asserted to X_prev and X, respectively)
		// Hence, these variables will contain the values at (timeStep+0.5) at the next iteration, as desired
		NodeSegmentConvertions::calcCaShell05Nodes(ca_shell, ca_shell_prev, v_shell_left, v_shell_right, ca_shell05_nodes);
		NodeSegmentConvertions::calcLMuscle05Seg(l_muscle, l_muscle_prev, l_muscle05);
		NodeSegmentConvertions::calcZ(segR, l_muscle05, z05);
		NodeSegmentConvertions::calcZNodes(x, segX, z05, z05_nodes);
		NodeSegmentConvertions::calcAxCell05Left (r, z_nodes_max, z05_nodes, a_x_cell05_left);
		NodeSegmentConvertions::calcAxCell05Right(r, z_nodes_max, z05_nodes, a_x_cell05_right);
		NodeSegmentConvertions::calcVCell05(segL, r, z_nodes_max, z05_nodes, v_cell05);

		/* Calculate particle motion at timeStep+1 and update particle's parameters accordingly */

		for (int i = 0; i < 2 * Params::N; i++) {
			v_lumen_prev[i] = v_lumen[i];
		}

		FlowFunctions::calcLumenArea(z_nodes_max, z_nodes, z_max, z, a_lumen, &print);
		FlowFunctions::calcLumenVolume(x, segL, z_nodes_max, z_nodes, z_max, z, v_lumen, &print);
		FlowFunctions::calcFlow(segL, v_lumen, v_lumen_prev, flow);
		FlowFunctions::calcVelocity(z_nodes_max, z_max, a_lumen, flow, 0, vel_max);
		FlowFunctions::calcVelocity(z_nodes_max, z_max, a_lumen, flow, 1, vel_mean);
		double v_coef_part1 = (Params::FLAT_META_CO == 0) ? 3.95 : 4.2,
			   v_coef_part2 = (Params::FLAT_META_CO == 0) ? 1.68 : 1.73,
			   v_coef_part3 = (Params::FLAT_META_CO == 0) ? 1.56 : 1.5,
			   v_coef_part4 = (Params::FLAT_META_CO == 0) ? 1.4 : 1.37,
			   v_coef_part5 = 1.4;
		double part3_t0 = 85;
		double part1_dt = 96.4;
		NodeSegmentConvertions::setZAll(z_nodes, z, z_all);
		if (t >= part3_t0 + part1_dt) {FlowFunctions::calcParticlePosition(timeStep, x_all, z_all, v_coef_part1, vel_mean, 1, &part1);}
		if (t >= part3_t0 + 45.5) {FlowFunctions::calcParticlePosition(timeStep, x_all, z_all, v_coef_part2, vel_mean, 1, &part2);}
		FlowFunctions::calcParticlePosition(timeStep, x_all, z_all, v_coef_part3, vel_mean, 1, &part3);
		FlowFunctions::calcParticlePosition(timeStep, x_all, z_all, v_coef_part4, vel_mean, 1, &part4);
		FlowFunctions::calcParticlePosition(timeStep, x_all, z_all, v_coef_part5, vel_mean, 1, &part5);
		/* Note: in order to use vel_mean (vel_max), set the next field to 1 (0), e.g.:
		FlowFunctions::calcParticlePosition(timeStep, x_all, z05_all, v_coef_part5, vel_max, 0, &part5);*/
	}
}

int main() {
	Params params;
	for (int i = 0; i < Params::NUM_GENES; i++) {
		params.genes[i] = GENE_RANGES[i].known;
	}
	runPharynxSimulator(params);
	std::cout << "Simulation finished successfully" << std::endl;
}
