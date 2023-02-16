#ifndef BASE_H_
#define BASE_H_

#include <array>
#include <cmath>
#include <iostream>

template <size_t N1>
using Vector = std::array<double, N1>;
template <size_t N1>
using boolVec = std::array<bool, N1>;

template <size_t N1>
std::ostream& operator<<(std::ostream& stream, const Vector<N1>& vec) {
	stream << "(";
	for (int i = 0; i < N1 - 1; i++) {
		stream << vec[i] << ", ";
	}
	stream << vec[N1 - 1] << ")";
	return stream;
}
using MomentsVec = std::array<double, 3>;
struct Distribution {
	double mean;
	double standard_deviation;
};
struct Particle {
	double diameter;	// [cm]
	double position;	// Position along the x-axis, of a particle that flows at the center of the lumen [cm]
	int seg;
};

constexpr double square(double x) { return x * x; }
constexpr double cube(double x) { return x * x * x; }
constexpr double quad(double x) { return square(square(x)); }

enum class GeneType : int {
	cca1_mInf_A2 = 0,
	cca1_mInf_v0 = 1,
	cca1_mInf_dv = 2,
	cca1_mTau_act_A1 = 3,
	cca1_mTau_act_A2 = 4,
	cca1_mTau_act_v0 = 5,
	cca1_mTau_act_dv = 6,
	cca1_mTau_deact = 7,
	cca1_hInf_v0 = 8,
	cca1_hInf_dv = 9,
	cca1_hTau_A1 = 10,
	cca1_hTau_A2 = 11,
	cca1_hTau_v0 = 12,
	cca1_hTau_dv = 13,
	egl19_mInf_v0 = 14,
	egl19_mInf_dv = 15,
	egl19_mTau_A1 = 16,
	egl19_mTau_A2 = 17,
	egl19_mTau_v0 = 18,
	egl19_mTau_dv = 19,
	egl19_fInf_A1 = 20,
	egl19_fInf_ca05 = 21,
	egl19_fInf_dca = 22,
	egl19_fTau_A1 = 23,
	egl19_fTau_A2 = 24,
	egl19_fTau_ca05 = 25,
	egl19_fTau_dca = 26,
	exp2_mInf_v0 = 27,
	exp2_mInf_dv = 28,
	exp2_mTau_act_A1 = 29,
	exp2_mTau_act_A2 = 30,
	exp2_mTau_act_v0 = 31,
	exp2_mTau_act_dv = 32,
	exp2_mTau_act_pow = 33,
	exp2_hInf_v0 = 34,
	exp2_hInf_dv = 35,
	exp2_hTau_inact_A1 = 36,
	exp2_hTau_inact_A2 = 37,
	exp2_hTau_inact_v0 = 38,
	exp2_hTau_inact_dv = 39,
	exp2_mTau_deact_A1 = 40,
	exp2_mTau_deact_A2 = 41,
	exp2_mTau_deact_v0 = 42,
	exp2_mTau_deact_dv = 43,
	exp2_hTau_deinact = 44,
	p_cca1_ca = 45,				//cca-1  maximal calcium permeability (or flow rate) [cm/s]
	g_cca1_na = 46,				//cca-1  maximal sodium conductance 				 [m Siemens/cm^2]
	p_egl19_ca = 47,			//egl-19 maximal calcium permeability (or flow rate) [cm/s]
	g_egl19_na = 48,			//egl-19 maximal sodium conductance 				 [m Siemens/cm^2]
	g_exp2 = 49,				//maximal potassium conductance 					 [m Siemens/cm^2]
	ca_in_tau = 50,				//Time constant of intracellular calcium decay		 [ms]
};

enum class GatePart0 {
	m,
	h,
	f,
};

enum class GatePart1 {
	cca1,
	egl19,
	exp2,
	x,
};

struct AlphaBeta {
	double a_prob;
	double b_prob;
};

#endif /* BASE_H_ */
