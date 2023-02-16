#ifndef TRI_MAT_OPERATIONS_H_
#define TRI_MAT_OPERATIONS_H_

#include "base.h"
#include "params.h"

struct FactoredTriMatN1 {
	Vector<Params::N> t_1;
	Vector<Params::N + 1> t_2;
	Vector<Params::N> t_3;
};

struct FactoredTriMatN {
	Vector<Params::N - 1> t_1;
	Vector<Params::N> 	  t_2;
	Vector<Params::N - 1> t_3;
};

class TriMatOperations {
 public:
	static void triMatFactorN1(const TriMatN1& in_mat, FactoredTriMatN1& out_mat);
	static void triMatFactorN (const TriMatN&  in_mat, FactoredTriMatN&  out_mat);
	static void triMatSolveN1(const FactoredTriMatN1& mat, const VectorN1& v, VectorN1& out_vec);
	static void triMatSolveN (const FactoredTriMatN&  mat, const VectorN&  v, VectorN&  out_vec);
	static void triMatSolve2N1(const TriMatN1& mat, const VectorN1& v, VectorN1& out_vec);
	static void triMatSolve2N (const TriMatN&  mat, const VectorN&  v, VectorN& out_vec);
	static void triMatVecMultN1(const TriMatN1& mat, const VectorN1& v, VectorN1& out_vec);
	static void triMatVecMultN (const TriMatN&  mat, const VectorN&  v, VectorN&  out_vec);
};

#endif	// TRI_MAT_OPERATIONS_H_
