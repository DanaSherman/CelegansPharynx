#include "tri_mat_operations.h"
#include <iostream>
#include <stdlib.h>

void TriMatOperations::triMatFactorN1(const TriMatN1& in_mat, FactoredTriMatN1& out_mat) {
	out_mat.t_2[0] = in_mat[0];
	out_mat.t_1[0] = in_mat[1];
	out_mat.t_3[0] = in_mat[2] / out_mat.t_2[0];
	for (int i = 1; i < Params::N; i++) {
		out_mat.t_2[i] = in_mat[3 * i] - out_mat.t_3[i - 1] * out_mat.t_1[i - 1];
		out_mat.t_1[i] = in_mat[3 * i + 1];
		out_mat.t_3[i] = in_mat[3 * i + 2] / out_mat.t_2[i];
	}
	out_mat.t_2[Params::N] = in_mat[3 * Params::N] - out_mat.t_3[Params::N - 1] * out_mat.t_1[Params::N - 1];
}

void TriMatOperations::triMatFactorN(const TriMatN& in_mat, FactoredTriMatN& out_mat) {
	out_mat.t_2[0] = in_mat[0];
	out_mat.t_1[0] = in_mat[1];
	out_mat.t_3[0] = in_mat[2] / out_mat.t_2[0];
	for (int i = 1; i < Params::N - 1; i++) {
		out_mat.t_2[i] = in_mat[3 * i] - out_mat.t_3[i - 1] * out_mat.t_1[i - 1];
		out_mat.t_1[i] = in_mat[3 * i + 1];
		out_mat.t_3[i] = in_mat[3 * i + 2] / out_mat.t_2[i];
	}
	out_mat.t_2[Params::N - 1] = in_mat[3 * (Params::N - 1)] - out_mat.t_3[Params::N - 2] * out_mat.t_1[Params::N - 2];
}

void TriMatOperations::triMatSolveN1(const FactoredTriMatN1& mat, const VectorN1& v, VectorN1& out_vec) {
	VectorN1 temp{};
	temp[0] = v[0];
	for (int i = 1; i < Params::N + 1; i++) {
		temp[i] = v[i] - mat.t_3[i - 1] * temp[i - 1];
	}
	out_vec[Params::N] = temp[Params::N] / mat.t_2[Params::N];
	for (int i = Params::N - 1; i >= 0; i--) {
		out_vec[i] = (temp[i] - mat.t_1[i] * out_vec[i + 1]) / mat.t_2[i];
	}
}

void TriMatOperations::triMatSolveN(const FactoredTriMatN& mat, const VectorN& v, VectorN& out_vec) {
	VectorN1 temp{};
	temp[0] = v[0];
	for (int i = 1; i < Params::N; i++) {
		temp[i] = v[i] - mat.t_3[i - 1] * temp[i - 1];
	}
	out_vec[Params::N - 1] = temp[Params::N - 1] / mat.t_2[Params::N - 1];
	for (int i = Params::N - 2; i >= 0; i--) {
		out_vec[i] = (temp[i] - mat.t_1[i] * out_vec[i + 1]) / mat.t_2[i];
	}
}

void TriMatOperations::triMatSolve2N1(const TriMatN1& mat, const VectorN1& v, VectorN1& out_vec) {
	FactoredTriMatN1 factored{};
	triMatFactorN1(mat, factored);
	triMatSolveN1(factored, v, out_vec);
}

void TriMatOperations::triMatSolve2N(const TriMatN& mat, const VectorN& v, VectorN& out_vec) {
	FactoredTriMatN factored{};
	triMatFactorN(mat, factored);
	triMatSolveN(factored, v, out_vec);
}

void TriMatOperations::triMatVecMultN1(const TriMatN1& mat, const VectorN1& v, VectorN1& out_vec) {
	int k = 3;
	out_vec[0] = mat[0] * v[0] + mat[1] * v[1];
	for (int i = 1; i < Params::N; i++) {
		out_vec[i] = mat[k - 1] * v[i - 1] + mat[k] * v[i] + mat[k + 1] * v[i + 1];
		k += 3;
	}
	out_vec[Params::N] = mat[k - 1] * v[Params::N - 1] + mat[k] * v[Params::N];
}

void TriMatOperations::triMatVecMultN(const TriMatN& mat, const VectorN& v, VectorN& out_vec) {
	int k = 3;
	out_vec[0] = mat[0] * v[0] + mat[1] * v[1];
	for (int i = 1; i < Params::N - 1; i++) {
		out_vec[i] = mat[k - 1] * v[i - 1] + mat[k] * v[i] + mat[k + 1] * v[i + 1];
		k += 3;
	}
	out_vec[Params::N - 1] = mat[k - 1] * v[Params::N - 2] + mat[k] * v[Params::N - 1];
}
