#include "build_finite_element_mesh.h"

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
using namespace std;

void BuildFiniteElementMesh::buildMesh(VectorN1& x_out, VectorN1& r_out) {
	static constexpr int inputLength = Params::XY_INPUT_LENGTH;
	Vector<inputLength> input_x{}, input_r{}, x{}, r{}, segL{};

	if (Params::EQUAL_VOLUME_ELEMENTS == 0) {
		std::string input_file_name = "PharynxGeometry_XsRs_png_dxdy0.1.csv";		//Params::FLATTEN_META_CO == 0
		if (Params::FLAT_META_CO == 1) {input_file_name = "PharynxGeometry_XsRs_png_dxdy0.1_flatMetaCo.csv";}
		std::ifstream geometry_file(input_file_name);
		if(!geometry_file.is_open()) throw std::runtime_error("Could not open file\n");
		std::string line;
		int i = 0, col;
		double val;
		while (std::getline(geometry_file, line) && i < inputLength){
	        std::stringstream ss(line);
			col = 0;
			while(ss >> val){						// Extract each double in the current line
				if (col == 0){input_x[i] = val;}
				if (col == 1){input_r[i] = val;}
				// If the next token is a comma, ignore it and move on
				if(ss.peek() == ',') ss.ignore();
				col++;
			}
			i++;
		}
	    geometry_file.close();

	    double delta_x_sum = 0;
	    for (int i = 1; i < inputLength; i++){
	    	delta_x_sum += input_x[i] - input_x[i - 1];
	    }
		double mean_segL = delta_x_sum / (inputLength - 1) * 1E-4;

		x_out[0] = input_x[0] * 1E-4;
		r_out[0] = input_r[0] * 1E-4;
	    for (int i = 1; i < inputLength; i++){
			//Convert from um to cm
			x_out[i] = x_out[i - 1] + mean_segL;
			r_out[i] = input_r[i] * 1E-4;
		}
	}
	else {
		std::string input_file_name = "PharynxGeometry_XsRs_png_dxdy0.1.csv";		//Params::FLATTEN_META_CO == 0
		if (Params::FLAT_META_CO == 1) {input_file_name = "PharynxGeometry_XsRs_png_dxdy0.1_flatMetaCo.csv";}
		std::ifstream geometry_file(input_file_name);
		if (!geometry_file.is_open()) throw std::runtime_error("Could not open file\n");
		std::string line;
		int i = 0, col;
		double val;
		while (std::getline(geometry_file, line) && i < inputLength){
			std::stringstream ss(line);
			col = 0;
			while(ss >> val){						// Extract each double in the current line
				if (col == 0){input_x[i] = val;}
				if (col == 1){input_r[i] = val;}
				// If the next token is a comma, ignore it and move on
				if(ss.peek() == ',') ss.ignore();
				col++;
			}
			i++;
		}
		geometry_file.close();

		int max_r_index = 3400;
		i = max_r_index;	//maximal r => minimal dx
		x[i] = input_x[i];
		r[i] = input_r[i];
		double dx = 0.5 * (input_x[i + 1] - input_x[i]);
		segL[i] = 2 * dx;
		double volume = M_PI * square(input_r[i]) * 2 * dx;
		double x_right = input_x[i] - dx;
		int j = i - 1;
		double divX = 100, tmpX, tmpR;
		while (x_right > input_x[0] && i >= 0 && j >= 0) {
			while (j > 0 && M_PI * square(input_r[j]) * 2 * (x_right - input_x[j]) < volume) {j--;}
			i--;
			// The desired x is between input_x[j] and input_x[j + 1]
			if (M_PI * square(input_r[j]) * 2 * (x_right - input_x[j]) == volume) {
				r[i]	= input_r[j];
				x[i]	= input_x[j];
				segL[i] = 2 * (x_right - x[i]);
				x_right = x[i] - (x_right - x[i]);
				j--;
			}
			else {
				// Eventually should find a point (x,r), ~on the line that passes through (input_x[j],input_r[j]) and (input_x[j+1],input_r[j+1]), i.e., r = mx + n,
				// that satisfies: PI * r^2 * 2(x_right - (r-n)/m) = volume.
				// Get the following eq.: r^3 - (m*x_right + n) r^2 + m*volume/(2*PI) = 0
				// i.e., a 3rd degree polynomial, whose solution is straightforward.
				double m = (input_r[j + 1] - input_r[j]) / (input_x[j + 1] - input_x[j]);
				double n = input_r[j + 1] - m * input_x[j + 1];
				double a = 1, b = - m * x_right - n, c = 0, d = (m * volume) / (2 * M_PI);
				double tmp1 = - cube(b) / (27 * cube(a)) + (b * c) / (6 * square(a)) - d / (2 * a);
				double tmp2 = cube(c / (3 * a) - square(b) / (9 * square(a)));
				if (square(tmp1) + tmp2 < 0){
					for (int p = 1; p < divX + 1; p++){
						tmpX = (divX - p) / divX * input_x[j + 1] + p / divX * input_x[j];;
						tmpR = m * (tmpX - input_x[j + 1]) + input_r[j + 1];
						if (M_PI * square(tmpR) * 2 * (x_right - tmpX) >= volume) {break;}
					}
					x[i]	= tmpX;
					r[i]	= tmpR;
					segL[i] = 2 * (x_right - x[i]);
					x_right = x[i] - (x_right - x[i]);
				}
				else {
					r[i] = std::cbrt(tmp1 + std::sqrt(square(tmp1) + tmp2)) + std::cbrt(tmp1 - std::sqrt(square(tmp1) + tmp2)) - b / (3 * a);
					if (m != 0){x[i] = (r[i] - n) / m;}
					else {
						for (int p = 1; p < divX + 1; p++){
							tmpX = (divX - p) / divX * input_x[j + 1] + p / divX * input_x[j];;
							if (M_PI * square(r[i]) * 2 * (x_right - tmpX) >= volume) {break;}
						}
						x[i] = tmpX;
					}
					if (x_right == x[i]){x_right = x[i] - 0.5 * (input_x[j + 1] - x[i]);}
					else {x_right = x[i] - (x_right - x[i]);}
					segL[i] = 2 * (x[i] - x_right);
				}
			}
		}
		int min_i = i;
		i = max_r_index;	//maximal r => minimal dx
		double x_left = input_x[i] + dx;
		j = i + 1;
		while (x_left < input_x[inputLength - 1] && i < inputLength && j < inputLength) {
			while (j < inputLength && M_PI * square(input_r[j]) * 2 * (input_x[j] - x_left) < volume) {j++;}
			// The desired x is between input_x[j - 1] and input_x[j]
			i++;
			if (M_PI * square(input_r[j]) * 2 * (input_x[j] - x_left) == volume) {
				r[i]	= input_r[j];
				x[i]	= input_x[j];
				segL[i] = 2 * (x[i] - x_left);
				x_left	= x[i] + (x[i] - x_left);
				j++;
			}
			else {
				double m = (input_r[j] - input_r[j - 1]) / (input_x[j] - input_x[j - 1]);
				double n = input_r[j] - m * input_x[j];
				double a = 1, b = - m * x_left - n, c = 0, d = - (m * volume) / (2 * M_PI);
				double tmp1 = - cube(b) / (27 * cube(a)) + (b * c) / (6 * square(a)) - d / (2 * a);
				double tmp2 = cube(c / (3 * a) - square(b) / (9 * square(a)));
				if (square(tmp1) + tmp2 < 0){
					for (int p = 1; p < divX + 1; p++){
						tmpX = (divX - p) / divX * input_x[j - 1] + p / divX * input_x[j];;
						tmpR = m * (tmpX - input_x[j]) + input_r[j];
						if (M_PI * square(tmpR) * 2 * (tmpX - x_left) >= volume) {break;}
					}
					x[i]	= tmpX;
					r[i]	= tmpR;
					segL[i] = 2 * (x[i] - x_left);
					x_left	= x[i] + (x[i] - x_left);
				}
				else {
					r[i] = std::cbrt(tmp1 + std::sqrt(square(tmp1) + tmp2)) + std::cbrt(tmp1 - std::sqrt(square(tmp1) + tmp2)) - b / (3 * a);
					if (m != 0){x[i] = (r[i] - n) / m;}
					else {
						for (int p = 1; p < divX + 1; p++){
							tmpX = (divX - p) / divX * input_x[j - 1] + p / divX * input_x[j];
							if (M_PI * square(r[i]) * 2 * (tmpX - x_left) >= volume) {break;}
						}
						x[i] = tmpX;
					}
					if (x_left == x[i]){x_left = x[i] + 0.5 * (x[i] - input_x[j - 1]);}
					else {x_left = x[i] + (x[i] - x_left);}
					segL[i] = 2 * (x_left - x[i]);
				}
			}
		}
		int max_i = i;
		for (int k = 0; k < max_i - min_i + 1; k++){
			//Convert from um to cm
			x_out[k]	= x[k + min_i] * 1.E-4;
			r_out[k]	= r[k + min_i] * 1.E-4;
		}
	}
}

void BuildFiniteElementMesh::setPerimeter(const VectorN1& r, VectorN1& p_out) {
	for (int i = 0; i < Params::N + 1; i++) {
		p_out[i] = 2.0 * M_PI * r[i];
	}
}

void BuildFiniteElementMesh::setCrossSection(const VectorN1& r, VectorN1& a_out) {
	for (int i = 0; i < Params::N + 1; i++) {
		a_out[i] = M_PI * square(r[i]);
	}
}
