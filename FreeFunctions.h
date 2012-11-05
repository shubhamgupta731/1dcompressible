/*
 * FreeFunctions.h
 *
 *  Created on: May 22, 2012
 *      Author: shubham
 */

#ifndef FREEFUNCTIONS_H_
#define FREEFUNCTIONS_H_

#include "Cell.h"
#include<fstream>
#include<string>
#include <cxxabi.h>
#include <cstdlib>
#include <sstream>
#include <Eigen/Dense>

#define LEFT false
#define RIGHT	true
#define SOLID_SIZE 4 
#define GAS_SIZE 4
#define small std::numeric_limits<double>::denorm_min()

template<class A>
std::string variableType(const A & a) {
	int status;
	char *name = abi::__cxa_demangle(typeid(a).name(), 0, 0, &status);
	std::string Name(name);
	return Name.substr(Name.find_last_of(':') + 1);
}

namespace CSE {

void printPhaseValues(std::string & flux, const std::vector<Cell> & cell);
template<class P, class Q> std::vector<double> solve(P A, Q b);

}

template<class P, class Q>
inline std::vector<double> CSE::solve(P A, Q b) {
	std::vector<double> x_ret;
	x_ret.resize(4);
	Eigen::Matrix4d A_;
	Eigen::Vector4d b_, x;
	for(unsigned i=0;i<4;++i){
		b_(i) = b[i];
		for(unsigned j=0;j<4;++j)
			A_(i,j) = A[i*4+j];
	}
	x = A_.ldlt().solve(b_);
	for(unsigned i=0;i<4;++i)
		x_ret[i] = x(i);
	return x_ret;
}


//Namespace CSE

#endif /* FREEFUNCTIONS_H_ */
