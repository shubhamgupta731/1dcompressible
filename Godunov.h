/*
 * Godunov.h
 *
 *  Created on: Apr 18, 2012
 *      Author: shubham
 */

#ifndef GODUNOV_H_
#define GODUNOV_H_

#include<vector>
#include<algorithm>
#include <fstream>
#include "FluxComputation.h"

namespace CSE {

class Godunov: public FluxComputation {
public:
	Godunov(std::vector<Cell> const & cell, double const & t_max) : //t_max: Max. duration of execution
			FluxComputation(cell, t_max) {
	}

	void Solve(double const & CFL);
	void SinglePhase(double const & CFL);
	void MultiPhase(double const & CFL, unsigned phase);
	void detachedRiemann(Cell const & leftCell, Cell const & rightCell,
			const unsigned i, unsigned phase);
	void multiphaseRiemann(Cell const & leftCell, Cell const & rightCell,
			const unsigned i, unsigned phase);
};

} /* namespace CSE */
#endif /* GODUNOV_H_ */
