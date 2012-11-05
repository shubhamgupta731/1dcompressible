/*
 * HLL.h
 *
 *  Created on: May 18, 2012
 *      Author: shubham
 */

#ifndef HLLC_H_
#define HLLC_H_

#include <cmath>
#include "RiemannSolver.h"
#include "FluxComputation.h"

namespace CSE {

class HLLC: public CSE::FluxComputation {
public:
	HLLC(std::vector<Cell> const & cell, double const & t_max) :
			FluxComputation(cell, t_max) {
	}

	void HLLC_Flux(Cell & lcell, Cell & rcell, double const & ustar,
			double const & pstar, RiemannSolver & riemann, const unsigned phase);
	void Solve(double const & CFL);

};

} /* namespace CSE */
#endif /* HLL_H_ */
