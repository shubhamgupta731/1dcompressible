/*
 * FluxComputation.h
 *
 *  Created on: May 18, 2012
 *      Author: shubham
 */

#ifndef FLUXCOMPUTATION_H_
#define FLUXCOMPUTATION_H_

#include<vector>
#include<algorithm>
#include <fstream>
#include "Cell.h"
#include "RiemannSolver.h"

namespace CSE {

class FluxComputation {
	double dt_;
	double t_max_;
	std::vector<Cell> cell_;
public:
	FluxComputation(const std::vector<Cell>& cell, const double& t_max) :
			dt_(0), t_max_(t_max), cell_(cell) {
	}

	void computeDT(const double& CFL);
	void BoundaryConditions();
	void timeIntegration();
	void copyResult(std::vector<Cell>& cell);

	std::vector<Cell> & setCell() {
		return cell_;
	}

	double & setDt() {
		return dt_;
	}

	double getMax() const {
		return t_max_;
	}

	std::vector<Cell> getCell() const {
		return cell_;
	}

	double getDt() const {
		return dt_;
	}

	virtual void Solve(double const & CFL)=0;
	virtual void SinglePhase(double const & CFL)=0;
	virtual void MultiPhase(double const & CFL, unsigned phase)=0;
	virtual ~FluxComputation() {
	}

};

} /* namespace CSE */
#endif /* FLUXCOMPUTATION_H_ */
