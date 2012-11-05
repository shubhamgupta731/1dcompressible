/*
 * FluxComputation.cpp
 *
 *  Created on: May 18, 2012
 *      Author: shubham
 */

#include "FluxComputation.h"

void CSE::FluxComputation::computeDT(const double& CFL) {
	double Smax = 0;
	for (unsigned phase = 0; phase < cell_[0].getPhase().size(); ++phase) {
		for (unsigned i = 0; i < cell_.size() - 2; ++i) {
			cell_[i].setPhase()[phase].computeC();
			if (Smax
					< std::abs(cell_[i].getPhase()[phase].getU())
							+ cell_[i].getPhase()[phase].getC())
				Smax = std::abs(cell_[i].getPhase()[phase].getU())
						+ cell_[i].getPhase()[phase].getC();
		}
	}
//	std::vector<Cell>::iterator it2 = std::min_element(cell_.begin(),
//			cell_.end(), MaxDX());
//	std::cout << cell_[0].getDx() << std::endl;
	dt_ = CFL * (cell_[0].getDx()) / Smax;
	std::cout << "Hello: " << dt_ << std::endl;
}

void CSE::FluxComputation::BoundaryConditions() {
	cell_[cell_.size() - 1] = cell_[0];

	cell_[cell_.size() - 2] = cell_[cell_.size() - 3];
}

void CSE::FluxComputation::timeIntegration() {
	for (unsigned i = 0; i < cell_.size(); ++i) {
		for (unsigned p = 0; p < cell_[i].getPhase().size(); ++p) {
			if (i == cell_.size() - 2)
				continue;
			for (unsigned j = 0; j < 3; ++j) {
				cell_[i % cell_.size()].setPhase()[p].setConsVar()[j] -=
						dt_ / (cell_[i % cell_.size()].getDx())
								* (cell_[i % cell_.size()].getPhase()[p].getInterCellFlux()[j]);
				cell_[(i + 1) % cell_.size()].setPhase()[p].setConsVar()[j] +=
						dt_ / (cell_[(i + 1) % cell_.size()].getDx())
								* (cell_[i % cell_.size()].getPhase()[p].getInterCellFlux()[j]);
			}
		}
	}
}

void CSE::FluxComputation::copyResult(std::vector<Cell>& cell) {
	cell = cell_;
}

/* namespace CSE */
