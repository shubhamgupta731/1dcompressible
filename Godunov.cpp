/*
 * Godunov.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: shubham
 */

#include "Godunov.h"

namespace CSE {

void Godunov::Solve(const double& CFL) {
	//dt_ = 0.00005;
	computeDT(CFL);
	if (getCell()[0].getPhase().size() == 1)
		SinglePhase(CFL);
	else
		MultiPhase(CFL, getCell()[0].getPhase().size());
	BoundaryConditions();
}

void Godunov::SinglePhase(const double & CFL) {
	double t = 0;
	while (t < getMax()) {
		BoundaryConditions();
		for (unsigned i = 0; i < getCell().size(); ++i) {
			if (i == getCell().size() - 2)
				continue;
			RiemannSolver riemann(getCell()[i % getCell().size()],
					getCell()[(i + 1) % getCell().size()], 1);
			riemann.solve();
			riemann.phase_star_[0].Sample(0, setCell()[i].setPhase()[0].setU_half(),
					setCell()[i].setPhase()[0].setP_half(),
					setCell()[i].setPhase()[0].setD_half(), riemann.phase_star_[0].pstar[0]
                                                          , riemann.phase_star_[0].ustar[0]);
			setCell()[i].setPhase()[0].InterCellFlux();
		}
		computeDT(CFL);
		timeIntegration();
		for (unsigned i = 0; i < getCell().size(); ++i) {
			setCell()[i].setPhase()[0].computeScalFromCons();
		}
		t += getDt();
	}
}

void Godunov::MultiPhase(const double & CFL, unsigned phase) {
	double t = 0;
	Cell leftCell;
	Cell rightCell;
	bool phaseInteraction;
	while (t < getMax()) {
		BoundaryConditions();
		for (unsigned i = 0; i < getCell().size(); ++i) {
            phaseInteraction = false;
			leftCell = getCell()[i % getCell().size()];
			rightCell = getCell()[(i + 1) % getCell().size()];
			if (i == getCell().size() - 2)
				continue;
			for (unsigned j = 0; j < phase; ++j)
				if (leftCell.getPhase()[j].getPhi()
						!= rightCell.getPhase()[j].getPhi()) {
					phaseInteraction = true;
					break;
				}
			if (!phaseInteraction) {
				for (unsigned j = 0; j < phase; ++j)
					detachedRiemann(leftCell, rightCell, i, j);
			}
		}
		computeDT(CFL);
		timeIntegration();
		for (unsigned i = 0; i < getCell().size(); ++i) {
			setCell()[i].computeScalFromCons();
		}
		t += getDt();
	}
}

void Godunov::detachedRiemann(Cell const & leftCell, Cell const & rightCell, const unsigned i,
		unsigned phase) {

	RiemannSolver riemann(leftCell, rightCell, phase);
	riemann.solve();

	riemann.phase_star_[0].Sample(0, setCell()[i].setPhase()[phase].setU_half(),
			                         setCell()[i].setPhase()[phase].setP_half(),
                                   	 setCell()[i].setPhase()[phase].setD_half(),
                                     riemann.phase_star_[0].pstar[0],
                                     riemann.phase_star_[0].ustar[0]);

	setCell()[i].setPhase()[phase].InterCellFlux();
}

void Godunov::multiphaseRiemann(const Cell& leftCell, const Cell& rightCell,
		const unsigned i, unsigned phase) {
	RiemannSolver riemann(leftCell, rightCell, phase);

}

/* namespace CSE */
}

