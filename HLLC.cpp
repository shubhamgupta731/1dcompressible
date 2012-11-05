/*
 * HLL.cpp
 *
 *  Created on: May 18, 2012
 *      Author: shubham
 */

#include "HLLC.h"

void CSE::HLLC::HLLC_Flux(Cell& lcell, Cell& rcell, const double& ustar,
		const double& pstar, RiemannSolver & riemann, const unsigned phase) {
	double S_l, Sstar, S_r, q_l, q_r;
	std::vector<double> Ustar_l(riemann.consLStarVar()), Ustar_r(
			riemann.consRStarVar());

	Sstar = ustar;
	q_l = (pstar <= lcell.getPhase()[phase].getP() ?
			1 :
			std::sqrt(
					1
							+ lcell.getPhase()[phase].getG2()
									* (pstar / lcell.getPhase()[phase].getP()
											- 1)));
	q_r = (pstar <= rcell.getPhase()[phase].getP() ?
			1 :
			std::sqrt(
					1
							+ rcell.getPhase()[phase].getG2()
									* (pstar / rcell.getPhase()[phase].getP()
											- 1)));

	lcell.setPhase()[phase].computeC();
	rcell.setPhase()[phase].computeC();

	S_l = lcell.getPhase()[phase].getU() - lcell.getPhase()[phase].getC() * q_l;
	S_r = rcell.getPhase()[phase].getU() + rcell.getPhase()[phase].getC() * q_r;

	lcell.setPhase()[phase].ComputeCellFlux();
	rcell.setPhase()[phase].ComputeCellFlux();

	if (S_l >= 0)
		for (unsigned i = 0; i < 3; ++i)
			lcell.setPhase()[phase].setInterCellFlux()[i] =
					lcell.getPhase()[phase].getFlux()[i];
	else {
		if (Sstar >= 0 && S_r >= 0)
			for (unsigned i = 0; i < 3; ++i)
				lcell.setPhase()[phase].setInterCellFlux()[i] =
						lcell.getPhase()[phase].getFlux()[i]
								+ S_l
										* (Ustar_l[i]
												- lcell.getPhase()[phase].getConsVar()[i]);
		if (Sstar < 0 && S_r >= 0)
			for (unsigned i = 0; i < 3; ++i)
				lcell.setPhase()[phase].setInterCellFlux()[i] =
						rcell.getPhase()[phase].getFlux()[i]
								+ S_r
										* (Ustar_r[i]
												- rcell.getPhase()[phase].getConsVar()[i]);
		if (S_r < 0)
			for (unsigned i = 0; i < 3; ++i)
				lcell.setPhase()[phase].setInterCellFlux()[i] =
						rcell.getPhase()[phase].getFlux()[i];
	}

}

void CSE::HLLC::Solve(const double& CFL) {
	double t = 0;
	//dt_ = 0.00005;
	computeDT(CFL);
	while (t < getMax()) {
		BoundaryConditions();
		for (unsigned i = 0; i < getCell().size(); ++i) {
			if (i == getCell().size() - 2)
				continue;
			RiemannSolver riemann(getCell()[i % getCell().size()],
					getCell()[(i + 1) % getCell().size()], 0);
			riemann.solve();
			HLLC_Flux(setCell()[i % getCell().size()],
					setCell()[(i + 1) % getCell().size()], riemann.getUstar(),
					riemann.getPstar(), riemann, 0);
		}
		computeDT(CFL);
		timeIntegration();
		for (unsigned i = 0; i < getCell().size(); ++i) {
			setCell()[i].setPhase()[0].computeScalFromCons();
		}
		t += getDt();
	};

	BoundaryConditions();
}
/* namespace CSE */
