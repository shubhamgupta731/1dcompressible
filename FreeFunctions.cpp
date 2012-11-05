/*
 * FreeFunctions.cpp
 *
 *  Created on: May 22, 2012
 *      Author: shubham
 */
#include "FreeFunctions.h"

void CSE::printPhaseValues(std::string & flux,
		const std::vector<CSE::Cell> & cell) {

	std::ofstream out;
	std::stringstream phase_type;
	for (unsigned phase = 0; phase < cell[0].getPhase().size(); ++phase) {
		phase_type.str("");
		phase_type << phase;
		std::string flux_var = flux;
		out.open(
				(((flux_var.append("_")).append(phase_type.str())).append(
						".data")).c_str());
		for (unsigned i = 0; i < cell.size() - 2; ++i)
			out << (i + 1) * cell[i].getDx() << '\t'
					<< cell[i].getPhase()[phase].getU() << '\t'
					<< cell[i].getPhase()[phase].getP() << '\t'
					<< cell[i].getPhase()[phase].getD() << std::endl;
		out.close();
	}
}
