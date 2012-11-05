/*
 * main.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: shubham
 */
//We have assumed that the value of gama of different phases 
//remains constant in the entire domain

#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<typeinfo>
#include <cxxabi.h>
#include <map>
#include "Cell.h"
#include "Godunov.h"
#include "HLLC.h"
#include "FreeFunctions.h"

int main() {
	double dx = 1.0 / 100;
    double Pinf_water = 1.1645*1000000000;
    double gama_water = 1.932;

	std::vector<CSE::Cell> cell(2);
     
	std::map<std::string, double> leftCell_phase1, leftCell_phase2, rightCell_phase1, rightCell_phase2;
	leftCell_phase1["gama"] = 1.4; 
	leftCell_phase1["density"] = 112.99;
	leftCell_phase1["pressure"] = 1e+07;
	leftCell_phase1["velocity"] = 0;
	leftCell_phase1["phi"] = 0.0;
    leftCell_phase2["Pinf"] = 0;

	leftCell_phase2["gama"] = gama_water;
	leftCell_phase2["density"] = 1025;
	leftCell_phase2["pressure"] = 1e+07;
	leftCell_phase2["velocity"] = 0;
	leftCell_phase2["phi"] = 1-leftCell_phase1["phi"];
    leftCell_phase2["Pinf"] = Pinf_water;

	cell[0].addPhase(leftCell_phase1);
	cell[0].addPhase(leftCell_phase2);

	rightCell_phase1["gama"] = 1.4;
	rightCell_phase1["density"] = 56.49;
	rightCell_phase1["pressure"] = 5e+06;
	rightCell_phase1["velocity"] = 0;
	rightCell_phase1["phi"] = 1.0;
    rightCell_phase2["Pinf"] = 0;

	rightCell_phase2["gama"] = gama_water;
	rightCell_phase2["density"] = 1025;
	rightCell_phase2["pressure"] = 5e+06;
	rightCell_phase2["velocity"] = 0;
	rightCell_phase2["phi"] = 1-rightCell_phase1["phi"];
    rightCell_phase2["Pinf"] = Pinf_water;

	cell[1].addPhase(rightCell_phase1);
	cell[1].addPhase(rightCell_phase2);
	//for (unsigned i = 0; i < 30; ++i)
	//	cell[i].addPhase(leftCell_phase1);
	//for (unsigned i = 30; i < 102; ++i)
	//	cell[i].addPhase(rightCell_phase1);
//	for (unsigned i = 0; i < 50; ++i)
//			cell[i].addPhase(1.4, 1, 0.4, -2.0, 0);
//    for (unsigned i = 50; i < 102; ++i)
//			cell[i].addPhase(1.4, 1, 0.4, 2.0, 0);
	for (unsigned i = 0; i < cell.size(); ++i)
		cell[i].computeConsVar();
//
//	CSE::FluxComputation *fluxMethod = new CSE::Godunov(cell, 0.2);
//	fluxMethod->Solve(0.2);
//	fluxMethod->copyResult(cell);
//	std::string fluxName;
//
//	fluxName = variableType<CSE::FluxComputation>(*fluxMethod);
//	CSE::printPhaseValues(fluxName, cell);

    CSE::RiemannSolver solver(cell[0], cell[1], 2, false);
	solver.multiphaseSolve();
    std::cout << "Over" << std::endl;
	return 1;
}
