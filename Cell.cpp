/*
 * Cell.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: shubham
 */

#include "Cell.h"

namespace CSE {

void Phase::computeC() {
	c_ = std::sqrt(gama_*(p_+pinf_)/d_);
}
void Phase::computeConsVar() {
	consVar_[0] = phi_*d_;
	consVar_[1] = phi_*d_*u_;
	consVar_[2] = phi_*d_*(0.5 * std::pow(u_, 2.0) + (p_+pinf_) / (gama_ - 1) / d_);
    consVar_[3] = phi_;
}
void Phase::computeE_half() {
	e_half = d_half
			* (0.5 * std::pow(u_half, 2.0) + (p_half + pinf_) / (gama_ - 1) / d_half);
}

void Phase::ComputeGamas() {
	g1_ = (gama_ - 1.0) / (2.0 * gama_); //g1 is Eq. 9.31(i)
	g2_ = (gama_ + 1.0) / (2.0 * gama_);
	g3_ = 2.0 * gama_ / (gama_ - 1.0);
	g4_ = 2.0 / (gama_ - 1.0);
	g5_ = 2.0 / (gama_ + 1.0);
	g6_ = (gama_ - 1.0) / (gama_ + 1.0);
	g7_ = (gama_ - 1.0) / 2.0;
	g8_ = gama_ - 1.0;
}

void Phase::computeScalFromCons() {
    phi_ = consVar_[3];
	if (phi_ != 0.0){
        d_ = consVar_[0]/phi_;
        u_ = consVar_[1]/d_/phi_;
        p_ = (consVar_[2]/d_/phi_-0.5*std::pow(u_, 2.0))*(gama_-1)*d_-pinf_;
    }
}

void Phase::InterCellFlux() {
	computeE_half();
	interCellFlux_[0] = phi_*d_half * u_half;
	interCellFlux_[1] = phi_*(d_half * std::pow((u_half), 2) + p_half);
	interCellFlux_[2] = phi_*u_half * (e_half + p_half);
    //if (phi_half_ != phi_) { 
    //    //interCellFlux_[3] = u_half*(phi_half_-phi_);
    //    std::cout << "Error: Solver will not be able to handle this problem" << std::endl;
    //}
}

void Cell::computeConsVar() {
	Phase_For_
		phase_[i].computeConsVar();
}

void Cell::computeScalFromCons() {
	Phase_For_
			phase_[i].computeScalFromCons();
}

void Cell::operator =(const Cell & tmp) {
	phase_.resize(tmp.getPhase().size());
	Phase_For_ {
		phase_[i].setP()    = tmp.getPhase()[i].getP();
		phase_[i].setPinf() = tmp.getPhase()[i].getPinf();
		phase_[i].setU()    = tmp.getPhase()[i].getU();
		phase_[i].setD()    = tmp.getPhase()[i].getD();
		phase_[i].setT()    = tmp.getPhase()[i].getT();
		phase_[i].setC()    = tmp.getPhase()[i].getC();
		dx_                 = tmp.getDx();
		phase_[i].setGama() = tmp.getPhase()[i].getGama();
		phase_[i].ComputeGamas();
		phase_[i].setInterCellFlux() = tmp.getPhase()[i].getInterCellFlux();
		phase_[i].setFlux() = tmp.getPhase()[i].getFlux();
		phase_[i].setConsVar() = tmp.getPhase()[i].getConsVar();
		phase_[i].setPhi() = tmp.getPhase()[i].getPhi();
	}
}

void Phase::ComputeCellFlux() {
    double e_;
	e_ = d_*(0.5*std::pow(u_, 2.0)+(p_+pinf_)/(gama_-1)/d_);
	Flux_[0] = phi_*d_*u_;
	Flux_[1] = phi_*(d_*std::pow((u_), 2)+p_);
	Flux_[2] = phi_*u_*(e_ + p_);
}

void Cell::addPhase(std::map<std::string, double> phaseData) {
	double gama, d, p, u, phi,pinf;
	gama = phaseData["gama"];
	d = phaseData["density"];
	p = phaseData["pressure"];
	u = phaseData["velocity"];
	phi = phaseData["phi"];
    pinf = phaseData["Pinf"];
	phase_.resize(phase_.size() + 1);
	unsigned phase = phase_.size() - 1;
	phase_[phase].setValues(gama, d, p, u, pinf);
	phase_[phase].setConsVar().resize(3);
	phase_[phase].setInterCellFlux().resize(3);
	phase_[phase].setFlux().resize(3);

	if(phaseData.find("Temperature") != phaseData.end())
		phase_[phase].setT() = phaseData["Temperature"];
    phase_[phase].computeC();
	phase_[phase].ComputeGamas();
	phase_[phase].setPhi() = phi;
}
}//namespace CSE

void CSE::Phase::setValues(double gama, double d, double p, double u,double Pinf) {
	gama_ = gama;
	d_    = d;
	p_    = p;
	u_    = u;
    pinf_ = Pinf;
}

/* namespace CSE */
