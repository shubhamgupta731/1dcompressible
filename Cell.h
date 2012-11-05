/*
 * Cell.h
 *
 *  Created on: Apr 17, 2012
 *      Author: shubham
 */
//TODO: The Cell constructors donot check if the values
//that have been given as input satisfy the EOS or not.

#ifndef CELL_H_
#define CELL_H_

#define Phase_For_ for(unsigned i=0;i<phase_.size();++i)
#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include <string>

namespace CSE {

struct Phase {
private:
	double gama_, d_, p_, u_, T_, e_i,pinf_;
    double u_half, d_half, p_half, c_, e_half;
	double phi_, phi_half_;
	double g1_, g2_, g3_, g4_, g5_, g6_, g7_, g8_;
	std::vector<double> consVar_;
	std::vector<double> interCellFlux_;
	std::vector<double> Flux_;

public:
	void computeC();
	void computeConsVar();
	void computeScalFromCons();
	void computeE_half();
	void ComputeGamas();
	void InterCellFlux();
	void ComputeCellFlux();
	void setValues(double gama, double d, double p, double u,double Pinf=0);
    void computeT(double p, double d, double R, double Pinf=0) {
        T_ = (p + Pinf)/d/R;
    }


	std::vector<double> getInterCellFlux() const {
		return interCellFlux_;
	}

	std::vector<double> & setInterCellFlux() {
		return interCellFlux_;
	}

	const std::vector<double> & getConsVar() const {
		return consVar_;
	}

	std::vector<double> & setFlux() {
		return Flux_;
	}

	const std::vector<double> & getFlux() const {
		return Flux_;
	}

	std::vector<double> & setConsVar() {
		return consVar_;
	}

	//Getters and Setters
	double& setGama() {
		return gama_;
	}

	const double& getGama() const {
		return gama_;
	}

	double& setD() {
		return d_;
	}

	const double& getD() const {
		return d_;
	}

	double& setP() {
		return p_;
	}

	const double& getP() const {
		return p_;
	}

	double& setU() {
		return u_;
	}

	const double& getU() const {
		return u_;
	}

	double& setPinf() {
		return pinf_;
	}

	const double& getPinf() const {
		return pinf_;
	}

	double& setT() {
		return T_;
	}

	const double& getT() const {
		return T_;
	}

	double& setG1() {
		return g1_;
	}

	const double& getG1() const {
		return g1_;
	}

	double& setG2() {
		return g2_;
	}

	const double& getG2() const {
		return g2_;
	}

	double& setG3() {
		return g3_;
	}

	const double& getG3() const {
		return g3_;
	}

	double& setG4() {
		return g4_;
	}

	const double& getG4() const {
		return g4_;
	}

	double& setG5() {
		return g5_;
	}

	const double& getG5() const {
		return g5_;
	}

	double& setG6() {
		return g6_;
	}

	const double& getG6() const {
		return g6_;
	}

	double& setG7() {
		return g7_;
	}

	const double& getG7() const {
		return g7_;
	}

	double& setG8() {
		return g8_;
	}

	const double& getG8() const {
		return g8_;
	}

	double& setP_half() {
		return p_half;
	}

	const double& getP_half() const {
		return p_half;
	}

	double& setD_half() {
		return d_half;
	}

	const double& getD_half() const {
		return d_half;
	}

	double& setU_half() {
		return u_half;
	}

	const double& getU_half() const {
		return u_half;
	}

	double& setC() {
		return c_;
	}

	const double& getC() const {
		return c_;
	}

	const double & getPhi() const {
		return phi_;
	}

	double & setPhi() {
		return phi_;
	}
};

class Cell {

	std::vector<Phase> phase_;
	double dx_;

public:
	Cell() {
//		phase_.resize(1);
		Phase_For_ {
			dx_ = 0;
			phase_[i].setConsVar().resize(3);
			phase_[i].setInterCellFlux().resize(3);
			phase_[i].setFlux().resize(3);
		}
	}

	Cell(double gama, double d, double p, double u, double dx) :
			dx_(dx) {
		phase_.resize(1);
		Phase_For_ {
			phase_[i].setValues(gama, d, p, u);
			phase_[i].setConsVar().resize(3);
			phase_[i].setInterCellFlux().resize(3);
			phase_[i].setFlux().resize(3);
			phase_[i].setT() = phase_[i].getP() / phase_[i].getD() / 8.314;
			phase_[i].setC() = std::sqrt(
					phase_[i].getGama() * phase_[i].getP() / phase_[i].getD());
			phase_[i].ComputeGamas();
		}
	}

	Cell(double gama, double d, double p, double u, double T, double dx) :
			dx_(dx) {
		phase_.resize(1);
		Phase_For_ {
			phase_[i].setValues(gama, d, p, u);
			phase_[i].setConsVar().resize(3);
			phase_[i].setInterCellFlux().resize(3);
			phase_[i].setFlux().resize(3);
			phase_[i].setT() = T;
			phase_[i].setC() = std::sqrt(
					phase_[i].getGama() * phase_[i].getP() / phase_[i].getD());
			phase_[i].ComputeGamas();
		}
	}

	Cell(double gama, double d, double p, double u, double dx, double phi,
			double Pinf, unsigned phase) :
			dx_(dx) {
		if (phase + 1 > phase_.size())
			phase_.resize(phase + 1);
		phase_[phase].setValues(gama, d, p, u, Pinf);
		phase_[phase].setConsVar().resize(4);
		phase_[phase].setInterCellFlux().resize(4);
		phase_[phase].setFlux().resize(4);
        phase_[phase].computeT(phase_[phase].getP(),
                               phase_[phase].getD(),
                               8.314,
                               phase_[phase].getPinf());
		phase_[phase].computeC();
		phase_[phase].ComputeGamas();
		phase_[phase].setPhi() = phi;
	}

	//Cell(double gama, double d, double p, double u, double dx, double T,
	//		double phi, double Pinf=0, unsigned phase) :
	//		dx_(dx) {
	//	if (phase > phase_.size())
	//		phase_.resize(phase);
	//	phase_[phase].setValues(gama, d, p, u, Pinf);
	//	phase_[phase].setConsVar().resize(3);
	//	phase_[phase].setInterCellFlux().resize(3);
	//	phase_[phase].setFlux().resize(3);
	//	phase_[phase].setT() = T;
	//	phase_[phase].setC() = std::sqrt(
	//			phase_[phase].getGama() * phase_[phase].getP()
	//					/ phase_[phase].getD());
	//	phase_[phase].ComputeGamas();
	//	phase_[phase].setPhi() = phi;
	//}

	void computeConsVar();
	void computeScalFromCons();
	void operator =(const Cell & tmp);

	const std::vector<Phase> & getPhase() const {
		return phase_;
	}

	std::vector<Phase> & setPhase() {
		return phase_;
	}

	const double & getDx() const {
		return dx_;
	}

	double & setDx() {
		return dx_;
	}

	void addPhase( std::map<std::string, double> phaseData);

};

struct MaxU {
	bool operator()(Cell const & cell1, Cell const & cell2) const {
		return cell1.getPhase()[0].getU() < cell2.getPhase()[0].getU();
	}
};

struct MaxDX {
	bool operator()(Cell const & cell1, Cell const & cell2) const {
		return cell1.getDx() < cell2.getDx();
	}
};

} /* namespace CSE */
#endif /* CELL_H_ */
