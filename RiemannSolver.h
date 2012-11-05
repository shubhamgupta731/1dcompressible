/*
 * Solver.h
 *
 *  Created on: Apr 17, 2012
 *      Author: shubham
 */

#ifndef RIEMANN_SOLVER_H_
#define RIEMANN_SOLVER_H_

#include "Cell.h"
#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Dense>
#include "FreeFunctions.h"

namespace CSE {

struct Phase_Riemann {
	double ul_, ur_, pl_, pr_, dl_, dr_, cl_, cr_, gama_, phi_l_, phi_r_;
    double pinf_;
    double g1_, g2_, g3_, g4_, g5_, g6_, g7_, g8_;
	std::vector<double> ustar, pstar, dstar;
	Phase_Riemann() {
	}

	Phase_Riemann(double const& ul, double const& ur, double const& pl,
			double const& pr, double const& dl, double const& dr,
			double const& cl, double const& cr, double const & phi_l,
			double const & phi_r, double const& gama, double const& pinf=0) {
		ul_ = ul;
		ur_ = ur;
		pl_ = pl;
		pr_ = pr;
		dl_ = dl;
		dr_ = dr;
		cl_ = cl;
		cr_ = cr;
		gama_ = gama;
		phi_l_ = phi_l;
		phi_r_ = phi_r;
        pinf_ = pinf;
	}

	void set(double const& ul, double const& ur, double const& pl,
			double const& pr, double const& dl, double const& dr,
			double const& cl, double const& cr, double const & phi_l,
			double const & phi_r, double const& gama, double const& pinf) {
		ul_ = ul;
		ur_ = ur;
		pl_ = pl;
		pr_ = pr;
		dl_ = dl;
		dr_ = dr;
		cl_ = cl;
		cr_ = cr;
		gama_ = gama;
		phi_l_ = phi_l;
		phi_r_ = phi_r;
        pinf_ = pinf;
	}

	inline double getC(unsigned const phase) {
		return std::sqrt(gama_ * (pstar[phase]+pinf_)/dstar[phase]);
	}

	void operator =(Phase_Riemann const & phase) {
		ustar.resize(phase.ustar.size());
		pstar.resize(phase.pstar.size());
		dstar.resize(phase.dstar.size());
		for (unsigned i = 0; i < ustar.size(); ++i) {
			ustar[i] = phase.ustar[i];
			pstar[i] = phase.pstar[i];
			dstar[i] = phase.dstar[i];
		}
		gama_ = phase.gama_;
        pinf_ = phase.pinf_;
	}

	double guessp();
	//void starpu(double& , double&);
	void Sample(double const & s, double & u, double & p, double & d, const double&, const double&);
	void initializeG(Cell const & cell, const unsigned phase);
	void prefun(double & f, double & fd, const double & p, const double & dk,
			const double & pk, const double & ck);

	void computeDensity();

};

class RiemannSolver {
	double ul_, ur_, pl_, pr_, dl_, dr_, cl_, cr_, gama_l_, gama_r_;
	double ustar, pstar, dstar_l, dstar_r;
	double CFL_;
    std::vector<double> g1_, g2_, g3_, g4_, g5_, g6_, g7_, g8_;
	unsigned phaseNumber_;
    unsigned l_,r_;
    bool detached_;
public:
	std::vector<Phase_Riemann> phase_star_;

	RiemannSolver(Cell const & leftCell, Cell const & rightCell,
			const unsigned phase, bool detached=true);

	void initializeG(Cell const & cell);

	void solve();

	void multiphaseSolve();

	double getUstar() const {
		return ustar;
	}

	double getPstar() const {
		return pstar;
	}

	void rarefactionWave(unsigned const left, unsigned const right,
			unsigned const phase, bool Side, bool Pressure);
	void shockWave(unsigned const left, unsigned const right,
			unsigned const phase, bool Side, bool Pressure);
	void SolidContact(unsigned const left_s, unsigned const right_s,
			unsigned const left_g, unsigned const right_g, bool Side);

	void Configuration2(unsigned const left, unsigned const right);

	void Configuration2dash();

	void Configuration1(unsigned const left, unsigned const right);

    void SpecialConfiguration();

	double Df_Shock(unsigned const left, unsigned const right,
			unsigned const phase, bool Side);

	double Df_Rarefaction(unsigned const left, unsigned const right,
			unsigned const phase, bool Side);

	void InitializeU0(double const & u0, double & u_p0);

	void InitializeP0(double const & p0, double & p_u0);

	void starpu_immiscible(double& , double&, double&, double&);

	std::vector<double> consLStarVar();
	std::vector<double> consRStarVar();
};

} /* namespace CSE */
#endif /* RIEMANN_SOLVER_H_ */
