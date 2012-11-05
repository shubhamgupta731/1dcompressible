/*
 * Solver.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: shubham
 */

#include "RiemannSolver.h"

namespace CSE {

//Constructor
RiemannSolver::RiemannSolver(const Cell& leftCell, const Cell& rightCell,
		const unsigned phase, bool detached) {
	// Phase 0 is s and phase 1 is g.
    detached_ = detached;
    if (!detached) {
        phase_star_.resize(phase);
        for (unsigned i = 0; i < phase; ++i) {
            phase_star_[i].dstar.resize(SOLID_SIZE);
            phase_star_[i].ustar.resize(SOLID_SIZE);
            phase_star_[i].pstar.resize(SOLID_SIZE);

            phase_star_[i].set(leftCell.getPhase()[i].getU(),
                    rightCell.getPhase()[i].getU(),
                    leftCell.getPhase()[i].getP(),
                    rightCell.getPhase()[i].getP(),
                    leftCell.getPhase()[i].getD(),
                    rightCell.getPhase()[i].getD(),
                    leftCell.getPhase()[i].getC(),
                    rightCell.getPhase()[i].getC(),
                    leftCell.getPhase()[i].getPhi(),
                    rightCell.getPhase()[i].getPhi(),
                    leftCell.getPhase()[i].getGama(),
                    leftCell.getPhase()[i].getPinf());
            phase_star_[i].initializeG(leftCell, i);
        }

        if (phase_star_[0].phi_l_ == 1.0){
            l_ = 0;
            r_ = 1;    
        } else {
            l_ = 1;
            r_ = 0;
        }
    }
    if (detached) {
        phase_star_.resize(1);
        phase_star_[0].dstar.resize(SOLID_SIZE);
        phase_star_[0].ustar.resize(SOLID_SIZE);
        phase_star_[0].pstar.resize(SOLID_SIZE);

        phase_star_[0].set(leftCell.getPhase()[phase].getU(),
                rightCell.getPhase()[phase].getU(),
                leftCell.getPhase() [phase].getP(),
                rightCell.getPhase()[phase].getP(),
                leftCell.getPhase() [phase].getD(),
                rightCell.getPhase()[phase].getD(),
                leftCell.getPhase() [phase].getC(),
                rightCell.getPhase()[phase].getC(),
                leftCell.getPhase() [phase].getPhi(),
                rightCell.getPhase()[phase].getPhi(),
                leftCell.getPhase() [phase].getGama(),
                leftCell.getPhase() [phase].getPinf());
        phase_star_[0].initializeG(leftCell, phase);
        l_ = 0;
        r_ = 0;
    }

}

void RiemannSolver::solve() {
    if (detached_) {
        phase_star_[0].pstar[0] = phase_star_[0].guessp();
        starpu_immiscible(phase_star_[0].pstar[0],phase_star_[0].pstar[0],
                          phase_star_[0].ustar[0],phase_star_[0].ustar[0]);
    } else{
        for(unsigned i=0; i<phase_star_.size();++i)
            phase_star_[i].pstar[0] = phase_star_[i].guessp();
        starpu_immiscible(phase_star_[0].pstar[0],phase_star_[1].pstar[0],
                          phase_star_[0].ustar[0],phase_star_[1].ustar[0]);
    }
}

void RiemannSolver::multiphaseSolve() {
	//Coniguration2
	//Configuration2dash();
    SpecialConfiguration();
}

//Initialize Gammas
void RiemannSolver::initializeG(const Cell& cell) {
    g1_.resize(phaseNumber_);
    g2_.resize(phaseNumber_);
    g3_.resize(phaseNumber_);
    g4_.resize(phaseNumber_);
    g5_.resize(phaseNumber_);
    g6_.resize(phaseNumber_);
    g7_.resize(phaseNumber_);
    for(unsigned i=0;i<phaseNumber_;++i){
        g1_[i] = cell.getPhase()[i].getG1();
        g2_[i] = cell.getPhase()[i].getG2();
        g3_[i] = cell.getPhase()[i].getG3();
        g4_[i] = cell.getPhase()[i].getG4();
        g5_[i] = cell.getPhase()[i].getG5();
        g6_[i] = cell.getPhase()[i].getG6();
        g7_[i] = cell.getPhase()[i].getG7();
        g8_[i] = cell.getPhase()[i].getG8();
    }
}

void Phase_Riemann::initializeG(const Cell& cell, const unsigned phase) {
    g1_ = cell.getPhase()[phase].getG1();
    g2_ = cell.getPhase()[phase].getG2();
    g3_ = cell.getPhase()[phase].getG3();
    g4_ = cell.getPhase()[phase].getG4();
    g5_ = cell.getPhase()[phase].getG5();
    g6_ = cell.getPhase()[phase].getG6();
    g7_ = cell.getPhase()[phase].getG7();
    g8_ = cell.getPhase()[phase].getG8();
}

//Initialize p in the star region based on the left and right cell data
double Phase_Riemann::guessp() {
	double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr, qmax, quser, um, pm;
	quser = 2.0;
	//Compute guess pressure from PVRS Riemann Solver
	cup = 0.25 * (dl_ + dr_) * (cl_ + cr_);
	ppv = 0.5 * (pl_ + pr_) + 0.5 * (ul_ - ur_) * cup; //Eq. 4.47(ii)
	ppv = std::max(0.0, ppv); //Eq. 4.47(i)
	pmin = std::min(pl_, pr_);
	pmax = std::max(pl_, pr_);
	qmax = pmax / pmin;
	if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
		pm = ppv; // select PVRS Riemann solver
	else {
		if (ppv < pmin) {
			// select Two-Rarefaction Riemann solver
			pq = std::pow(pl_ / pr_, g1_);
			um = (pq * ul_ / cl_ + ur_ / cr_ + g4_ * (pq - 1.0))
					/ (pq / cl_ + 1.0 / cr_); //Eq. 9.35
			ptl = 1.0 + g7_ * (ul_ - um) / cl_; //Eq. 9.36
			ptr = 1.0 + g7_ * (um - ur_) / cr_; //Eq. 9.36
			pm = 0.5 * (std::pow(pl_ * ptl, g3_) + std::pow(pr_ * ptr, g3_)); //Eq. 9.36
		} else {
			// select Two-Shock Riemann solver with PVRS as estimate
			gel =
					std::sqrt(
							(float) (g5_ / dl_) / (g6_ * pl_ + ppv)); //Eq. 9.31 & 9.41
			ger =
					std::sqrt(
							(float) (g5_ / dr_) / (g6_ * pr_ + ppv));

			pm = (gel * pl_ + ger * pr_ - (ur_ - ul_)) / (gel + ger); //Eq. 9.42
		}
	}
	return pm;
}

//Relevant part of code which is used to compute p_star
void RiemannSolver::starpu_immiscible(double& ptemp1, double& ptemp2
                                    , double& utemp1, double& utemp2) {
	const int nriter = 200;
	const double tolpre = 1.0e-6;
	double change, fl, fld, fr, frd, pold, udiff;
    double ul, ur, pl, pr, dl, dr, cl, cr;
	unsigned i = 1;
    double ptemp, utemp;

    ptemp = (ptemp1 + ptemp2)/2;
    utemp = (utemp1 + utemp2)/2;
    
    ul = phase_star_[l_].ul_;
    ur = phase_star_[r_].ur_;
    dl = phase_star_[l_].dl_;
    dr = phase_star_[r_].dr_;
    pl = phase_star_[l_].pl_;
    pr = phase_star_[r_].pr_;
    cl = phase_star_[l_].cl_;
    cr = phase_star_[r_].cr_;
    

	udiff = ur - ul;
	//	std::cout << "----------------------------------------\n"
	//			<< "   Iteration number     Change\n"
	//			<< "----------------------------------------" << std::endl;
	pold = ptemp;
    std::cout << "ptemp: " << ptemp << l_ << r_ << std::endl;
        std::cout << pold << '\t' << dl << '\t'<< pl << '\t'<< cl << std::endl;
	for (; i <= nriter; i++) {
		phase_star_[l_].prefun(fl, fld, pold, dl, pl, cl);
		phase_star_[r_].prefun(fr, frd, pold, dr, pr, cr);

		ptemp = pold - (fl + fr + udiff) / (fld + frd); //Eq. 4.44. Note that
		//fd=fld+frd as is given
		//in the explaination below Eq. 4.37
		//Also f=fl+fr+udiff as given in Eq. 4.36
		change = 2.0 * std::abs((ptemp - pold) / (ptemp + pold));
		//if (change > 0)
			std::cout << '\t' << i << "\t\t" << change << "\t\t" << fl << "\t\t" << fr << std::endl;

		if (change <= tolpre)
			break;

		if (ptemp < 0.0)
			ptemp = tolpre;

		pold = ptemp;
	}
	if (i > nriter)
		std::cout << "divergence in Newton-Raphson iteration" << std::endl;

	// compute velocity in star region
	utemp = 0.5 * (ul + ur + fr - fl);
	//	std::cout << "----------------------------------------\n"
	//			<< "     Pressure           Velocity\n"
	//			<< "----------------------------------------\n" << "     " << pstar
	//			<< "\t\t" << ustar << '\n'
	//			<< "----------------------------------------" << std::endl;
    ptemp1 = ptemp;
    ptemp2 = ptemp;
    utemp1 = utemp;
    utemp2 = utemp;
}

void Phase_Riemann::prefun(double& f, double& fd, const double& p,
		const double& dk, const double& pk, const double& ck) {
	// purpose: to evaluate the pressure functions
	//          fl and fr in exact Riemann solver
	//          and their first derivatives
	double ak, bk, pratio, qrt;
	if (p <= pk) {
		// rarefaction wave
		pratio = (p+pinf_)/(pk+pinf_);
		f = g4_ * ck * (std::pow(pratio, g1_) - 1.0); //Eq. 4.35(Next)
		fd = (1.0 / (dk * ck)) * std::pow(pratio, -g2_); //First Derivative
	} else {
		//  shock wave
		ak = g5_/dk;
		bk = g6_*pk+gama_*pinf_/(gama_+1);
		qrt = std::sqrt((float) (ak / (bk + p)));
		f = (p - pk) * qrt; //Eq. 4.31(Next)
		fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt; //First Derivative
	}
}

//Sampling of Solution in the relevant part of the code
//Phi is not being sampled beause its still not a feature
//of the solver.

void Phase_Riemann::Sample(const double& s, double& u, double& p
                         , double& d, const double& ptemp, const double& utemp) {
	double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;
	if (s <= utemp) {
		// sampling point lies to the left of the contact discontinuity
		if (ptemp <= pl_) {
			// left rarefaction
			shl = ul_ - cl_;
			if (s <= shl) {
				// sampled point is left data state
				d = dl_;
				u = ul_;
				p = pl_;
			} else {
				cml = cl_ * pow((ptemp+pinf_)/(pl_+pinf_), g1_);
				stl = utemp - cml;
				if (s > stl) {
					// sampled point is star left state
					d = dl_ * std::pow((ptemp+pinf_)/(pl_+pinf_), 1.0 / gama_);
					u = utemp;
					p = ptemp;
				} else {
					// sampled point is inside left fan
					u = g5_ * (cl_ + g7_ * ul_ + s);
					c = g5_ * (cl_ + g7_ * (ul_ - s));
					d = dl_ * std::pow(c / cl_, g4_);
					p = (pl_+pinf_)*std::pow(c / cl_, g3_)-pinf_;
				}
			}

		} else {
			// left shock
			pml = ptemp / pl_;
			sl = ul_ - cl_ * std::sqrt(g2_ * pml + g1_);
			if (s <= sl) {
				// sampled point is left data state
				d = dl_;
				u = ul_;
				p = pl_;
			} else {
				// sampled point is star left state
				d = dl_ * (pml+g6_+gama_/(gama_+1)*pinf_/pl_) / (pml*g6_+1.0+gama_/(gama_+1)*pinf_/pl_);
				u = utemp;
				p = ptemp;
			}
		}

	} else {
		// sampling point lies to the right of the contact discontinuity
		if (ptemp > pr_) {
			// right shock
			pmr = ptemp / pr_;
			sr = ur_ + cr_ * std::sqrt(g2_ * pmr + g1_);
			if (s >= sr) {
				// sampled point is right data state
				d = dr_;
				u = ur_;
				p = pr_;
			} else {
				// sampled point is star right state
				d = dr_ * (pmr+g6_+gama_/(gama_+1)*pinf_/pr_)/(pmr*g6_+1.0+gama_/(gama_+1)*pinf_/pr_);
				u = utemp;
				p = ptemp;
			}
		} else {
			// right rarefaction
			shr = ur_ + cr_;
			if (s >= shr) {
				// sampled point is right data state
				d = dr_;
				u = ur_;
				p = pr_;
			} else {
				cmr = cr_ * std::pow((ptemp+pinf_)/(pr_+pinf_), g1_);
				str = utemp + cmr;
				if (s <= str) {
					// sampled point is star right state
					d = dr_ * std::pow((ptemp+pinf_)/(pr_+pinf_), 1.0 / gama_);
					u = utemp;
					p = ptemp;
				} else {
					// sampled point is inside left fan
					u = g5_ * (-cr_ + g7_ * ur_ + s);
					c = g5_ * (cr_ - g7_ * (ur_ - s));
					d = dr_ * std::pow(c / cr_, g4_);
					p = (pr_+pinf_) * std::pow(c / cr_, g3_)-pinf_;
				}
			}

		}

	}

}

void Phase_Riemann::computeDensity() {
	if (pstar[0] <= pl_)
		dstar[0] = dl_ * std::pow(pstar[0] / pl_, 1.0 / gama_);
	else
		dstar[0] = dl_ * (pstar[0] / pl_ + g6_) / (pstar[0] / pl_ * g6_ + 1.0);

	if (pstar[0] > pr_)
		dstar[1] = dr_ * (pstar[0] / pr_ + g6_) / (pstar[0] / pr_ * g6_ + 1.0);
	else
		dstar[1] = dr_ * std::pow(pstar[0] / pr_, 1.0 / gama_);
}

std::vector<double> RiemannSolver::consLStarVar() {
//	computeDensity();
//	std::vector<double> consStarVar(3);
//	consStarVar[0] = dstar_l;
//	consStarVar[1] = dstar_l * ustar;
//	consStarVar[2] = dstar_l
//			* (0.5 * std::pow(ustar, 2.0) + pstar / (gama_l_ - 1) / dstar_l);
//	return consStarVar;
}

//Pressure is True
void RiemannSolver::rarefactionWave(const unsigned left, const unsigned right,
		const unsigned phase, bool Side, bool Pressure) {
	Phase_Riemann phase_star = phase_star_[phase];
	//	double a_s;
	double F_s, G_s;
	double ds, p, ps, as, us;
	double u, y;
	//Right is True
	if (Pressure) {
		if (Side) {
			//d = phase_star.dstar[left];
			ds = phase_star.dstar[right];
			p = phase_star.pstar[left];
			ps = phase_star.pstar[right];
			as = phase_star.getC(right);
			us = phase_star.ustar[right];
			F_s = phase_star.g4_ * as * (std::pow(p / ps, phase_star.g1_) - 1);
			G_s = ds * std::pow(p / ps, 1 / phase_star.gama_);
			phase_star_[phase].ustar[left] = us + F_s;
			phase_star_[phase].dstar[left] = G_s;
		} else {
			ds = phase_star.dstar[left];
			//d = phase_star.dstar[right];
			ps = phase_star.pstar[left];
			p = phase_star.pstar[right];
			as = phase_star.getC(left);
			us = phase_star.ustar[left];
			//A_s = g5_ / ds;
			//B_s = g6_ * ps;
			F_s = phase_star.g4_ * as * (std::pow(p / ps, phase_star.g1_) - 1);
			G_s = ds * std::pow(p / ps, 1 / phase_star.gama_);
			phase_star_[phase].ustar[right] = us - F_s;
			phase_star_[phase].dstar[right] = G_s;
		}
	} else {
		if (Side) {
			//d = phase_star.dstar[left];
			ds = phase_star.dstar[right];
			ps = phase_star.pstar[right];
			as = phase_star.getC(right);
			us = phase_star.ustar[right];
			u = phase_star.ustar[left];
			y = u - us;
			phase_star_[phase].pstar[left] = ps
					* std::pow(1 + (y * phase_star.g7_ / as), phase_star.g3_);
			G_s = ds
					* std::pow(phase_star_[phase].pstar[left] / ps,
							1 / phase_star.gama_);
			phase_star_[phase].dstar[left] = G_s;
		} else {
			ds = phase_star.dstar[left];
			//d = phase_star.dstar[right];
			ps = phase_star.pstar[left];
			as = phase_star.getC(left);
			us = phase_star.ustar[left];
			u = phase_star.ustar[right];
			y = us - u;
			phase_star_[phase].pstar[right] = ps
					* std::pow(1 + (y * phase_star.g7_ / as), phase_star.g3_);
			G_s = ds
					* std::pow(phase_star_[phase].pstar[right] / ps,
							1 / phase_star.gama_);
			phase_star_[phase].dstar[right] = G_s;
		}
	}

}

void RiemannSolver::shockWave(const unsigned left, const unsigned right,
		const unsigned phase, bool Side, bool Pressure) {
	Phase_Riemann phase_star = phase_star_[phase];
	double A_s, B_s; //, a_s;
	double F_s;
	double ds, p, ps, us;
	double y, u;
	//Right is True
	if (Pressure) {
		if (Side) {
			//		d = phase_star.dstar[left];
			ds = phase_star.dstar[right];
			p = phase_star.pstar[left];
			ps = phase_star.pstar[right];
			//as = phase_star.getC(right);
			us = phase_star.ustar[right];
			A_s = phase_star.g5_ / ds;
			B_s = phase_star.g6_ * ps;
			F_s = (p - ps) * std::sqrt(A_s / (p + B_s));
			phase_star_[phase].dstar[left] = ds
					* ((phase_star.g8_ * ps + (phase_star.gama_ + 1) * p)
							/ (phase_star.g8_ * p + (phase_star.gama_ + 1) * ps));
			phase_star_[phase].ustar[left] = us + F_s;
		} else {
			ds = phase_star.dstar[left];
			//		d = phase_star.dstar[right];
			ps = phase_star.pstar[left];
			p = phase_star.pstar[right];
			//		as = phase_star.getC(left);
			us = phase_star.ustar[left];
			A_s = phase_star.g5_ / ds;
			B_s = phase_star.g6_ * ps;
			F_s = (p - ps) * std::sqrt(A_s / (p + B_s));
			phase_star_[phase].dstar[right] = ds
					* ((phase_star.g8_ * ps + (phase_star.gama_ + 1) * p)
							/ (phase_star.g8_ * p + (phase_star.gama_ + 1) * ps));
			phase_star_[phase].ustar[right] = us - F_s;
		}
	} else {
		p = 0;
//		p0 = p;
		//Right is True
		if (Side) {
			//		d = phase_star.dstar[left];
			ds = phase_star.dstar[right];
			ps = phase_star.pstar[right];
			//as = phase_star.getC(right);
			us = phase_star.ustar[right];
			A_s = phase_star.g5_ / ds;
			B_s = phase_star.g6_ * ps;
			u = phase_star_[phase].ustar[left];
			y = u - us;
			p = ps
					+ (std::pow(y, 2)
							+ std::sqrt(
									std::pow(y, 4)
											+ 4 * (ps + B_s) * A_s
													* std::pow(y, 2))) / 2
							/ A_s;
			phase_star_[phase].pstar[left] = p;
			phase_star_[phase].dstar[left] = ds
					* ((phase_star.g8_ * ps + (phase_star.gama_ + 1) * p)
							/ (phase_star.g8_ * p + (phase_star.gama_ + 1) * ps));
		} else {
			ds = phase_star.dstar[left];
			ps = phase_star.pstar[left];
			us = phase_star.ustar[left];
			A_s = phase_star.g5_ / ds;
			B_s = phase_star.g6_ * ps;
			u = phase_star_[phase].ustar[right];
			y = u - us;
			p = ps
					+ (std::pow(y, 2)
							+ std::sqrt(
									std::pow(y, 4)
											+ 4 * (ps + B_s) * A_s
													* std::pow(y, 2))) / 2
							/ A_s;
			phase_star_[phase].pstar[right] = p;
			phase_star_[phase].dstar[right] = ds
					* ((phase_star.g8_ * ps + (phase_star.gama_ + 1) * p)
							/ (phase_star.g8_ * p + (phase_star.gama_ + 1) * ps));
		}
	}

}

void RiemannSolver::SolidContact(const unsigned left_s, const unsigned right_s,
		const unsigned left_g, const unsigned right_g, bool Side) {
//	double v_1[2], p_1[2], v_2[2], alpha_1[2], alpha_2[2], rho;
//	double p1, p0;
//	double norm_err = 1000, gamma;
//	double left[2];
//	left[0] = left_s;
//	left[1] = left_g;
////	right[0] = right_s;
////	right[1] = right_g;
//	Eigen::Matrix3d Df;
//	gamma = phase_star_[0].gama_;
//	Eigen::Vector3d x, x0, f, dx;
//	for (unsigned i = 0; i < 2; ++i) {
//		v_1[i] = phase_star_[i].ustar[left[i]];
//		if (i == 1)
//			rho = phase_star_[i].dstar[left[i]];
//
//		if (i == 0)
//			v_2[i] = v_1[i];
//
//		p_1[i] = phase_star_[i].pstar[left[i]];
//		//p_2[i] = phase_star_[i].pstar[right];
//		if (!Side) {
//			alpha_1[i] = phase_star_[i].phi_l_;
//			alpha_2[i] = phase_star_[i].phi_r_;
//		} else {
//			alpha_1[i] = phase_star_[i].phi_r_;
//			alpha_2[i] = phase_star_[i].phi_l_;
//		}
//	}
//
//	x << phase_star_[1].ustar[right_g], phase_star_[1].pstar[right_g], phase_star_[0].pstar[right_g];
//	x0 = x;
//	while (std::abs(norm_err) > 0.001) {
//		p1 = x(1);
//		p0 = x(0);
//		f(0) = -(alpha_2[1] * std::pow((p1 / p_1[1]), (1 / gamma))
//				* (x(0) - v_1[0]) - alpha_1[1] * (v_1[1] - v_1[0]));
//		f(1) = -(alpha_2[0] * x(2) + alpha_2[1] * x(1) - alpha_1[0] * p_1[0]
//				- alpha_1[1] * p_1[1]
//				+ alpha_1[1] * rho * (v_1[1] - v_1[0]) * (x(0) - v_1[1]));
//		f(2) = -((gamma * x(1) / (gamma - 1) / rho)
//				* std::pow((p_1[1] / p1), (1 / gamma))
//				+ 0.5 * std::pow(p0 - v_2[0], 2)
//				- (gamma * p_1[1] / ((gamma - 1) * rho))
//				- 0.5 * std::pow(v_1[1] - v_1[0], 2));
//		Df(0, 0) = alpha_2[1] * std::pow(p1 / p_1[1], (1 / gamma));
//		Df(0, 1) = alpha_2[1] / gamma * std::pow(p1 / p_1[1], 1 / gamma) / x(1)
//				* (x(0) - v_1[0]);
//		Df(0, 2) = 0;
//		Df(1, 0) = alpha_1[1] * rho * (v_1[1] - v_1[0]);
//		Df(1, 1) = alpha_2[1];
//		Df(1, 2) = alpha_2[0];
//		Df(2, 0) = x(0) - v_1[0];
//		Df(2, 1) = std::pow(p_1[1] / p1, 1 / gamma) / rho;
//		Df(2, 2) = 0;
//		dx = Df.fullPivLu().solve(f);
////		if (dx(1) != dx(1))
////			dx = Df.partialPivLu().solve(f);
//		x = x0 + dx;
//		norm_err = x.norm() - x0.norm();
//		x0 = x;
//	}
//	std::cout << x << std::endl;
//	phase_star_[1].ustar[right_g] = x(0);
//	phase_star_[1].pstar[right_g] = x(1);
//	phase_star_[0].pstar[right_s] = x(2);
//	p1 = x(1);
//	phase_star_[1].dstar[right_g] = std::pow(
//			p1 / p_1[1] * std::pow(phase_star_[1].dstar[left_g], gamma),
//			1 / gamma);
}

void RiemannSolver::Configuration2(const unsigned left, const unsigned right) {
//	double ext_err = 1000, int_err = 1000;
//	double Df_u;
//	double u0;
//	phase_star_[1].pstar[1] = phase_star_[1].pstar[GAS_SIZE - 1];
//	phase_star_[0].ustar[1] = phase_star_[0].ustar[SOLID_SIZE - 1];
//	while (ext_err > 0.00001) {
//		if (phase_star_[1].pstar[1] > phase_star_[1].pstar[0])
//			shockWave(0, 1, 1, LEFT, true);
//		else
//			rarefactionWave(0, 1, 1, LEFT, true);
////		p0 = phase_star_[1].pstar[1];
//		u0 = phase_star_[0].ustar[1];
//		while (int_err < 0.00001) {
//			if (phase_star_[0].ustar[1] < phase_star_[0].ustar[0])
//				shockWave(0, 1, 0, LEFT, false);
//			else
//				rarefactionWave(0, 1, 0, LEFT, false);
////			SolidContact(1, 2, true);
//			phase_star_[0].ustar[3] = phase_star_[0].ustar[1];
//			if (phase_star_[0].ustar[3] > phase_star_[0].ustar[SOLID_SIZE - 1])
//				shockWave(3, 4, 0, RIGHT, false);
//			else
//				rarefactionWave(3, 4, 0, RIGHT, false);
//			Df_u = phase_star_[1].phi_l_ * phase_star_[1].dstar[1]
//					* (phase_star_[1].ustar[2] - phase_star_[1].ustar[1])
//					/ (phase_star_[0].phi_r_);
//			if (phase_star_[0].ustar[1] > phase_star_[0].ustar[SOLID_SIZE - 1])
//				Df_u += Df_Shock(3, 4, 0, true);
//			else
//				Df_u += Df_Rarefaction(3, 4, 0, true);
//
//			phase_star_[0].ustar[1] = u0
//					- (phase_star_[0].pstar[3] - phase_star_[0].pstar[2])
//							/ Df_u;
//		}
//	}
}

void RiemannSolver::Configuration1(const unsigned left, const unsigned right) {
//	double v_1[2], p_1[2], v_2[2], alpha_1[2], alpha_2[2], rho;
//	double norm_err = 1000, gamma;
//	Eigen::Matrix3d Df;
//	gamma = phase_star_[0].gama_;
//	Eigen::Vector3d x, x0, f, dx;
//	for (unsigned i = 0; i < 2; ++i) {
//		v_1[i] = phase_star_[i].ustar[left];
//		if (i == 1)
//			rho = phase_star_[i].dstar[right];
//
//		if (i == 0)
//			v_2[i] = v_1[i];
//
//		p_1[i] = phase_star_[i].pstar[left];
//		//p_2[i] = phase_star_[i].pstar[right];
//		alpha_1[i] = phase_star_[i].phi_l_;
//		alpha_2[i] = phase_star_[i].phi_r_;
//	}
//	x << phase_star_[1].ustar[right], phase_star_[1].pstar[right], phase_star_[0].pstar[right];
//	x0 = x;
//	while (std::abs(norm_err) > 0.001) {
//		double p1 = x(1);
//		double p0 = x(0);
//		f(0) = -(alpha_2[1] * (x(0) - v_1[0])
//				- alpha_1[1] * std::pow((p_1[1] / p1), (1 / gamma))
//						* (v_1[1] - v_1[0]));
//		f(1) = -(alpha_2[0] * x(2) + alpha_2[1] * x(1) - alpha_1[0] * p_1[0]
//				- alpha_1[1] * p_1[1]
//				+ alpha_2[1] * rho * (x(0) - v_1[0]) * (x(0) - v_1[1]));
//		f(2) = -((gamma * x(1) / (gamma - 1) / rho)
//				+ 0.5 * std::pow(p0 - v_2[0], 2)
//				- (gamma * p_1[1] / ((gamma - 1) * rho))
//						* std::pow((p1 / p_1[1]), (1 / gamma))
//				- 0.5 * std::pow(v_1[1] - v_1[0], 2));
//		Df(0, 0) = alpha_2[1];
//		Df(0, 1) = -alpha_1[1] / gamma * std::pow(p_1[1] / p1, 1 / gamma) / x(1)
//				* (v_1[1] - v_1[0]);
//		Df(0, 2) = 0;
//		Df(1, 0) = 2 * alpha_2[1] * rho * x(0)
//				- alpha_2[1] * rho * (v_1[0] + v_1[1]);
//		Df(1, 1) = alpha_2[1];
//		Df(1, 2) = alpha_2[0];
//		Df(2, 0) = x(0) - v_1[0];
//		Df(2, 1) = gamma / (gamma - 1) / rho
//				- p_1[1] / (gamma - 1) / rho * std::pow(p1 / p_1[1], 1 / gamma)
//						/ x(1);
//		Df(2, 2) = 0;
//		dx = Df.fullPivLu().solve(f);
//		x = x0 + dx;
//		norm_err = x.norm() - x0.norm();
//		x0 = x;
//	}
//	phase_star_[1].ustar[right] = x(0);
//	phase_star_[1].pstar[right] = x(1);
//	phase_star_[0].pstar[right] = x(2);
}

double RiemannSolver::Df_Shock(const unsigned left, const unsigned right,
		const unsigned phase, bool Side) {
//	Phase_Riemann phase_star = phase_star_[phase];
//	double A_s, B_s; //, a_s;
////	double F_s;
//	double ds, ps, us;
//	double y, u;
//	double Df;
////	p = 0;
////	p0 = p;
//	//Right is True
//	if (Side) {
//		//		d = phase_star.dstar[left];
//		ds = phase_star.dstar[right];
//		ps = phase_star.pstar[right];
//		//as = phase_star.getC(right);
//		us = phase_star.ustar[right];
//		A_s = g5_ / ds;
//		B_s = g6_ * ps;
//		u = phase_star_[phase].ustar[left];
//		y = u - us;
//	} else {
//		ds = phase_star.dstar[left];
//		ps = phase_star.pstar[left];
//		us = phase_star.ustar[left];
//		A_s = g5_ / ds;
//		B_s = g6_ * ps;
//		u = phase_star_[phase].ustar[right];
//		y = u - us;
//	}
//	Df = (2 * y
//			+ (4 * std::pow(y, 3) + 8 * (ps + B_s) * A_s * y)
//					/ (2
//							* std::sqrt(
//									std::pow(y, 4)
//											+ 4 * (ps + B_s) * A_s
//													* std::pow(y, 2))))
//			/ (2 * A_s);
//	return Df;
}

double RiemannSolver::Df_Rarefaction(const unsigned left, const unsigned right,
		const unsigned phase, bool Side) {
//	Phase_Riemann phase_star = phase_star_[phase];
//	//	double a_s;
////	double F_s;
//	double ps, as, us;
//	double u, y;
//	double Df;
//	if (Side) {
//		//d = phase_star.dstar[left];
////		ds = phase_star.dstar[right];
//		ps = phase_star.pstar[right];
//		as = phase_star.getC(right);
//		us = phase_star.ustar[right];
//		u = phase_star.ustar[left];
//		y = u - us;
//		Df = phase_star_[0].gama_ * ps / as * std::pow(1 + y / as * g7_, g2_);
//	} else {
////		ds = phase_star.dstar[left];
//		//d = phase_star.dstar[right];
//		ps = phase_star.pstar[left];
////		p = phase_star.pstar[right];
//		u = phase_star.ustar[right];
//		as = phase_star.getC(left);
//		us = phase_star.ustar[left];
//		y = us - u;
//		Df = -phase_star_[0].gama_ * ps / as * std::pow(1 + y / as * g7_, g2_);
//	}
//	return Df;
}

void RiemannSolver::InitializeU0(const double& u0, double& u_p0) {
//	phase_star_[0].ustar[2] = u0;
//	phase_star_[0].ustar[1] = u0;
//	if (phase_star_[0].ustar[2] > phase_star_[0].ustar[3])
//		shockWave(2, 3, 0, RIGHT, false);
//	else
//		rarefactionWave(2, 3, 0, RIGHT, false);
//
//	phase_star_[1].pstar[2] = phase_star_[1].pstar[3] - 0.1;
//	phase_star_[1].ustar[2] = phase_star_[1].ustar[3] - 0.1;
//	phase_star_[0].pstar[1] = phase_star_[0].pstar[2] - 0.1;
//
//	//left_s, right_s, left_g, right_g
//	std::cout << "ds: " << phase_star_[0].dstar[2] << std::endl;
//	SolidContact(2, 1, 3, 2, RIGHT);
//	phase_star_[0].ustar[4] = phase_star_[0].ustar[2];
//	if (phase_star_[0].ustar[4] < phase_star_[0].ustar[0])
//		shockWave(0, 4, 0, LEFT, false);
//	else
//		rarefactionWave(0, 4, 0, LEFT, false);
//
//	u_p0 = phase_star_[0].pstar[4] - phase_star_[0].pstar[1];
//	phase_star_[0].ustar[2] = 0;
}

void RiemannSolver::InitializeP0(const double& p0, double& p_u0) {
//	double u0, u_p0, int_err = 1000, u1;
//	u0 = 0.0;
//	u1 = 0.01;
//	phase_star_[1].pstar[GAS_SIZE - 3] = p0;
//
//	if (phase_star_[1].pstar[GAS_SIZE - 3] > phase_star_[1].pstar[GAS_SIZE - 2])
//		shockWave(3, 4, 1, RIGHT, true);
//	else
//		rarefactionWave(3, 4, 1, RIGHT, true);
//
//	InitializeU0(u0, u_p0);
//
//	phase_star_[0].ustar[2] = u1;
//
//	if (phase_star_[0].ustar[2] > phase_star_[0].ustar[3])
//		shockWave(2, 3, 0, RIGHT, false);
//	else
//		rarefactionWave(2, 3, 0, RIGHT, false);
//
//	phase_star_[0].pstar[1] = phase_star_[0].pstar[2] - 0.1;
//
//	while (int_err > 0.00001) {
//		phase_star_[0].ustar[1] = phase_star_[0].ustar[2];
//
//		std::cout << "ur_: " << phase_star_[0].ustar[3] << std::endl;
//
//		if (phase_star_[0].ustar[2] > phase_star_[0].ustar[3])
//			shockWave(2, 3, 0, RIGHT, false);
//		else
//			rarefactionWave(2, 3, 0, RIGHT, false);
//		//left_s, right_s, left_g, right_g
//		std::cout << "ds: " << phase_star_[0].dstar[2] << std::endl;
//		SolidContact(2, 1, 3, 2, RIGHT);
//		phase_star_[0].ustar[4] = phase_star_[0].ustar[2];
//		if (phase_star_[0].ustar[4] < phase_star_[0].ustar[0])
//			shockWave(0, 4, 0, LEFT, false);
//		else
//			rarefactionWave(0, 4, 0, LEFT, false);
//		phase_star_[0].ustar[2] = u1
//				- (phase_star_[0].pstar[4] - phase_star_[0].pstar[1])
//						/ (phase_star_[0].pstar[4] - phase_star_[0].pstar[1]
//								- u_p0) * (u1 - u0); //(phase_star_[0].pstar[4] - phase_star_[0].pstar[1])
//				//	/ Df_u;
//		std::cout << "u_s2: " << phase_star_[0].ustar[2] << std::endl;
//
//		u_p0 = phase_star_[0].pstar[4] - phase_star_[0].pstar[1];
//		int_err = std::abs(u_p0); //std::abs((phase_star_[0].ustar[2] - u0) / u0);
//		std::cout << "Error: " << int_err << std::endl;
//		u0 = u1;
//		u1 = phase_star_[0].ustar[2];
//	}
//	phase_star_[1].pstar[1] = phase_star_[1].pstar[2];
//	if (phase_star_[1].pstar[1] < phase_star_[1].pstar[0])
//		shockWave(0, 1, 1, LEFT, true);
//	else
//		rarefactionWave(0, 1, 1, LEFT, true);
//	p_u0 = phase_star_[1].ustar[1] - phase_star_[1].ustar[2];
}

void RiemannSolver::Configuration2dash() {
//	double ext_err = 1000, p0, u0, u_p0, int_err = 1000, u1, p1, p_u0;
//	p0 = 0.95;
//	p1 = 0.8707;
//	u0 = 0;
//	u1 = 0.01;
//	InitializeP0(p0, p_u0);
//	phase_star_[1].pstar[GAS_SIZE - 3] = p1;
//	phase_star_[1].pstar[2] = 0.5;
//	phase_star_[1].ustar[2] = 0;
//	phase_star_[0].pstar[1] = p1 - 0.1;
//	while (ext_err > 0.00001) {
//		std::cout << "Exterior Iteration" << std::endl;
//		if (phase_star_[1].pstar[GAS_SIZE - 3]
//				> phase_star_[1].pstar[GAS_SIZE - 2])
//			shockWave(3, 4, 1, RIGHT, true);
//		else
//			rarefactionWave(3, 4, 1, RIGHT, true);
////		u0 = 0.01;
////		u1 = 0;
//		InitializeU0(u0, u_p0);
//		phase_star_[0].ustar[2] = u1;
//		int_err = 1000;
//		while (int_err > 0.0001) {
//			phase_star_[0].ustar[1] = phase_star_[0].ustar[2];
//
//			std::cout << "ur_: " << phase_star_[0].ustar[3] << std::endl;
//
//			if (phase_star_[0].ustar[2] > phase_star_[0].ustar[3])
//				shockWave(2, 3, 0, RIGHT, false);
//			else
//				rarefactionWave(2, 3, 0, RIGHT, false);
//			//left_s, right_s, left_g, right_g
//			std::cout << "ds: " << phase_star_[0].dstar[2] << std::endl;
//			SolidContact(2, 1, 3, 2, RIGHT);
//			phase_star_[0].ustar[4] = phase_star_[0].ustar[2];
//			if (phase_star_[0].ustar[4] < phase_star_[0].ustar[0])
//				shockWave(0, 4, 0, LEFT, false);
//			else
//				rarefactionWave(0, 4, 0, LEFT, false);
//			phase_star_[0].ustar[2] = u1
//					- (phase_star_[0].pstar[4] - phase_star_[0].pstar[1])
//							/ (phase_star_[0].pstar[4] - phase_star_[0].pstar[1]
//									- u_p0) * (u1 - u0); //(phase_star_[0].pstar[4] - phase_star_[0].pstar[1])
//					//	/ Df_u;
//			std::cout << "u_s2: " << phase_star_[0].ustar[2] << std::endl;
//
//			u_p0 = phase_star_[0].pstar[4] - phase_star_[0].pstar[1];
//			int_err = std::abs(u_p0); //std::abs((phase_star_[0].ustar[2] - u0) / u0);
//			std::cout << "Error: " << int_err << std::endl;
//			u0 = u1;
//			u1 = phase_star_[0].ustar[2];
//		}
//		phase_star_[1].pstar[1] = phase_star_[1].pstar[2];
//		if (phase_star_[1].pstar[1] < phase_star_[1].pstar[0])
//			shockWave(0, 1, 1, LEFT, true);
//		else
//			rarefactionWave(0, 1, 1, LEFT, true);
//		phase_star_[1].pstar[GAS_SIZE - 3] = p1
//				- (phase_star_[1].ustar[1] - phase_star_[1].ustar[2])
//						/ (phase_star_[1].ustar[1] - phase_star_[1].ustar[2]
//								- p_u0) * (p1 - p0);
//		std::cout << "pg: " << phase_star_[1].pstar[GAS_SIZE - 3] << std::endl;
//		p_u0 = phase_star_[1].ustar[1] - phase_star_[1].ustar[2];
//		ext_err = std::abs(p_u0);
//		std::cout << "Ext. error: " << ext_err << std::endl;
//		p0 = p1;
//		p1 = phase_star_[1].pstar[GAS_SIZE - 3];
//	}
}

std::vector<double> RiemannSolver::consRStarVar() {

//	std::vector<double> consStarVar(3);
//	consStarVar[0] = dstar_r;
//	consStarVar[1] = dstar_r * ustar;
//	consStarVar[2] = dstar_r
//			* (0.5 * std::pow(ustar, 2.0) + pstar / (gama_r_ - 1) / dstar_r);
//	return consStarVar;
}

void RiemannSolver::SpecialConfiguration() {
    //if (phase_star_[0].phi_l_==1){
    //    phase_star_[1].ul_ = phase_star_[0].ul_;
    //    phase_star_[1].pl_ = phase_star_[0].pl_;
    //    phase_star_[1].cl_ = phase_star_[0].cl_;
    //} else {
    //    phase_star_[0].ul_ = phase_star_[1].ul_;
    //    phase_star_[0].pl_ = phase_star_[1].pl_;
    //    phase_star_[0].cl_ = phase_star_[1].cl_;
    //}

    //if (phase_star_[0].phi_r_==1){
    //    phase_star_[1].ur_ = phase_star_[0].ur_;
    //    phase_star_[1].pr_ = phase_star_[0].pr_;
    //    phase_star_[1].cr_ = phase_star_[0].cr_;
    //} else {
    //    phase_star_[0].ur_ = phase_star_[1].ur_;
    //    phase_star_[0].pr_ = phase_star_[1].pr_;
    //    phase_star_[0].cr_ = phase_star_[1].cr_;
    //}
    double u, p, d;
    solve();
    std::cout << phase_star_[0].pstar[0] << "\t" << phase_star_[0].ustar[0] << std::endl;
    std::cout << phase_star_[1].pstar[0] << "\t" << phase_star_[1].ustar[0] << std::endl;
    std::cout << r_ << '\t' 
              << phase_star_[r_].dr_
               *(phase_star_[r_].pstar[0]
                /phase_star_[r_].pr_
                +phase_star_[r_].g6_
                +phase_star_[r_].gama_
               /(phase_star_[r_].gama_+1)
                *phase_star_[r_].pinf_
                /phase_star_[r_].pr_) / 
                (phase_star_[r_].pstar[0]/
                 phase_star_[r_].pr_
                *phase_star_[r_].g6_+1.0
                +phase_star_[r_].gama_/
                (phase_star_[r_].gama_+1)
                *phase_star_[r_].pinf_
                /phase_star_[r_].pr_) << std::endl;
    phase_star_[l_].Sample(0,u,p,d,phase_star_[l_].pstar[0],phase_star_[l_].ustar[0]);
    std::cout << "Density: " << d << std::endl;
 

//	double ext_err = 1000, int_err = 1000;
//	double Df_u;
//	double u0;
//	phase_star_[1].pstar[1] = phase_star_[1].pstar[GAS_SIZE - 1];
//	phase_star_[0].ustar[1] = phase_star_[0].ustar[SOLID_SIZE - 1];
//	while (ext_err > 0.00001) {
//		if (phase_star_[1].pstar[1] > phase_star_[1].pstar[0])
//			shockWave(0, 1, 1, LEFT, true);
//		else
//			rarefactionWave(0, 1, 1, LEFT, true);
////		p0 = phase_star_[1].pstar[1];
//		u0 = phase_star_[0].ustar[1];
//		while (int_err < 0.00001) {
//			if (phase_star_[0].ustar[1] < phase_star_[0].ustar[0])
//				shockWave(0, 1, 0, LEFT, false);
//			else
//				rarefactionWave(0, 1, 0, LEFT, false);
////			SolidContact(1, 2, true);
//			phase_star_[0].ustar[3] = phase_star_[0].ustar[1];
//			if (phase_star_[0].ustar[3] > phase_star_[0].ustar[SOLID_SIZE - 1])
//				shockWave(3, 4, 0, RIGHT, false);
//			else
//				rarefactionWave(3, 4, 0, RIGHT, false);
//			Df_u = phase_star_[1].phi_l_ * phase_star_[1].dstar[1]
//					* (phase_star_[1].ustar[2] - phase_star_[1].ustar[1])
//					/ (phase_star_[0].phi_r_);
//			if (phase_star_[0].ustar[1] > phase_star_[0].ustar[SOLID_SIZE - 1])
//				Df_u += Df_Shock(3, 4, 0, true);
//			else
//				Df_u += Df_Rarefaction(3, 4, 0, true);
//
//			phase_star_[0].ustar[1] = u0
//					- (phase_star_[0].pstar[3] - phase_star_[0].pstar[2])
//							/ Df_u;
//		}
//	}
}
/* namespace CSE */

}

