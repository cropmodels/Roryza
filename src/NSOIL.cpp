//----------------------------------------------------------------------//
//  SUBROUTINE NSOIL                                                    //
//  Used in ORYZA2000 model version 1.0                                 //
//  Date   : December 2001                                              //
//  Author : B.A.M. Bouman                                              //
//  Purpose: This module calculates the supply of N from the soil.      //
//----------------------------------------------------------------------//

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include "model.h"
#include "oryzaUtil.h"

// Expand compact FERTIL events (day, amount, ...) into a pulsed AFGEN table
// so LINT2/AFGEN returns the amount on the event day and 0 elsewhere.
// If the table already contains zero amounts (FORTRAN padded form), keep it.
static std::vector<double> normalize_fertil(std::vector<double> fert) {
	if (fert.size() < 2) {
		return {0., 0., 366., 0.};
	}
	if (fert.size() % 2 != 0) {
		fert.pop_back();
	}

	bool has_zero_amt = false;
	for (size_t i = 1; i < fert.size(); i += 2) {
		if (fert[i] == 0.0) {
			has_zero_amt = true;
			break;
		}
	}
	if (has_zero_amt) {
		return fert;
	}

	std::vector<std::pair<double, double>> pts;
	pts.emplace_back(0., 0.);
	for (size_t i = 0; i + 1 < fert.size(); i += 2) {
		double d = fert[i];
		double a = fert[i + 1];
		if (d > 0.) {
			pts.emplace_back(d - 1., 0.);
		}
		pts.emplace_back(d, a);
		pts.emplace_back(d + 1., 0.);
	}
	pts.emplace_back(366., 0.);
	std::sort(pts.begin(), pts.end(),
		[](const std::pair<double, double> &a, const std::pair<double, double> &b) {
			return a.first < b.first;
		});

	std::vector<double> out;
	out.reserve(pts.size() * 2);
	for (size_t i = 0; i < pts.size(); ++i) {
		// Keep later duplicate X (event day over preceding zero at same X if any)
		if (i + 1 < pts.size() && pts[i].first == pts[i + 1].first) {
			continue;
		}
		out.push_back(pts[i].first);
		out.push_back(pts[i].second);
	}
	return out;
}


void nsoil_initialize(oryza_model &m) {
	m.FERTIL_TB = normalize_fertil(m.control.FERTIL);
	if (m.control.RECNIT.size() < 2) {
		m.control.RECNIT = {0., 0.30, 0.2, 0.35, 0.4, 0.50, 0.8, 0.75, 1.0, 0.75, 2.5, 0.75};
	}
	m.crop.TNSOIL = 0.;
	m.NFERTP = 0.;
	m.XFERT = 0.;
	m.crop.NACR = 0.;
}


void nsoil_rate(oryza_model &m) {
	double fert = 0.;
	if (m.FERTIL_TB.size() >= 2) {
		fert = AFGEN(m.FERTIL_TB, m.crop.DAE);
	}
	double recov = 0.;
	if (m.control.RECNIT.size() >= 2) {
		recov = AFGEN(m.control.RECNIT, m.crop.DVS);
	}
	m.XFERT = fert * recov;
}


void nsoil_state(oryza_model &m) {
	double soilsp = m.control.SOILSP;
	double nacr = m.crop.NACR;
	// INTGRL(NFERTP, XFERT - MAX(0, NACR - SOILSP), DELT)
	m.NFERTP = m.NFERTP + (m.XFERT - std::max(0., nacr - soilsp)) * m.DELT;
	m.crop.TNSOIL = m.NFERTP + soilsp;
}


// Legacy ITASK wrapper (unused by driver; kept for header compatibility)
double NSOIL(int ITASK, int /*IUNITD*/, int /*IUNITL*/, std::string /*FILEIT*/,
	double /*OUTPUT*/, double DELT, double DAE, double DVS, double NACR) {
	(void)DELT; (void)DAE; (void)DVS; (void)NACR;
	if (ITASK == 1) return 0.;
	return 0.;
}
