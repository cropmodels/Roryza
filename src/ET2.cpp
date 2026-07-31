//----------------------------------------------------------------------//
//  SUBROUTINE ET2                                                      //
//  Potential soil evaporation and crop transpiration (Penman/Makkink/PT)
//----------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <string>
#include "model.h"

namespace {
	double ETDCUM1 = 0, EVSCCUM1 = 0, TRCCUM1 = 0;
	double ETDCUM2 = 0, EVSCCUM2 = 0, TRCCUM2 = 0;
	double ETDCUM3 = 0, EVSCCUM3 = 0, TRCCUM3 = 0;
}

void ET2_initialize() {
	ETDCUM1 = EVSCCUM1 = TRCCUM1 = 0.;
	ETDCUM2 = EVSCCUM2 = TRCCUM2 = 0.;
	ETDCUM3 = EVSCCUM3 = TRCCUM3 = 0.;
}


void ET2_rate(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
              double FAOF, double WL0, const std::vector<double> &WCLQT, const std::vector<double> &WCST, double LAI,
              double &EVSC, double &ETD, double &TRC) {

	int ISURF;
	double ALB, DT, ETAE = 0, ETRD = 0, RF, RFS;

	if (WL0 > 5.) {
		ALB = 0.05;
		RFS = ALB;
	} else {
		ALB = 0.25;
		double wc = WCLQT.empty() ? 0.3 : WCLQT[0];
		double ws = WCST.empty() ? 0.3 : WCST[0];
		RFS = ALB * (1. - 0.5 * wc / ws);
	}

	RF = RFS * exp(-0.5 * LAI) + 0.25 * (1. - exp(-0.5 * LAI));

	if (ETMOD == "PENMAN") {
		if (CROPSTA < 3) {
			ISURF = (WL0 > 5.) ? 1 : 2;
		} else {
			ISURF = 3;
		}
		std::vector<double> setpmd = SETPMD(IDOY, LAT, ISURF, RF, ANGA, ANGB, 0., RDD, TMDA, WN, VP);
		ETD = setpmd[0];
		ETRD = setpmd[1];
		ETAE = setpmd[2];
		DT = setpmd[3];
		(void)DT;
	} else if (ETMOD == "MAKKINK") {
		ETD = SETMKD(RDD, TMDA);
		ETRD = 0.75 * ETD;
		ETAE = ETD - ETRD;
	} else if (ETMOD == "PRIESTLEY TAYLOR") {
		ETD = SETPTD(IDOY, LAT, RF, RDD, TMDA);
		ETRD = 0.75 * ETD;
		ETAE = ETD - ETRD;
	} else {
		ETD = 0.;
		ETRD = 0.;
		ETAE = 0.;
	}

	ETD  = ETD  * FAOF;
	ETRD = ETRD * FAOF;
	ETAE = ETAE * FAOF;

	EVSC = exp(-0.5 * LAI) * (ETRD + ETAE);
	EVSC = std::max(EVSC, 0.);
	// Bas, June 2006: transpiration also before transplanting
	TRC = ETRD * (1. - exp(-0.5 * LAI)) + ETAE * std::min(2.0, LAI);
}


void ET2_state(double DELT, int CROPSTA, const std::string &ESTAB, double ETD, double EVSC, double TRC) {
	ETDCUM1  += ETD * DELT;
	EVSCCUM1 += EVSC * DELT;
	TRCCUM1  += TRC * DELT;
	if (CROPSTA >= 1) {
		ETDCUM2  += ETD * DELT;
		EVSCCUM2 += EVSC * DELT;
		TRCCUM2  += TRC * DELT;
	}
	if (ESTAB == "TRANSPLANT" && CROPSTA >= 3) {
		ETDCUM3  += ETD * DELT;
		EVSCCUM3 += EVSC * DELT;
		TRCCUM3  += TRC * DELT;
	}
}
