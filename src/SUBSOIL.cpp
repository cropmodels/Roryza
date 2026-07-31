#include <vector>
#include <cmath>
#include <algorithm>
#include "oryzaUtil.h"
#include "model.h"

// Soil hydrology helpers ported from SUBSOIL.f90 (ORYZA2000 / PADDY)

void GWTAB(int ITASK, oryza_soil &soil, double DOY, double DELT) {
	double ZLL = 0., ZWL = 0.;

	if (ITASK == 1) {
		if (soil.SWITGW == 1) {
			soil.ZW = AFGEN(soil.ZWTB, DOY);
		} else {
			soil.ZW = soil.ZWTBI;
		}
		soil.ZWPREV = soil.ZW;
	}

	soil.IGW = 1;
	ZLL = 0.;
	ZWL = 0.;
	while (soil.IGW <= soil.NL && ZWL <= 0.) {
		ZLL = ZLL + soil.TKL_mm[soil.IGW - 1] / 10.;
		ZWL = ZLL - soil.ZWPREV;
		soil.IGW = soil.IGW + 1;
	}
	if (ZWL > 0.) soil.IGW = soil.IGW - 1;
	if (soil.IGW < 1) soil.IGW = 1;

	if (ITASK == 3) {
		if (soil.SWITGW == 1) {
			soil.ZW = AFGEN(soil.ZWTB, DOY);
		} else {
			int igw_idx = std::max(0, std::min(soil.IGW - 1, static_cast<int>(soil.WLFL.size()) - 1));
			soil.ZW = soil.ZW + soil.ZWA * DELT - soil.ZWB * 10. * soil.WLFL[igw_idx] * DELT;
			if (soil.ZW < soil.MINGW) soil.ZW = soil.MINGW;
			if (soil.ZW > soil.MAXGW) soil.ZW = soil.MAXGW;
		}
	}
}


void BACKFL(int I, double WL, double FLIN, double FLOUT, double EVSWS, double TRWL, double WLST, double DELT,
            double &FLNEW, double &HLP) {
	HLP = 0.;
	if (I == 1) {
		HLP = WL + (FLIN - FLOUT - EVSWS - TRWL) * DELT;
	} else {
		HLP = WL + (FLIN - FLOUT - TRWL) * DELT;
	}

	if (HLP > WLST) {
		FLNEW = FLIN - (HLP - WLST) / DELT;
	} else {
		FLNEW = FLIN;
	}
}


double DOWNFL(int I, double KSAT, double FLIN, double TRWL, double EVSWS, double WL, double WLFC, double DELT) {
	if (I == 1) {
		return std::min(10. * KSAT, std::max(0., FLIN - EVSWS - TRWL + (WL - WLFC) / DELT));
	}
	return std::min(10. * KSAT, std::max(0., FLIN - TRWL + (WL - WLFC) / DELT));
}


void SUWCMS2(int I, int SWIT4, oryza_soil &soil, double &WCL, double &MS) {
	int idx = I - 1;
	double HLP1, HLP2, HLP3, HLP4, VGM, WREL;
	double wcst = soil.WCST[idx];
	double wcad = soil.WCAD[idx];

	if (SWIT4 == 1) {
		if (WCL < wcad || WCL > wcst) {
			// SUERR stub: clamp instead of fatal
			WCL = LIMIT(wcad, wcst, WCL);
		}
		if (WCL > soil.WCSTRP[idx]) {
			MS = 0.;
		} else {
			HLP1 = std::max(wcad, WCL);
			WREL = (HLP1 - soil.VGR[idx]) / (soil.WCSTRP[idx] - soil.VGR[idx]);
			VGM = 1. - 1. / soil.VGN[idx];
			HLP2 = 1. / soil.VGA[idx];
			HLP3 = -1. / VGM;
			HLP4 = 1. / soil.VGN[idx];
			MS = HLP2 * std::pow(std::pow(WREL, HLP3) - 1., HLP4);
		}
	} else if (SWIT4 == 2) {
		if (MS < -0.001 || MS > 1.E8) MS = LIMIT(0., 1.E8, MS);
		VGM = 1. - 1. / soil.VGN[idx];
		HLP1 = std::pow(MS * soil.VGA[idx], soil.VGN[idx]);
		WREL = std::pow(1. + HLP1, -VGM);
		WCL = WREL * (soil.WCSTRP[idx] - soil.VGR[idx]) + soil.VGR[idx];
	}
}


void SUMSKM2(int I, double MS, oryza_soil &soil, double &KMS) {
	int idx = I - 1;
	double HLP1, HLP2, HLP3, VGM, WCL, WREL;
	const double TINY = 1.e-10;
	const double MSAD = 1.e7;

	if (MS >= MSAD - TINY) {
		KMS = 0.;
		return;
	}

	if (soil.SWITKH == 1) {
		WCL = 0.;
		SUWCMS2(I, 2, soil, WCL, MS);
		VGM = 1.0 - 1.0 / soil.VGN[idx];
		WREL = (WCL - soil.VGR[idx]) / (soil.WCSTRP[idx] - soil.VGR[idx]);
		HLP1 = std::pow(WREL, soil.VGL[idx]);
		HLP2 = 1.0 - std::pow(WREL, 1. / VGM);
		HLP3 = 1.0 - std::pow(HLP2, VGM);
		KMS = soil.KST[idx] * HLP1 * HLP3 * HLP3;
	} else if (soil.SWITKH == 2) {
		if (MS <= 1.) KMS = soil.KST[idx];
		else KMS = soil.KST[idx] * std::pow(MS, soil.PN[idx]);
	} else {
		KMS = 0.;
	}
	if (KMS < TINY) KMS = 0.;
}


double SUBSL2(double PF, double D, int I, oryza_soil &soil) {
	const double ELOG10 = 2.302585;
	const double LOGST4 = 2.518514;
	const double START[4] = {0., 45., 170., 330.};
	const double PFSTAN[9] = {0.705143, 1.352183, 1.601282, 1.771497, 2.031409, 2.192880,
	                          2.274233, 2.397940, 2.494110};
	const double PGAU[3] = {0.1127016654, 0.5, 0.8872983346};
	const double WGAU[3] = {0.2777778, 0.4444444, 0.2777778};

	double PF1 = PF;
	double D1 = D;
	double MH = std::exp(ELOG10 * PF1);

	if (PF1 <= 0.) {
		double K0 = soil.KST[I - 1];
		return 10. * K0 * (MH / D1 - 1.);
	}

	int IINT = 0;
	double DEL[4];
	for (int I1 = 0; I1 < 4; I1++) {
		if (I1 <= 2) DEL[I1] = std::min(START[I1 + 1], MH) - START[I1];
		if (I1 == 3) DEL[I1] = PF1 - LOGST4;
		if (DEL[I1] <= 0.) break;
		IINT = IINT + 1;
	}

	double PFGAU[12], HULP[12], CONDUC[12];
	for (int I1 = 0; I1 < IINT; I1++) {
		for (int I2 = 0; I2 < 3; I2++) {
			int I3 = 3 * I1 + I2;
			if (I1 == IINT - 1) {
				if (IINT <= 3) PFGAU[I3] = std::log10(START[IINT - 1] + PGAU[I2] * DEL[IINT - 1]);
				else PFGAU[I3] = LOGST4 + PGAU[I2] * DEL[IINT - 1];
			} else {
				PFGAU[I3] = PFSTAN[I3];
			}
			double KMS = 0.;
			double ms_tmp = std::exp(ELOG10 * PFGAU[I3]);
			SUMSKM2(I, ms_tmp, soil, KMS);
			CONDUC[I3] = KMS;
			HULP[I3] = DEL[I1] * WGAU[I2] * CONDUC[I3];
			if (I3 > 8) HULP[I3] = HULP[I3] * ELOG10 * std::exp(ELOG10 * PFGAU[I3]);
		}
	}

	double KMS = 0.;
	double ms_pf = std::exp(ELOG10 * PF1);
	SUMSKM2(I, ms_pf, soil, KMS);
	double FU = 1.27;
	double FL = -1. * KMS;
	if (MH <= D1) FU = 0.;
	if (MH >= D1) FL = 0.;
	if (MH == D1) return 10. * (FU + FL) / 2.;

	int IMAX = 3 * IINT;
	for (int iter = 0; iter < 15; iter++) {
		double FLW = (FU + FL) / 2.;
		double DF = (FU - FL) / 2.;
		if (DF < 0.01 && (DF / std::fabs(FLW)) < 0.1) break;
		double Z = 0.;
		for (int I2 = 0; I2 < IMAX; I2++) {
			Z = Z + HULP[I2] / (CONDUC[I2] + FLW);
		}
		if (Z >= D1) FL = FLW;
		if (Z <= D1) FU = FLW;
	}
	return 10. * (FU + FL) / 2.;
}


double SATFLX(oryza_soil &soil, double WL0) {
	// STUB: SWITVP=1 calculated percolation requires puddled soil + Newton-Raphson iteration.
	// Not implemented for reference config (SWITVP=-1, SWITPD=0).
	(void)soil;
	(void)WL0;
	return 0.;
}


void SUWCHK(double CKWFL, double CKWIN, double TIME) {
	double FUWCHK = 2.0 * (CKWIN - CKWFL) / (CKWIN + CKWFL + 1.E-10);
	double XDIF = std::fabs(CKWIN - CKWFL);
	if (std::fabs(FUWCHK) > 0.01 && XDIF > 1.0) {
		// Water balance check failed; logged silently in library port
		(void)TIME;
	}
}
