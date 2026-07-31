//----------------------------------------------------------------------//
//  SUBROUTINE WSTRESS                                                  //
//----------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "oryzaUtil.h"
#include "model.h"

namespace {
	double g_ULLS = 74.13, g_LLLS = 794.33;
	double g_ULDL = 630.95, g_LLDL = 1584.89;
	double g_ULLE = 1.45, g_LLLE = 1404.;
	double g_ULRT = 74.13, g_LLRT = 1584.89;
	std::string g_SWIRTR = "DATA";
	double g_SWIRTRF = 0.003297;
	int g_IESTAB = 1;
	const double TINY = 0.0000000001;
}

void WSTRESS_initialization(double DELT, double TRC, double ZRT, const std::vector<double> &TKL, int NL, int CROPSTA,
                            const std::vector<double> &WCLQT, const std::vector<double> &WCWP, const std::vector<double> &WCAD,
                            const std::vector<double> &MSKPA, const oryza_crop &crop, const std::string &ESTAB,
                            double &TRW, std::vector<double> &TRWL, double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW) {
	(void)DELT; (void)TRC; (void)ZRT; (void)TKL; (void)NL; (void)CROPSTA;
	(void)WCLQT; (void)WCWP; (void)WCAD; (void)MSKPA;

	g_ULLS = crop.ULLS;
	g_LLLS = crop.LLLS;
	g_ULDL = crop.ULDL;
	g_LLDL = crop.LLDL;
	g_ULLE = crop.ULLE;
	g_LLLE = crop.LLLE;
	g_ULRT = crop.ULRT;
	g_LLRT = crop.LLRT;
	g_SWIRTR = crop.SWIRTR;
	g_SWIRTRF = crop.SWIRTRF;
	g_IESTAB = (ESTAB == "TRANSPLANT") ? 2 : 1;

	LESTRS = 1.;
	PCEW = 1.;
	CPEW = 1.;
	LRSTRS = 1.;
	LDSTRS = 1.;
	for (int i = 0; i < NL; i++) TRWL[i] = 0.;
	TRW = 0.;
}


void WSTRESS_rate(double DELT, double TRC, double ZRT, const std::vector<double> &TKL, int NL, int CROPSTA,
                  const std::vector<double> &WCLQT, const std::vector<double> &WCWP, const std::vector<double> &WCAD,
                  const std::vector<double> &MSKPA, double &TRW, std::vector<double> &TRWL,
                  double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW) {

	double TRRM, ZRTL, ZLL, LRAV, LEAV, LDAV;
	TRWL.assign(10, 0.);
	std::vector<double> LR(10, 0.), LE(10, 0.), LD(10, 0.), TRR(10, 0.), WLA(10, 0.);

	if (CROPSTA > g_IESTAB) {
		TRRM = TRC / (ZRT + 1.0E-10);
		TRW = 0.;
		ZLL = 0.;
		LRAV = 0.;
		LEAV = 0.;
		LDAV = 0.;

		for (int i = 0; i < NL; i++) {
			ZRTL = std::min(TKL[i], std::max(ZRT - ZLL, 0.0));

			LR[i] = (std::log10(MSKPA[i] + TINY) - std::log10(g_LLLS)) / (std::log10(g_ULLS) - std::log10(g_LLLS));
			LR[i] = LIMIT(0., 1., LR[i]);
			LRAV = LRAV + (ZRTL / (ZRT + TINY)) * LR[i];

			LE[i] = (std::log10(MSKPA[i] + TINY) - std::log10(g_LLLE)) / (std::log10(g_ULLE) - std::log10(g_LLLE));
			LE[i] = LIMIT(0., 1., LE[i]);
			LEAV = LEAV + (ZRTL / (ZRT + TINY)) * LE[i];

			LD[i] = (std::log10(MSKPA[i] + TINY) - std::log10(g_LLDL)) / (std::log10(g_ULDL) - std::log10(g_LLDL));
			LD[i] = LIMIT(0., 1., LD[i]);
			LDAV = LDAV + (ZRTL / (ZRT + TINY)) * LD[i];

			if (MSKPA[i] >= 10000.) {
				TRR[i] = 0.;
			} else {
				if (g_SWIRTR == "DATA") {
					TRR[i] = (std::log10(MSKPA[i] + TINY) - std::log10(g_LLRT)) / (std::log10(g_ULRT) - std::log10(g_LLRT));
					TRR[i] = LIMIT(0., 1., TRR[i]);
				} else {
					TRR[i] = 2. / (1. + std::exp(g_SWIRTRF * MSKPA[i]));
				}
			}
			TRR[i] = LIMIT(0., 1., TRR[i]);
			WLA[i] = std::max(0.0, (WCLQT[i] - WCWP[i]) * ZRTL * 1000.);
			TRWL[i] = std::min(TRR[i] * ZRTL * TRRM, WLA[i] / DELT);
			TRW = TRW + TRWL[i];
			ZLL = ZLL + TKL[i];
		}

		for (int i = 0; i < NL; i++) {
			if (TRW < TRC) {
				if (TRR[i] >= 1 && TRWL[i] < WLA[i] / DELT) {
					TRWL[i] = std::min(WLA[i] / DELT, TRWL[i] + (TRC - TRW));
				}
				TRW = 0.;
				for (int j = 0; j < NL; j++) TRW = TRW + TRWL[j];
			}
		}

		PCEW = (TRC <= 1.0E-10) ? 1. : LIMIT(0., 1., TRW / TRC);
		LRSTRS = LRAV;
		LDSTRS = LDAV;
		LESTRS = LEAV;
	} else {
		PCEW = 1.;
		LRSTRS = 1.;
		LDSTRS = 1.;
		LESTRS = 1.;
	}
	CPEW = LESTRS;
}
