#include <cmath>
#include <algorithm>
#include <vector>
#include "model.h"
#include "oryzaUtil.h"

namespace {
	const int MNL = 10;
	const double TINY = 1.0E-5;

	void resize_soil_arrays(oryza_soil &s, int nl) {
		auto resize10 = [nl](std::vector<double> &v, double fill) {
			v.assign(MNL, fill);
			if (static_cast<int>(v.size()) < nl) v.resize(nl, fill);
		};
		resize10(s.WCL, 0.3);
		resize10(s.WCLQT, 0.3);
		resize10(s.MSKPA, 0.);
		resize10(s.WL, 0.);
		resize10(s.WLFC, 0.);
		resize10(s.WLAD, 0.);
		resize10(s.WLST, 0.);
		resize10(s.KSAT, 0.);
		resize10(s.CAPRI, 0.);
		resize10(s.GWFILL, 0.);
		resize10(s.WLCH, 0.);
		resize10(s.ZL, 0.);
		resize10(s.MS, 0.);
		resize10(s.MSUC, 0.);
		s.WLFL.assign(nl + 1, 0.);
	}

	void default_paddyin_soil(oryza_soil &s) {
		// Fill only missing fields so R-supplied soil can override
		if (s.SCODE.empty()) s.SCODE = "PADDY";
		if (s.NL <= 0) s.NL = 10;
		if (s.TKL.empty()) s.TKL = {0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.15, 0.15, 0.15, 0.15};
		if (s.ZRTMS <= 0) s.ZRTMS = 1.0;
		if (s.WL0MX <= 0) s.WL0MX = 250.;
		if (s.FIXPERC <= 0) s.FIXPERC = 3.0; // realistic default; paddyin.dat used 10000
		if (s.ZWTB.empty()) s.ZWTB = {1., 200., 366., 200.};
		if (s.KST.empty()) s.KST.assign(s.NL, 42.462);
		if (s.WCST.empty()) s.WCST.assign(s.NL, 0.432);
		if (s.VGA.empty()) s.VGA.assign(s.NL, 0.04106);
		if (s.VGL.empty()) s.VGL.assign(s.NL, -3.447);
		if (s.VGN.empty()) s.VGN.assign(s.NL, 1.17);
		if (s.VGR.empty()) s.VGR.assign(s.NL, 0.01);
		if (s.WCLI.empty()) s.WCLI.assign(s.NL, 0.30);
		if (s.RIWCLI.empty()) s.RIWCLI = "NO";
	}

	double intgrl(double state, double rate, double delt) {
		return state + rate * delt;
	}
}

void sync_soil_to_crop(oryza_model &m) {
	oryza_soil &s = m.soil;
	oryza_crop &c = m.crop;
	m.NL = s.NL;
	c.WL0 = s.WL0;
	c.ZRTMS = s.ZRTMS;
	c.TKLT = s.TKLT;
	c.TKL = s.TKL;
	s.WCLQT = s.WCL;
	c.WCLQT = s.WCL;
	c.WCST = s.WCST;
	c.WCFC = s.WCFC;
	c.WCWP = s.WCWP;
	c.WCAD = s.WCAD;
	c.MSKPA = s.MSKPA;
}


void paddy_initialize(oryza_model &m) {
	oryza_soil &s = m.soil;
	default_paddyin_soil(s);

	s.PUDDLD = (s.SWITPD == 1);
	s.GRWAT = (s.SWITGW == 1 || s.SWITGW == 2);
	s.RWCLI = (s.RIWCLI == "YES");

	int nl = s.NL;
	if (nl > MNL) nl = MNL;
	s.NL = nl;
	resize_soil_arrays(s, nl);

	if (static_cast<int>(s.KST.size()) < nl) s.KST.assign(nl, 42.462);
	if (static_cast<int>(s.WCST.size()) < nl) s.WCST.assign(nl, 0.432);
	if (static_cast<int>(s.VGA.size()) < nl) s.VGA.assign(nl, 0.04106);
	if (static_cast<int>(s.VGL.size()) < nl) s.VGL.assign(nl, -3.447);
	if (static_cast<int>(s.VGN.size()) < nl) s.VGN.assign(nl, 1.17);
	if (static_cast<int>(s.VGR.size()) < nl) s.VGR.assign(nl, 0.01);
	if (static_cast<int>(s.WCLI.size()) < nl) s.WCLI.assign(nl, 0.30);

	s.WCFC.assign(nl, 0.);
	s.WCWP.assign(nl, 0.);
	s.WCAD.assign(nl, 0.);
	s.WCSTRP.assign(nl, 0.);

	for (int i = 0; i < nl; i++) {
		s.KSAT[i] = s.KST[i];
		if (!s.PUDDLD) s.WCSTRP[i] = s.WCST[i];
	}

	if (s.SWITPF == 1) {
		for (int i = 0; i < nl; i++) {
			double wcl = 0., ms = 100.;
			SUWCMS2(i + 1, 2, s, wcl, ms);
			s.WCFC[i] = wcl;
			wcl = 0.; ms = 1.6E4;
			SUWCMS2(i + 1, 2, s, wcl, ms);
			s.WCWP[i] = wcl;
			wcl = 0.; ms = 1.0E7;
			SUWCMS2(i + 1, 2, s, wcl, ms);
			s.WCAD[i] = wcl;
		}
	} else if (s.SWITPF == 0) {
		// DATA mode: WCFC, WCWP, WCAD must be supplied
		if (static_cast<int>(s.WCFC.size()) < nl) s.WCFC.assign(nl, 0.48);
		if (static_cast<int>(s.WCWP.size()) < nl) s.WCWP.assign(nl, 0.21);
		if (static_cast<int>(s.WCAD.size()) < nl) s.WCAD.assign(nl, 0.01);
	}

	s.TKL_mm.assign(nl, 0.);
	s.TKLT = 0.;
	for (int i = 0; i < nl; i++) {
		s.TKL_mm[i] = 1000. * s.TKL[i];
		s.TKLT = s.TKLT + s.TKL[i];
		s.WLFC[i] = s.WCFC[i] * s.TKL_mm[i];
		s.WLAD[i] = s.WCAD[i] * s.TKL_mm[i];
		s.WLST[i] = s.WCST[i] * s.TKL_mm[i];
		s.WL[i] = s.WCLI[i] * s.TKL_mm[i];
	}

	for (int i = 0; i < nl; i++) {
		if (i == 0) s.ZL[i] = 0.;
		else s.ZL[i] = s.ZL[i - 1] + s.TKL_mm[i - 1] / 10.;
	}

	if (s.GRWAT) GWTAB(1, s, m.DOY, m.DELT);

	s.WL0 = s.WL0I;
	s.WCUMI = 0.;
	for (int i = 0; i < nl; i++) {
		s.WCL[i] = s.WCLI[i];
		s.WCUMI = s.WCUMI + s.WL[i];
	}
	s.WCUM = s.WCUMI;
	s.DSPW = 1.;
	s.CRACKS = false;
	s.PERC = 0.;

	s.WCUMCO = 0.;
	s.WL0CO = 0.;
	s.WL0FCUM1 = 0.;
	s.UPRICUM1 = 0.;
	s.GWCUM1 = 0.;
	s.PERCCUM1 = 0.;
	s.CAPTOTCUM1 = 0.;
	s.RAINCUM1 = 0.;
	s.IRCUM1 = 0.;
	s.RUNOFCUM1 = 0.;
	s.EVSWCUM1 = 0.;
	s.TRWCUM1 = 0.;
	s.DRAICUM1 = 0.;

	sync_soil_to_crop(m);
}


void paddy_rate(oryza_model &m) {
	oryza_soil &s = m.soil;
	oryza_crop &c = m.crop;
	const std::string &ESTAB = m.control.ESTAB;
	int nl = s.NL;
	double DELT = m.DELT;
	double RAIN = m.atm.RAIN;
	double EVSC = c.EVSC;
	double TRW = c.TRW;
	std::vector<double> &TRWL = c.TRWL;
	double IR = c.IR;
	int CROPSTA = c.CROPSTA;

	if (static_cast<int>(TRWL.size()) < nl) TRWL.assign(MNL, 0.);

	if (s.RWCLI) {
		if ((ESTAB == "TRANSPLANT" && CROPSTA == 3) ||
		    (ESTAB == "DIRECT-SEED" && CROPSTA == 1)) {
			s.WL0 = std::min(s.WL0I, s.WL0MX);
			for (int i = 0; i < nl; i++) {
				s.WCL[i] = s.WCST[i];
				s.WL[i] = s.WCL[i] * s.TKL_mm[i];
			}
		}
	}

	s.WL0CH = 0.;
	s.WCUMCH = 0.;
	s.RUNOF = 0.;
	s.EVSW = 0.;
	s.EVSWS = 0.;
	s.CAPTOT = 0.;
	s.GWTOT = 0.;
	s.DRAIN = 0.;
	s.PERC = 0.;

	for (int i = 0; i < nl; i++) {
		s.WLFL[i] = 0.;
		s.WLCH[i] = 0.;
		s.CAPRI[i] = 0.;
		s.GWFILL[i] = 0.;
	}
	s.WLFL[nl] = 0.;

	if (s.WL0 >= TINY) {
		s.DSPW = 1.;
		if (s.WL0 / DELT + RAIN + IR >= EVSC + TRW) {
			s.WL0CH = RAIN + IR - EVSC - TRW;
			for (int i = 0; i < nl; i++) TRWL[i] = 0.;
			s.EVSW = EVSC;
			if (s.SWITVP == -1) s.PERCOL = s.FIXPERC;
			else if (s.SWITVP == 0 && !s.PERTB.empty()) s.PERCOL = AFGEN(s.PERTB, s.ZW);
			else if (s.SWITVP == 2 && !s.PTABLE.empty()) s.PERCOL = AFGEN(s.PTABLE, s.TIME);
			else s.PERCOL = 0.;
			if (s.WL0 / DELT + s.WL0CH >= s.PERCOL) s.PERC = s.PERCOL;
			else s.PERC = s.WL0 / DELT + s.WL0CH;
			s.WL0CH = s.WL0CH - s.PERC;
			if (s.WL0 + s.WL0CH * DELT >= s.WL0MX) {
				s.RUNOF = (s.WL0 + s.WL0CH * DELT - s.WL0MX) / DELT;
				s.WL0CH = s.WL0CH - s.RUNOF;
			}
			for (int i = 0; i <= nl; i++) s.WLFL[i] = s.PERC;
		} else if (s.WL0 / DELT + RAIN + IR >= EVSC &&
		           s.WL0 / DELT + RAIN + IR < EVSC + TRW) {
			s.WL0CH = -s.WL0 / DELT;
			for (int i = 0; i < nl; i++) {
				if (TRW > 1.E-10) {
					TRWL[i] = ((TRW + EVSC - RAIN - IR - s.WL0 / DELT) / TRW) * TRWL[i] * DELT;
				}
			}
			s.EVSW = EVSC;
		} else {
			s.WL0CH = -s.WL0 / DELT;
			s.WLFL[0] = RAIN + IR;
			s.EVSWS = std::min(EVSC + s.WL0CH, s.WL[0] / DELT - s.WLAD[0] / DELT + RAIN + IR);
			s.EVSW = s.WL0 / DELT + s.EVSWS;
		}
	} else {
		double EVSH = std::min(EVSC, std::max(0., (s.WL[0] - s.WLAD[0]) / DELT + RAIN + IR));
		double EVSD = std::min(EVSC, 0.6 * EVSC * (std::sqrt(s.DSPW) - std::sqrt(s.DSPW - 1.)) + RAIN + IR);
		s.EVSW = INSW(s.DSPW - 1.1, EVSH, EVSD);
		s.EVSW = std::min(s.EVSW, std::max(0., RAIN + IR + (s.WL[0] - s.WLAD[0]) / DELT));
		s.EVSWS = s.EVSW;
		s.DSPW = s.DSPW + 1.;
		s.WLFL[0] = RAIN + IR;
		for (int i = 0; i < nl; i++) {
			s.WLFL[i + 1] = DOWNFL(i + 1, s.KSAT[i], s.WLFL[i], TRWL[i], s.EVSWS, s.WL[i], s.WLFC[i], DELT);
		}
		double REST = 0., FLNEW = 0., HLP = 0.;
		for (int i = nl - 1; i >= 0; i--) {
			BACKFL(i + 1, s.WL[i], s.WLFL[i], s.WLFL[i + 1], s.EVSWS, TRWL[i], s.WLST[i], DELT, FLNEW, HLP);
			s.WLFL[i] = FLNEW;
			if (i == 0) REST = HLP;
		}
		s.WL0CH = std::max(0., (REST - s.WLST[0]) / DELT);
	}

	for (int i = 0; i < nl; i++) {
		if (i == 0) {
			s.WLCH[i] = s.WLFL[i] - s.WLFL[i + 1] - TRWL[i] - s.EVSWS + s.CAPRI[i] + s.GWFILL[i];
		} else {
			s.WLCH[i] = s.WLFL[i] - s.WLFL[i + 1] - TRWL[i] + s.CAPRI[i] + s.GWFILL[i];
		}
		s.WCUMCH = s.WCUMCH + s.WLCH[i];
	}
}


void paddy_state(oryza_model &m) {
	oryza_soil &s = m.soil;
	int nl = s.NL;
	double DELT = m.DELT;

	if (s.GRWAT) {
		s.ZWPREV = s.ZW;
		GWTAB(3, s, m.DOY, DELT);
	}

	s.WL0 = intgrl(s.WL0, s.WL0CH + s.WL0FILL, DELT);

	for (int i = 0; i < nl; i++) {
		s.WCL[i] = intgrl(s.WCL[i], s.WLCH[i] / s.TKL_mm[i], DELT);
		s.WCL[i] = std::min(s.WCST[i], std::max(s.WCL[i], s.WCAD[i]));
		s.WL[i] = s.WCL[i] * s.TKL_mm[i];
	}

	for (int i = 0; i < nl; i++) {
		if (s.SWITPF == 1) {
			double wcl = s.WCL[i];
			SUWCMS2(i + 1, 1, s, wcl, s.MSUC[i]);
		} else {
			double FACT = 0.;
			if (s.WCL[i] >= s.WCFC[i]) {
				FACT = std::max(0., std::min(1., (s.WCST[i] - s.WCL[i]) / (s.WCST[i] - s.WCFC[i])));
				s.MSUC[i] = std::pow(10., FACT * 2.0);
			} else if (s.WCL[i] >= s.WCWP[i]) {
				FACT = std::max(0., std::min(1., (s.WCL[i] - s.WCWP[i]) / (s.WCFC[i] - s.WCWP[i])));
				s.MSUC[i] = std::pow(10., 4.2 - FACT * 2.2);
			} else {
				FACT = std::max(0., std::min(1., (s.WCL[i] - s.WCAD[i]) / (s.WCWP[i] - s.WCAD[i])));
				s.MSUC[i] = std::pow(10., 7.0 - FACT * 2.8);
			}
		}
		s.MSKPA[i] = s.MSUC[i] / 10.;
	}

	s.IRCUM1 = s.IRCUM1 + m.crop.IR * DELT;
	s.RAINCUM1 = s.RAINCUM1 + m.atm.RAIN * DELT;
	s.TRWCUM1 = s.TRWCUM1 + m.crop.TRW * DELT;
	s.EVSWCUM1 = s.EVSWCUM1 + s.EVSW * DELT;
	s.PERCCUM1 = s.PERCCUM1 + s.PERC * DELT;
	s.WCUM = s.WCUM + s.WCUMCH * DELT;
	s.TIME = s.TIME + DELT;
	sync_soil_to_crop(m);
}
