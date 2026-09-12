//----------------------------------------------------------------------//
// SUBROUTINE NCROP2                                                    //
// Port of ORYZA2000 Ncrop2.f90 (Bouman, Dec 2001 / Aug 2003)            //
// Nitrogen dynamics in the rice crop and effects on growth.            //
//----------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <string>
#include "model.h"
#include "oryzaUtil.h"


void SUBNBC(double CHKIN, double CKCFL, double TIME, double &NBCHK, bool &TERMINAL) {
	(void)TIME;
	NBCHK = 2.0 * (CHKIN - CKCFL) / (CHKIN + CKCFL + 1.E-10);
	if (std::fabs(NBCHK) > 0.001) {
		TERMINAL = true;
	}
}


void ncrop2_initialize(oryza_model &m) {
	oryza_crop &c = m.crop;

	m.ANLV = 0.;
	m.ANSO = 0.;
	m.ANST = 0.;
	m.ANLD = 0.;
	m.ANCR = 0.;
	m.ANLVA = 0.;
	m.ANSTA = 0.;
	m.ANCRF = 0.;
	m.NALVS = 0.;
	m.NASTS = 0.;
	m.NASOS = 0.;
	m.NACRS = 0.;
	m.NTRTS = 0.;
	m.NALV = 0.;
	m.NAST = 0.;
	m.NASO = 0.;
	c.NACR = 0.;
	m.NLV = 0.;
	m.NST = 0.;
	m.NSO = 0.;
	m.NLDLV = 0.;
	m.NLVAN = 0.;
	m.NSTAN = 0.;
	m.NTRT = 0.;

	m.FNLV = c.FNLVI;
	m.FNST = 0.5 * c.FNLVI;
	m.FNSO = 0.;
	c.NFLV = c.NFLVI;
	c.NSLLV = 1.;
	c.RNSTRS = 1.;

	m.NMINSO = AFGEN(c.NMINSOT, m.ANCRF);
	m.NMAXL = AFGEN(c.NMAXLT, c.DVS);
	m.NMINL = AFGEN(c.NMINLT, c.DVS);
	m.FNLV = c.FNLVI;
	c.NFLV = c.NFLVI;
}


void ncrop2_rate(oryza_model &m) {
	oryza_crop &c = m.crop;
	const double DELT = m.DELT;

	m.NMINSO = AFGEN(c.NMINSOT, m.ANCRF);
	m.NMAXL = AFGEN(c.NMAXLT, c.DVS);
	m.NMINL = AFGEN(c.NMINLT, c.DVS);

	if (c.CROPSTA < 4) {
		return;
	}

	// Potential leaf N content on LAI basis (kept for FORTRAN parity; unused in NFLV path)
	(void)AFGEN(c.NFLVTB, c.DVS);

	// Maximum N demand of leaves, stems, storage organs
	double NDEML = (m.NMAXL * (c.WLVG + c.GLV * DELT) - m.ANLV) / DELT;
	if (NDEML < 0.) NDEML = 0.;

	double NDEMS = (m.NMAXL * 0.5 * (c.WST + c.GST * DELT) - m.ANST) / DELT;
	if (NDEMS < 0.) NDEMS = 0.;

	double NDEMSX = c.NMAXSO * c.GSO;
	if (NDEMSX < 0.) NDEMSX = 0.;

	double NDEMSN = m.NMINSO * c.GSO;
	if (NDEMSN < 0.) NDEMSN = 0.;

	// Translocation to storage organs (after DVS 0.95)
	double ATNLV = 0., ATNST = 0., ATNRT = 0., ATN = 0., NTSO = 0.;
	if (c.DVS >= 0.95) {
		ATNLV = std::max(0., m.ANLV - c.WLVG * c.RFNLV);
		ATNST = std::max(0., m.ANST - c.WST * c.RFNST);
		ATNRT = (ATNLV + ATNST) * c.FNTRT;
		ATN = ATNLV + ATNST + ATNRT;
		NTSO = ATN / c.TCNTRF;
		NTSO = LIMIT(NDEMSN, NDEMSX, NTSO);
	}

	double NTLV = NTSO * ATNLV / NOTNUL(ATN);
	double NTST = NTSO * ATNST / NOTNUL(ATN);
	m.NTRT = NTSO * ATNRT / NOTNUL(ATN);

	double NUPP = std::min(c.NMAXUP, c.TNSOIL);
	if (NUPP < 0.) NUPP = 0.;

	double NDEMC = (NDEML + NTLV) + (NDEMS + NTST) + (NDEMSX - NTSO);

	m.NALV = std::max(0., std::min(NDEML + NTLV, NUPP * ((NDEML + NTLV) / NOTNUL(NDEMC))));
	m.NAST = std::max(0., std::min(NDEMS + NTST, NUPP * ((NDEMS + NTST) / NOTNUL(NDEMC))));
	m.NASO = std::max(0., std::min(NDEMSX - NTSO, NUPP * ((NDEMSX - NTSO) / NOTNUL(NDEMC))));
	c.NACR = m.NALV + m.NAST + m.NASO;

	double NSHKLV = m.ANLV * (1. - c.PLTR);
	double NSHKST = m.ANST * (1. - c.PLTR);
	m.NLDLV = (c.LLV + c.DLDR) * c.RFNLV;
	m.NLDLV = std::max(m.NLDLV, m.ANLV - m.NMAXL * (c.WLVG + c.GLV - c.LLV - c.DLDR));

	m.NLV = m.NALV - NTLV - m.NLDLV - NSHKLV;
	m.NST = m.NAST - NTST - NSHKST;
	m.NSO = NTSO + m.NASO;

	if (c.DVS < 1.) {
		m.NSTAN = m.NST;
		m.NLVAN = m.NLV;
	} else {
		m.NSTAN = 0.;
		m.NLVAN = 0.;
	}
}


void ncrop2_state(oryza_model &m) {
	oryza_crop &c = m.crop;
	const double DELT = m.DELT;

	m.ANSO = m.ANSO + m.NSO * DELT;
	m.ANLV = m.ANLV + m.NLV * DELT;
	m.ANST = m.ANST + m.NST * DELT;
	m.ANLD = m.ANLD + m.NLDLV * DELT;
	m.ANCR = m.ANSO + m.ANLV + m.ANLD + m.ANST;

	m.ANLVA = m.ANLVA + m.NLVAN * DELT;
	m.ANSTA = m.ANSTA + m.NSTAN * DELT;
	m.ANCRF = m.ANSTA + m.ANLVA;

	m.NALVS = m.NALVS + m.NALV * DELT;
	m.NASTS = m.NASTS + m.NAST * DELT;
	m.NASOS = m.NASOS + m.NASO * DELT;
	m.NACRS = m.NALVS + m.NASTS + m.NASOS;

	m.NTRTS = m.NTRTS + m.NTRT * DELT;

	double NBCHK = 0.;
	SUBNBC(m.ANCR, m.NACRS + m.NTRTS, static_cast<double>(m.time), NBCHK, m.TERMINAL);
	if (m.TERMINAL && std::fabs(NBCHK) > 0.001) {
		m.messages.push_back("Error in Nitrogen Balance (NCROP2)");
		m.fatalError = true;
	}

	if (c.CROPSTA < 4) {
		m.FNLV = c.FNLVI;
		m.FNST = 0.5 * c.FNLVI;
		m.FNSO = 0.;
		c.NFLV = c.NFLVI;
		c.NSLLV = 1.;
		c.RNSTRS = 1.;
	} else {
		m.FNLV = m.ANLV / NOTNUL(c.WLVG);
		m.FNST = m.ANST / NOTNUL(c.WST);
		m.FNSO = m.ANSO / NOTNUL(c.WSO);

		if (c.LAI == 0.) {
			c.NFLV = c.NFLVI;
		} else {
			// INTGR2 without forcing returns NFLV1 = FNLV/(10*SLA)
			c.NFLV = m.FNLV / (10. * c.SLA);
		}

		double ANCRPT = c.WLVG * m.NMAXL + c.WST * m.NMAXL * 0.5 + c.WSO * c.NMAXSO;
		double NSTRES = (m.ANCR == 0.) ? 2. : (ANCRPT / m.ANCR);
		if (NSTRES < 1.) NSTRES = 1.;
		if (NSTRES > 2.) NSTRES = 2.;
		c.NSLLV = AFGEN(c.NSLLVT, NSTRES);

		c.RNSTRS = (m.FNLV - 0.9 * m.NMAXL) / (m.NMAXL - 0.9 * m.NMAXL);
		if (c.RNSTRS > 1.) c.RNSTRS = 1.;
		if (c.RNSTRS < 0.) c.RNSTRS = 0.;
	}

	if (c.LAI > 1. && m.FNLV <= 0.5 * m.NMINL) {
		m.messages.push_back("Leaf N < 0.5*MINIMUM; simulation stopped");
		m.TERMINAL = true;
	}
}
