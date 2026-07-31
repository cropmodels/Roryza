#include "model.h"
#include "oryzaUtil.h"

namespace {
	int WL0CNT = 0;
	double IRC = 0.;
	double IRCU = 0.;
}

void irrig_initialize(oryza_model &m) {
	(void)m;
	m.crop.IR = 0.;
	IRC = 0.;
	IRCU = 0.;
	WL0CNT = 0;
	if (m.control.ISTAGET.size() < 15) {
		m.control.ISTAGET = std::vector<double>(15, 10.);
	}
}


void irrig_rate(oryza_model &m) {
	oryza_control &c = m.control;
	oryza_crop &crp = m.crop;
	oryza_soil &s = m.soil;

	crp.IR = 0.;

	if (crp.DVS >= c.DVSIMAX) return;

	int SWITIR = c.SWITIR;

	if (SWITIR == 0) {
		crp.IR = 0.;
	} else if (SWITIR == 1) {
		if (!c.RIRRIT.empty()) {
			crp.IR = AFGEN(c.RIRRIT, m.DOY);
		}
	} else if (SWITIR == 2) {
		if (s.WL0 <= c.WL0MIN) crp.IR = c.IRRI;
	} else if (SWITIR == 3) {
		int sl = std::max(1, std::min(c.SLMIN, s.NL)) - 1;
		if (s.MSKPA[sl] >= c.KPAMIN) crp.IR = c.IRRI;
	} else if (SWITIR == 4) {
		int sl = std::max(1, std::min(c.SLMIN, s.NL)) - 1;
		if (s.WCL[sl] <= c.WCMIN) crp.IR = c.IRRI;
	} else if (SWITIR == 5) {
		if (s.WL0 <= 1.) {
			if (WL0CNT == c.WL0DAY) {
				crp.IR = c.IRRI;
				WL0CNT = 0;
			} else {
				crp.IR = 0.;
				WL0CNT = WL0CNT + static_cast<int>(m.DELT);
			}
		} else {
			crp.IR = 0.;
		}
	} else if (SWITIR == 6) {
		const std::vector<double> &st = c.ISTAGET;
		int sl = std::max(1, std::min(c.SLMIN, s.NL)) - 1;
		for (int p = 0; p < 5; p++) {
			int base = p * 3;
			if (base + 2 >= static_cast<int>(st.size())) break;
			if (crp.DVS > st[base] && crp.DVS <= st[base + 1]) {
				if (s.MSKPA[sl] >= st[base + 2]) crp.IR = c.IRRI;
				break;
			}
		}
	}
}


void irrig_state(oryza_model &m) {
	IRCU = IRCU + m.crop.IR * m.DELT;
	if (m.crop.CROPSTA >= 3) {
		IRC = IRC + m.crop.IR * m.DELT;
	}
}
