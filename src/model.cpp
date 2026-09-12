#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include "model.h"


bool oryza_model::weather_step() {
	if (time >= wth.tmin.size()) {
		fatalError = true;
		messages.push_back("reached end of weather data");
		return false;
	}
	if (std::isnan(wth.tmin[time]) || std::isnan(wth.tmax[time]) ||
	    std::isnan(wth.prec[time]) || std::isnan(wth.srad[time]) ||
	    std::isnan(wth.vapr[time]) || std::isnan(wth.wind[time])) {
		fatalError = true;
		messages.push_back("missing value in weather data");
		return false;
	}

	atm.TMMN = wth.tmin[time];
	atm.TMMX = wth.tmax[time];
	atm.TMDA = (atm.TMMN + atm.TMMX) / 2.;
	// RWofost weather: srad in kJ m-2 d-1; ORYZA uses J m-2 d-1
	crop.RDD = wth.srad[time] * 1000.;
	atm.WN = wth.wind[time];
	atm.VP = wth.vapr[time]; // kPa
	atm.RAIN = wth.prec[time]; // mm

	if (wth.date.size() > time) {
		DOY = doy_from_days(wth.date[time]);
		IDOY = static_cast<int>(DOY);
	} else {
		IDOY = static_cast<int>(time + 1);
		DOY = IDOY;
	}
	atm.latitude = control.latitude;

	double solcon, angot, dsinb, dsinbe, sinld, cosld;
	SASTRO(IDOY, atm.latitude, solcon, angot, crop.DAYL, crop.DAYLP, dsinb, dsinbe, sinld, cosld);
	return true;
}


void oryza_model::model_output() {
	if (control.output_option == "BATCH") {
		output.values.push_back(crop.WSO);
	} else if (control.NITROENV == "NITROGEN BALANCE") {
		output.values.insert(output.values.end(), {
			double(step), crop.DVS, crop.LAI,
			crop.WRT, crop.WLV, crop.WST, crop.WSO, crop.WRR14,
			crop.TRC, crop.EVSC, crop.DAE, double(crop.CROPSTA),
			crop.TNSOIL, crop.NACR, XFERT, NFERTP
		});
	} else if (control.PRODENV == "WATER BALANCE") {
		double msk1 = crop.MSKPA.empty() ? 0. : crop.MSKPA[0];
		output.values.insert(output.values.end(), {
			double(step), crop.DVS, crop.LAI,
			crop.WRT, crop.WLV, crop.WST, crop.WSO, crop.WRR14,
			crop.TRC, crop.TRW, crop.EVSC, crop.WL0, crop.IR, msk1,
			crop.DAE, double(crop.CROPSTA), crop.PCEW, crop.LESTRS
		});
	} else {
		output.values.insert(output.values.end(), {
			double(step), crop.DVS, crop.LAI,
			crop.WRT, crop.WLV, crop.WST, crop.WSO, crop.WRR14,
			crop.TRC, crop.EVSC, crop.DAE, double(crop.CROPSTA)
		});
	}
}


void oryza_model::update_cropsta() {
	if (crop.CROPSTA == 3) crop.CROPSTA = 4;

	if (crop.CROPSTA == 2) {
		if (crop.DAE == double(crop.SBDUR)) crop.CROPSTA = 3;
		if (control.ESTAB == "DIRECT-SEED") crop.CROPSTA = 4;
	}
	if (crop.CROPSTA == 1) {
		if (control.ESTAB == "TRANSPLANT") {
			crop.CROPSTA = 2;
		} else if (control.ESTAB == "DIRECT-SEED") {
			crop.CROPSTA = 4;
		}
	}
}


void oryza_model::run() {
	step = 1;
	time = 0;
	TERMINAL = false;
	fatalError = false;
	messages.clear();
	output.values.clear();

	if (control.output_option == "BATCH") {
		output.names = {"WSO"};
	} else if (control.nitrogen_limited) {
		output.names = {"step", "DVS", "LAI", "WRT", "WLV", "WST", "WSO", "WRR14",
			"TRC", "EVSC", "DAE", "CROPSTA", "TNSOIL", "NACR", "XFERT", "NFERTP"};
	} else if (control.water_limited) {
		output.names = {"step", "DVS", "LAI", "WRT", "WLV", "WST", "WSO", "WRR14",
			"TRC", "TRW", "EVSC", "WL0", "IR", "MSKPA1", "DAE", "CROPSTA", "PCEW", "LESTRS"};
	} else {
		output.names = {"step", "DVS", "LAI", "WRT", "WLV", "WST", "WSO", "WRR14", "TRC", "EVSC", "DAE", "CROPSTA"};
	}

	if (wth.tmin.empty()) {
		messages.push_back("no weather data");
		fatalError = true;
		return;
	}
	if (control.modelstart != 0 && !wth.date.empty()) {
		if (control.modelstart < wth.date[0] || control.modelstart > wth.date.back()) {
			messages.push_back("modelstart outside weather date range");
			fatalError = true;
			return;
		}
		time = 0;
		while (time < wth.date.size() && wth.date[time] < control.modelstart) {
			time++;
		}
	}

	if (control.water_limited) {
		control.PRODENV = "WATER BALANCE";
	} else {
		control.PRODENV = "POTENTIAL";
	}
	if (control.nitrogen_limited) {
		control.NITROENV = "NITROGEN BALANCE";
	} else {
		control.NITROENV = "POTENTIAL";
	}

	model_initialize();
	if (fatalError) return;

	unsigned maxstep = static_cast<unsigned>(control.max_duration);
	// cropstart is days after modelstart (0 = emerge on first day)
	unsigned emerge_step = 1 + control.cropstart;

	while ((!TERMINAL) && (step <= maxstep)) {
		if (!weather_step()) break;

		// Emergence day: CROPSTA 0 -> 1
		if (crop.CROPSTA == 0 && step >= emerge_step) {
			crop.CROPSTA = 1;
		}

		model_rate();
		model_state();
		model_output();

		time++;
		step++;
		if (fatalError) break;
	}
}


void oryza_model::model_initialize() {
	DELT = 1.;
	TERMINAL = false;
	fatalError = false;
	NL = crop.NLXM;

	crop.ANGA = control.ANGSTA;
	crop.ANGB = control.ANGSTB;
	crop.FAOF = control.FAOF;
	crop.CO2 = control.CO2;
	crop.TMPSB = control.TMPSB;
	crop.SBDUR = control.SBDUR;
	//crop.TMCTB = control.TMCTB;
	atm.latitude = control.latitude;

	if (control.ESTAB == "DIRECT-SEED") {
		crop.SBDUR = 0;
	}

	crop.CROPSTA = 0;
	crop.RAINCU = 0.;

	ET2_initialize();

	if (control.PRODENV == "WATER BALANCE" && control.WATBAL == "PADDY") {
		paddy_initialize(*this);
		irrig_initialize(*this);
		WSTRESS_initialization(DELT, crop.TRC, crop.ZRT, crop.TKL, NL, crop.CROPSTA,
		                       crop.WCLQT, crop.WCWP, crop.WCAD, crop.MSKPA, crop, control.ESTAB,
		                       crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
	} else if (control.PRODENV == "POTENTIAL") {
		WNOSTRESS(NL, crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
		crop.TKLT = 100.;
		crop.ZRTMS = 100.;
		crop.WL0 = 0.;
		NL = crop.NLXM;
		for (int i = 0; i < NL; i++) {
			crop.WCLQT[i] = 0.3;
			crop.WCST[i] = 0.3;
		}
	}

	oryza_initialize();

	if (control.NITROENV == "POTENTIAL") {
		NNOSTRESS2_initialization(crop.NFLVI, crop.NMAXLT, crop.NFLVTB, DELT, crop.CROPSTA,
		                          crop.DVS, crop.WLVG, crop.LAI, crop.SLA, crop.NFLV, crop.NSLLV, crop.RNSTRS);
	} else if (control.NITROENV == "NITROGEN BALANCE") {
		nsoil_initialize(*this);
		ncrop2_initialize(*this);
	}
}


void oryza_model::model_rate() {
	atm.TMDA = (atm.TMMX + atm.TMMN) / 2.;

	if (control.PRODENV == "WATER BALANCE" && control.WATBAL == "PADDY") {
		sync_soil_to_crop(*this);
	}

	ET2_rate(crop.ANGA, crop.ANGB, crop.RDD, atm.TMDA, atm.VP, atm.WN, atm.latitude, IDOY,
	         control.ETMOD, crop.CROPSTA, crop.FAOF, crop.WL0, crop.WCLQT, crop.WCST, crop.LAI,
	         crop.EVSC, crop.ETD, crop.TRC);

	if (control.PRODENV == "WATER BALANCE" && control.WATBAL == "PADDY") {
		WSTRESS_rate(DELT, crop.TRC, crop.ZRT, crop.TKL, NL, crop.CROPSTA,
		             crop.WCLQT, crop.WCWP, crop.WCAD, crop.MSKPA,
		             crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
	} else if (control.PRODENV == "POTENTIAL") {
		WNOSTRESS(NL, crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
	}

	oryza_rate();

	if (control.PRODENV == "WATER BALANCE" && control.WATBAL == "PADDY") {
		irrig_rate(*this);
		paddy_rate(*this);
	}

	if (control.NITROENV == "POTENTIAL") {
		NNOSTRESS2_rate(crop.NFLVI, crop.NMAXLT, crop.NFLVTB, DELT, crop.CROPSTA,
		                crop.DVS, crop.WLVG, crop.LAI, crop.SLA, crop.NFLV, crop.NSLLV, crop.RNSTRS);
	} else if (control.NITROENV == "NITROGEN BALANCE") {
		// MODELS order: NCROP2 (uses TNSOIL → NACR) then NSOIL rate
		ncrop2_rate(*this);
		nsoil_rate(*this);
	}

	if (control.PRODENV == "POTENTIAL") {
		crop.TKLT = 100.;
		crop.ZRTMS = 100.;
		crop.WL0 = 0.;
		NL = crop.NLXM;
		for (int i = 0; i < NL; i++) {
			crop.WCLQT[i] = 0.3;
			crop.WCST[i] = 0.3;
		}
	}
}


void oryza_model::model_state() {
	crop.RAINCU = crop.RAINCU + atm.RAIN;

	if (control.PRODENV == "WATER BALANCE" && control.WATBAL == "PADDY") {
		paddy_state(*this);
		irrig_state(*this);
	}

	oryza_state();
	ET2_state(DELT, crop.CROPSTA, control.ESTAB, crop.ETD, crop.EVSC, crop.TRC);
	if (control.NITROENV == "NITROGEN BALANCE") {
		ncrop2_state(*this);
		nsoil_state(*this);
	}
	update_cropsta();
}
