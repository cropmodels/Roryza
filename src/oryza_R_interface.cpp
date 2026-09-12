/*
Roryza one-shot R interface
*/

#include <Rcpp.h>
using namespace Rcpp;
#include "R_interface_util.h"
#include "model.h"


static std::vector<double> flatVec(List lst, const char *s) {
	if (!lst.containsElementNamed(s)) return {};
	SEXP x = lst[s];
	return as<std::vector<double>>(x);
}


static void set_crop_from_list(oryza_crop &crp, List crop) {
	crp.TBD = valueFromList<double>(crop, "TBD");
	crp.TBLV = valueFromList<double>(crop, "TBLV");
	crp.TMD = valueFromList<double>(crop, "TMD");
	crp.TOD = valueFromList<double>(crop, "TOD");
	crp.DVRJ = valueFromList<double>(crop, "DVRJ");
	crp.DVRI = valueFromList<double>(crop, "DVRI");
	crp.DVRP = valueFromList<double>(crop, "DVRP");
	crp.DVRR = valueFromList<double>(crop, "DVRR");
	crp.MOPP = valueFromList<double>(crop, "MOPP");
	crp.PPSE = valueFromList<double>(crop, "PPSE");
	crp.SHCKD = valueFromList<double>(crop, "SHCKD");
	crp.COLDMIN = valueFromList<double>(crop, "COLDMIN");
	crp.COLDEAD = valueFromList<double>(crop, "COLDEAD");
	crp.RGRLMX = valueFromList<double>(crop, "RGRLMX");
	crp.RGRLMN = valueFromList<double>(crop, "RGRLMN");
	crp.SHCKL = valueFromList<double>(crop, "SHCKL");
	crp.SWISLA = valueFromList<std::string>(crop, "SWISLA");
	crp.ASLA = valueFromList<double>(crop, "ASLA");
	crp.BSLA = valueFromList<double>(crop, "BSLA");
	crp.CSLA = valueFromList<double>(crop, "CSLA");
	crp.DSLA = valueFromList<double>(crop, "DSLA");
	crp.SLAMAX = valueFromList<double>(crop, "SLAMAX");
	crp.SLATB = TableFromList(crop, "SLATB");
	crp.SSGATB = TableFromList(crop, "SSGATB");
	crp.FRPAR = valueFromList<double>(crop, "FRPAR");
	crp.SCP = valueFromList<double>(crop, "SCP");
	crp.CO2REF = valueFromList<double>(crop, "CO2REF");
	crp.CO2 = valueFromList<double>(crop, "CO2");
	crp.KDFTB = TableFromList(crop, "KDFTB");
	crp.KNFTB = TableFromList(crop, "KNFTB");
	crp.EFFTB = TableFromList(crop, "EFFTB");
	crp.REDFTT = TableFromList(crop, "REDFTT");
	crp.NFLVTB = TableFromList(crop, "NFLVTB");
	crp.MAINLV = valueFromList<double>(crop, "MAINLV");
	crp.MAINST = valueFromList<double>(crop, "MAINST");
	crp.MAINSO = valueFromList<double>(crop, "MAINSO");
	crp.MAINRT = valueFromList<double>(crop, "MAINRT");
	crp.TREF = valueFromList<double>(crop, "TREF");
	crp.Q10 = valueFromList<double>(crop, "Q10");
	crp.CRGLV = valueFromList<double>(crop, "CRGLV");
	crp.CRGST = valueFromList<double>(crop, "CRGST");
	crp.CRGSO = valueFromList<double>(crop, "CRGSO");
	crp.CRGRT = valueFromList<double>(crop, "CRGRT");
	crp.CRGSTR = valueFromList<double>(crop, "CRGSTR");
	crp.LRSTR = valueFromList<double>(crop, "LRSTR");
	crp.FSTR = valueFromList<double>(crop, "FSTR");
	crp.TCLSTR = valueFromList<double>(crop, "TCLSTR");
	crp.SPGF = valueFromList<double>(crop, "SPGF");
	crp.WGRMX = valueFromList<double>(crop, "WGRMX");
	crp.FSHTB = TableFromList(crop, "FSHTB");
	crp.FLVTB = TableFromList(crop, "FLVTB");
	crp.FSTTB = TableFromList(crop, "FSTTB");
	crp.FSOTB = TableFromList(crop, "FSOTB");
	crp.DRLVT = TableFromList(crop, "DRLVT");
	crp.FCLV = valueFromList<double>(crop, "FCLV");
	crp.FCST = valueFromList<double>(crop, "FCST");
	crp.FCSO = valueFromList<double>(crop, "FCSO");
	crp.FCRT = valueFromList<double>(crop, "FCRT");
	crp.FCSTR = valueFromList<double>(crop, "FCSTR");
	crp.GZRT = valueFromList<double>(crop, "GZRT");
	crp.ZRTMCW = valueFromList<double>(crop, "ZRTMCW");
	crp.ZRTMCD = valueFromList<double>(crop, "ZRTMCD");
	crp.NFLVI = valueFromList<double>(crop, "NFLVI");
	crp.NMAXLT = TableFromList(crop, "NMAXLT");
	crp.FNLVI = valueFromListDefault<double>(crop, "FNLVI", crp.FNLVI);
	crp.NMAXUP = valueFromListDefault<double>(crop, "NMAXUP", crp.NMAXUP);
	crp.NMAXSO = valueFromListDefault<double>(crop, "NMAXSO", crp.NMAXSO);
	crp.RFNLV = valueFromListDefault<double>(crop, "RFNLV", crp.RFNLV);
	crp.RFNST = valueFromListDefault<double>(crop, "RFNST", crp.RFNST);
	crp.TCNTRF = valueFromListDefault<double>(crop, "TCNTRF", crp.TCNTRF);
	crp.FNTRT = valueFromListDefault<double>(crop, "FNTRT", crp.FNTRT);
	{
		auto t = flatVec(crop, "NMINLT"); if (!t.empty()) crp.NMINLT = t;
		t = flatVec(crop, "NMINSOT"); if (!t.empty()) crp.NMINSOT = t;
		t = flatVec(crop, "NSLLVT"); if (!t.empty()) crp.NSLLVT = t;
	}
	crp.LAPE = valueFromList<double>(crop, "LAPE");
	crp.DVSI = valueFromList<double>(crop, "DVSI");
	crp.WLVGI = valueFromList<double>(crop, "WLVGI");
	crp.WSTI = valueFromList<double>(crop, "WSTI");
	crp.WRTI = valueFromList<double>(crop, "WRTI");
	crp.WSOI = valueFromList<double>(crop, "WSOI");
	crp.ZRTI = valueFromList<double>(crop, "ZRTI");
	crp.ZRTTR = valueFromList<double>(crop, "ZRTTR");
	crp.NH = valueFromList<double>(crop, "NH");
	crp.NPLH = valueFromList<double>(crop, "NPLH");
	crp.NPLSB = valueFromList<double>(crop, "NPLSB");
	crp.NPLDS = valueFromList<double>(crop, "NPLDS");
	crp.ULLS = valueFromListDefault<double>(crop, "ULLS", crp.ULLS);
	crp.LLLS = valueFromListDefault<double>(crop, "LLLS", crp.LLLS);
	crp.ULDL = valueFromListDefault<double>(crop, "ULDL", crp.ULDL);
	crp.LLDL = valueFromListDefault<double>(crop, "LLDL", crp.LLDL);
	crp.ULLE = valueFromListDefault<double>(crop, "ULLE", crp.ULLE);
	crp.LLLE = valueFromListDefault<double>(crop, "LLLE", crp.LLLE);
	crp.ULRT = valueFromListDefault<double>(crop, "ULRT", crp.ULRT);
	crp.LLRT = valueFromListDefault<double>(crop, "LLRT", crp.LLRT);
	if (crop.containsElementNamed("SWIRTR")) crp.SWIRTR = as<std::string>(crop["SWIRTR"]);
	crp.SWIRTRF = valueFromListDefault<double>(crop, "SWIRTRF", crp.SWIRTRF);
}


static void set_soil_from_list(oryza_soil &sol, List soil) {
	if (soil.size() == 0) return;
	if (soil.containsElementNamed("SCODE")) sol.SCODE = as<std::string>(soil["SCODE"]);
	if (soil.containsElementNamed("SWITPD")) sol.SWITPD = valueFromListDefault<int>(soil, "SWITPD", sol.SWITPD);
	if (soil.containsElementNamed("SWITGW")) sol.SWITGW = valueFromListDefault<int>(soil, "SWITGW", sol.SWITGW);
	if (soil.containsElementNamed("SWITPF")) sol.SWITPF = valueFromListDefault<int>(soil, "SWITPF", sol.SWITPF);
	if (soil.containsElementNamed("SWITVP")) sol.SWITVP = valueFromListDefault<int>(soil, "SWITVP", sol.SWITVP);
	if (soil.containsElementNamed("SWITKH")) sol.SWITKH = valueFromListDefault<int>(soil, "SWITKH", sol.SWITKH);
	if (soil.containsElementNamed("NL")) sol.NL = valueFromListDefault<int>(soil, "NL", sol.NL);
	if (soil.containsElementNamed("ZRTMS")) sol.ZRTMS = valueFromListDefault<double>(soil, "ZRTMS", sol.ZRTMS);
	if (soil.containsElementNamed("WL0MX")) sol.WL0MX = valueFromListDefault<double>(soil, "WL0MX", sol.WL0MX);
	if (soil.containsElementNamed("WL0I")) sol.WL0I = valueFromListDefault<double>(soil, "WL0I", sol.WL0I);
	if (soil.containsElementNamed("FIXPERC")) sol.FIXPERC = valueFromListDefault<double>(soil, "FIXPERC", sol.FIXPERC);
	if (soil.containsElementNamed("RIWCLI")) sol.RIWCLI = as<std::string>(soil["RIWCLI"]);
	auto tkl = flatVec(soil, "TKL"); if (!tkl.empty()) sol.TKL = tkl;
	auto kst = flatVec(soil, "KST"); if (!kst.empty()) sol.KST = kst;
	auto wcst = flatVec(soil, "WCST"); if (!wcst.empty()) sol.WCST = wcst;
	auto vga = flatVec(soil, "VGA"); if (!vga.empty()) sol.VGA = vga;
	auto vgl = flatVec(soil, "VGL"); if (!vgl.empty()) sol.VGL = vgl;
	auto vgn = flatVec(soil, "VGN"); if (!vgn.empty()) sol.VGN = vgn;
	auto vgr = flatVec(soil, "VGR"); if (!vgr.empty()) sol.VGR = vgr;
	auto wcli = flatVec(soil, "WCLI"); if (!wcli.empty()) sol.WCLI = wcli;
	auto wcfc = flatVec(soil, "WCFC"); if (!wcfc.empty()) sol.WCFC = wcfc;
	auto wcwp = flatVec(soil, "WCWP"); if (!wcwp.empty()) sol.WCWP = wcwp;
	auto wcad = flatVec(soil, "WCAD"); if (!wcad.empty()) sol.WCAD = wcad;
	auto zwtb = flatVec(soil, "ZWTB"); if (!zwtb.empty()) sol.ZWTB = zwtb;
	auto pertb = flatVec(soil, "PERTB"); if (!pertb.empty()) sol.PERTB = pertb;
	auto ptable = flatVec(soil, "PTABLE"); if (!ptable.empty()) sol.PTABLE = ptable;
}


static void set_control_from_list(oryza_control &cntr, List control) {
	cntr.modelstart = valueFromList<long>(control, "modelstart");
	cntr.cropstart = valueFromList<unsigned>(control, "cropstart");
	cntr.output_option = valueFromListDefault<std::string>(control, "output", "");
	cntr.latitude = valueFromList<double>(control, "latitude");
	cntr.elevation = valueFromListDefault<double>(control, "elevation", 0);
	cntr.CO2 = valueFromListDefault<double>(control, "CO2", 340);
	cntr.ANGSTA = valueFromListDefault<double>(control, "ANGSTA", 0.29);
	cntr.ANGSTB = valueFromListDefault<double>(control, "ANGSTB", 0.45);
	cntr.FAOF = valueFromListDefault<double>(control, "FAOF", 1.0);
	cntr.water_limited = valueFromListDefault<bool>(control, "water_limited", false);
	cntr.nitrogen_limited = valueFromListDefault<bool>(control, "nitrogen_limited", false);
	cntr.max_duration = valueFromListDefault<int>(control, "max_duration", 365);
	cntr.ESTAB = valueFromListDefault<std::string>(control, "ESTAB", "DIRECT-SEED");
	cntr.ETMOD = valueFromListDefault<std::string>(control, "ETMOD", "PENMAN");
	cntr.RICETYPE = valueFromListDefault<std::string>(control, "RICETYPE", "LOWLAND");
	cntr.SBDUR = valueFromListDefault<int>(control, "SBDUR", 0);
	cntr.TMPSB = valueFromListDefault<double>(control, "TMPSB", 0);
	cntr.WATBAL = valueFromListDefault<std::string>(control, "WATBAL", "PADDY");
	cntr.SWITIR = valueFromListDefault<int>(control, "SWITIR", 0);
	cntr.DVSIMAX = valueFromListDefault<double>(control, "DVSIMAX", 2.0);
	cntr.IRRI = valueFromListDefault<double>(control, "IRRI", 75.);
	cntr.SLMIN = valueFromListDefault<int>(control, "SLMIN", 3);
	cntr.KPAMIN = valueFromListDefault<double>(control, "KPAMIN", 5.);
	cntr.WCMIN = valueFromListDefault<double>(control, "WCMIN", 0.30);
	cntr.WL0DAY = valueFromListDefault<int>(control, "WL0DAY", 5);
	cntr.WL0MIN = valueFromListDefault<double>(control, "WL0MIN", 10.);
	if (control.containsElementNamed("TMCTB")) {
		SEXP tm = control["TMCTB"];
		if (Rf_isMatrix(tm)) {
			cntr.TMCTB = TableFromList(control, "TMCTB");
		} else {
			cntr.TMCTB = as<std::vector<double>>(tm);
		}
	}
	if (control.containsElementNamed("RIRRIT")) {
		cntr.RIRRIT = as<std::vector<double>>(control["RIRRIT"]);
	}
	if (control.containsElementNamed("ISTAGET")) {
		cntr.ISTAGET = as<std::vector<double>>(control["ISTAGET"]);
	}
	if (control.containsElementNamed("FERTIL")) {
		cntr.FERTIL = as<std::vector<double>>(control["FERTIL"]);
	}
	if (control.containsElementNamed("RECNIT")) {
		cntr.RECNIT = as<std::vector<double>>(control["RECNIT"]);
	}
	if (control.containsElementNamed("SOILSP")) {
		cntr.SOILSP = valueFromListDefault<double>(control, "SOILSP", cntr.SOILSP);
	}
}


void setSoil(oryza_model *m, List soil) {
	set_soil_from_list(m->soil, soil);
}


// [[Rcpp::export(.oryza)]]
NumericMatrix oryza(List crop, DataFrame weather, List soil, List control) {
	oryza_model m;
	set_control_from_list(m.control, control);
	set_crop_from_list(m.crop, crop);
	set_soil_from_list(m.soil, soil);

	m.wth.tmin = vectorFromDF<double>(weather, "tmin");
	m.wth.tmax = vectorFromDF<double>(weather, "tmax");
	m.wth.srad = vectorFromDF<double>(weather, "srad");
	m.wth.prec = vectorFromDF<double>(weather, "prec");
	m.wth.vapr = vectorFromDF<double>(weather, "vapr");
	m.wth.wind = vectorFromDF<double>(weather, "wind");
	m.wth.date = vectorFromDF<long>(weather, "date");

	m.run();

	if (m.fatalError) {
		for (size_t i = 0; i < m.messages.size(); i++) {
			Rcout << m.messages[i] << std::endl;
		}
	}

	size_t nc = m.output.names.size();
	size_t nr = nc == 0 ? 0 : m.output.values.size() / nc;
	NumericMatrix mat(nr, nc);
	CharacterVector cnames = wrap(m.output.names);
	colnames(mat) = cnames;
	size_t k = 0;
	for (size_t i = 0; i < nr; i++) {
		for (size_t j = 0; j < nc; j++) {
			mat(i, j) = m.output.values[k++];
		}
	}
	return mat;
}
