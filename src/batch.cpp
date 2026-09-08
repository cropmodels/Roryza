/*
Spatial / batch ORYZA runs (RWofost-style)
*/

#include <cmath>
#include <vector>
#include "model.h"
#include "Rcpp.h"

std::vector<double> oryza_model::run_batch(
	std::vector<double> tmin, std::vector<double> tmax, std::vector<double> srad,
	std::vector<double> prec, std::vector<double> vapr, std::vector<double> wind,
	std::vector<long> date, std::vector<long> mstart, std::vector<int> soilindex,
	OryzaSoilCollection soils, std::vector<double> depth,
	std::vector<double> elevation, std::vector<double> latitude) {

	bool watlim = control.water_limited;

	size_t sz = date.size();
	if (sz == 0 || tmin.size() < sz) {
		return {};
	}
	size_t nc = tmin.size() / sz;
	size_t nsim = mstart.size();
	std::vector<double> out(nc * nsim, NAN);

	int nsoils = static_cast<int>(soils.soils.size());
	if (nsoils == 0) {
		Rcpp::Rcout << "bad soil data" << std::endl;
		return out;
	}

	bool varsoils = false;
	if (watlim && (nsoils > 1)) {
		if (nc != soilindex.size()) {
			Rcpp::Rcout << "bad soil index data" << std::endl;
			return out;
		}
		varsoils = true;
	}

	control.output_option = "BATCH";
	wth.date = date;

	auto slice = [&](const std::vector<double> &src, size_t offset) {
		if (src.size() >= offset + sz) {
			return std::vector<double>(src.begin() + offset, src.begin() + offset + sz);
		}
		return std::vector<double>(sz, 0.0);
	};

	for (size_t i = 0; i < nc; i++) {
		if (i < latitude.size()) control.latitude = latitude[i];
		if (i < elevation.size()) control.elevation = elevation[i];

		size_t offset = sz * i;
		if (std::isnan(tmin[offset])) {
			continue;
		}

		int sidx = 0;
		if (varsoils) {
			sidx = soilindex[i] - 1;
			if ((sidx < 0) || (sidx >= nsoils)) {
				continue;
			}
		}
		// Always copy from the template so water-balance state does not carry over cells
		soil = soils.soils[sidx];
		if (watlim && i < depth.size() && depth[i] >= 0) {
			soil.ZRTMS = depth[i];
		}

		wth.tmin = slice(tmin, offset);
		wth.tmax = slice(tmax, offset);
		wth.srad = slice(srad, offset);
		wth.prec = slice(prec, offset);
		wth.vapr = slice(vapr, offset);
		wth.wind = slice(wind, offset);

		for (size_t j = 0; j < nsim; j++) {
			control.modelstart = mstart[j];
			double yield = NAN;
			try {
				run();
				if (!output.values.empty()) {
					yield = output.values.back();
				}
			} catch (...) {
				yield = NAN;
			}
			out[j * nc + i] = yield;
		}
	}
	return out;
}
