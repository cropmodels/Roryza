#ifndef ORYZA_MODEL_H_

#define ORYZA_MODEL_H_



#include <algorithm>

#include <vector>

#include <cmath>

#include <string>

struct oryza_control;
struct oryza_crop;
struct oryza_soil;
struct oryza_model;

void ET(int ITASK, double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, double IDOY, std::string ETMOD,

                  double CROPSTA, double NL, double FAOF, double WL0, std::vector<double> WCLQT, std::vector<double> WCST, double LAI,

                  double &EVSC, double &ETD, double &TRC);

void ET2_initialize();

void ET2_rate(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,

              double FAOF, double WL0, const std::vector<double> &WCLQT, const std::vector<double> &WCST, double LAI,

              double &EVSC, double &ETD, double &TRC);

void ET2_state(double DELT, int CROPSTA, const std::string &ESTAB, double ETD, double EVSC, double TRC);

void GPPARGET(double xGAI, double xGAID, double xAmaxIn, double xEffIn, double &xAmaxOut, double &xEffOut);

void GPPARSET( double xCO2, double xKNF, double xNFLV, double xREDFT );



void SUBNBC( double CHKIN, double CKCFL, double TIME, double &NBCHK, bool &TERMINAL );

void NNOSTRESS2_initialization( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS );

void NNOSTRESS2_rate( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS );

double NSOIL( int ITASK, int IUNITD, int IUNITL, std::string FILEIT, double OUTPUT, double DELT, double DAE, double DVS, double NACR );

void nsoil_initialize(oryza_model &m);
void nsoil_rate(oryza_model &m);
void nsoil_state(oryza_model &m);

void ncrop2_initialize(oryza_model &m);
void ncrop2_rate(oryza_model &m);
void ncrop2_state(oryza_model &m);

std::vector<double> PHENOL( double DVS, double DVRJ, double DVRI, double DVRP, double DVRR, double HU, double DAYL, double MOPP, double PPSE, double TS, double SHCKD, int CROPSTA);

void SASTRO(double IDOY, double LAT, double &SOLCON, double &ANGOT, double &DAYL, double &DAYLP, double &DSINB, double &DSINBE,

            double &SINLD, double &COSLD);

double SETMKD(double RDD, double TMDA);

std::vector<double> SETPMD(int IDOY, double LAT, int ISURF, double RF, double ANGA, double ANGB, double TMDI, double RDD, double TMDA, double WN,

                      double VP);

double SETPTD( int IDOY, double LAT, double RF, double RDD, double TMDA );

std::vector<double> SGPC1( double CSLV, double AMAX, double EFF, double ECPDF, double GAI, double SINB, double RDPDR, double RDPDF );

std::vector<double> SGPC2( double CSLV, double AMAX, double EFF, double ECPDF, double GAI, double SINB, double RDPDR, double RDPDF );

void SGPCDT(int IACC, int IDOY, double LAT, double RDD, double FRPAR, double CSLV, double AMAX, double EFF, double ECPDF,

            double GAI, double &DAYL, double &DAYLP, double &GPCDT, double &RAPCDT );

void SGPL( double CSLV, double AMAX1, double EFF1, double ECPDF, double GAI, double GAID, double SINB, double &RDPDR, double &RDPDF,

          double &GPL, double &RAPL);

std::vector<double> SRDPRF( double GAID, double CSLV, double SINB, double ECPDF, double RDPDR, double RDPDF );

void SSKYC( double HOUR, double SOLCON, double FRPAR, double DSINBE, double SINLD, double COSLD, double RDD, double &SINB, double &RDPDR, double &RDPDF );

double SUBCD(int CROPSTA, double TAV, double TIME);

void SUBCD2(double COLDMIN, int CROPSTA, double TAV, double &NCOLD);

double SUBDD( double TMAX, double TMIN, double TBD, double TOD, double TMD );

void SUBGRN( double GCR, double CROPSTA, double LRSTRS, double DVS, double SF1, double SF2, double SPGF,

            double TAV, double TMAX, double NSP, double TIME, double &GNSP, double &GNGR, double &SPFERT, bool &GRAINS);

void SUBLAI2(int CROPSTA, double RWLVG, double DLDR, double TSLV, double HULV, double SHCKL, double LESTRS, double SLA,

             double NH, double NPLH, double NPLSB, double DVS, double LAI, double RGRLMX, double RGRLMN, std::string ESTAB, double &GAI, double &RGRL);

void SUBLAI2(int CROPSTA, double RGRLMX, double RGRLMN, double TSLV, double HULV, double SHCKL, double LESTRS, double RNSTRS, double SLA,

             double NH, double NPLH, double NPLSB, double DVS, double LAI, std::string ESTAB, double RWLVG, double DLDR, double WLVG, double &GLAI, double &RGRL);

void SUBLAI3(int CROPSTA,double RGRLMX, double RGRLMN, double TSLV, double HULV, double SHCKL, double LESTRS, double RNSTRS,

             double SLA, double NH, int NPLH, double NPLSB, double DVS, double LAI, std::string ESTAB, double RWLVG,double DLDR,

             double WLVG, double &GLAI, double &RGRL);

std::vector<double> SVPS1(double TMA);

void WNOSTRESS( int NL, double &TRW, std::vector<double> &TRWL, double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW );

void WSTRESS_initialization(double DELT, double TRC, double ZRT, const std::vector<double> &TKL, int NL, int CROPSTA,

                            const std::vector<double> &WCLQT, const std::vector<double> &WCWP, const std::vector<double> &WCAD,

                            const std::vector<double> &MSKPA, const oryza_crop &crop, const std::string &ESTAB,

                            double &TRW, std::vector<double> &TRWL, double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW);

void WSTRESS_rate(double DELT, double TRC, double ZRT, const std::vector<double> &TKL, int NL, int CROPSTA,

                  const std::vector<double> &WCLQT, const std::vector<double> &WCWP, const std::vector<double> &WCAD,

                  const std::vector<double> &MSKPA, double &TRW, std::vector<double> &TRWL,

                  double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW);



void SUBCBC( double CKCIN, double CKCFL, double TIME, double &CBCHK, bool &TERMNL );



// Soil water balance helpers (SUBSOIL.f90)

void GWTAB(int ITASK, oryza_soil &soil, double DOY, double DELT);

void BACKFL(int I, double WL, double FLIN, double FLOUT, double EVSWS, double TRWL, double WLST, double DELT,

            double &FLNEW, double &HLP);

double DOWNFL(int I, double KSAT, double FLIN, double TRWL, double EVSWS, double WL, double WLFC, double DELT);

void SUWCMS2(int I, int SWIT4, oryza_soil &soil, double &WCL, double &MS);

void SUMSKM2(int I, double MS, oryza_soil &soil, double &KMS);

double SUBSL2(double PF, double D, int I, oryza_soil &soil);

double SATFLX(oryza_soil &soil, double WL0); // stub unless SWITVP=1 puddled

void SUWCHK(double CKWFL, double CKWIN, double TIME);



void irrig_initialize(oryza_model &m);

void irrig_rate(oryza_model &m);

void irrig_state(oryza_model &m);



void paddy_initialize(oryza_model &m);

void paddy_rate(oryza_model &m);

void paddy_state(oryza_model &m);



void sync_soil_to_crop(oryza_model &m);





inline unsigned doy_from_days(long z) {

	z += 719468;

	const long era = (z >= 0 ? z : z - 146096) / 146097;

	const unsigned doe = static_cast<unsigned>(z - era * 146097);

	const unsigned yoe = (doe - doe/1460 + doe/36524 - doe/146096) / 365;

	const long y = static_cast<long>(yoe) + era * 400;

	const unsigned doy = doe - (365*yoe + yoe/4 - yoe/100);

	const unsigned mp = (5*doy + 2)/153;

	const unsigned m = mp + (mp < 10 ? 3 : -9);

	long yy = y + (m <= 2);

	bool isleap = yy % 4 == 0 && (yy % 100 != 0 || yy % 400 == 0);

	return (doy + 59 + isleap) % (365 + isleap) + 1;

}





struct oryza_control {

	long modelstart = 0;

	unsigned cropstart = 0;

	int max_duration = 365;

	bool water_limited = false;

	bool nitrogen_limited = false;

	double latitude = 0;

	double elevation = 0;

	double CO2 = 340;

	double ANGSTA = 0.29;

	double ANGSTB = 0.45;

	double FAOF = 1.0;

	std::string output_option;

	std::string RICETYPE = "LOWLAND";

	std::string RUNMODE = "EXPERIMENT";

	std::string ESTAB = "DIRECT-SEED";

	std::string ETMOD = "PENMAN";

	std::string PRODENV = "POTENTIAL";

	std::string NITROENV = "POTENTIAL";

	std::string WATBAL = "PADDY";

	int SBDUR = 0;

	double TMPSB = 0;

	std::vector<double> TMCTB = {0., 0., 366., 0.};

	// Soil nitrogen (experiment / control). FERTIL is a compact event list:
	// day, amount, day, amount, ... (kg N ha-1 d-1 on those DAE days; 0 elsewhere).
	// A padded FORTRAN AFGEN table (with explicit zeros) is also accepted.
	std::vector<double> FERTIL;
	std::vector<double> RECNIT = {0., 0.30, 0.2, 0.35, 0.4, 0.50, 0.8, 0.75, 1.0, 0.75, 2.5, 0.75};
	double SOILSP = 0.8;

	// Irrigation management (experiment file / control)

	int SWITIR = 0;

	double DVSIMAX = 2.0;

	double IRRI = 75.;

	int SLMIN = 3;

	double KPAMIN = 5.;

	double WCMIN = 0.30;

	int WL0DAY = 5;

	double WL0MIN = 10.;

	std::vector<double> RIRRIT;

	std::vector<double> ISTAGET;

};





struct oryza_crop {

	double RDD = 0;



	bool DLEAF = false, DROUT = false, GRAINS = false;



	int IMX = 40;

	int ILDRLV = 0, ILEFFT = 0, ILFLVT = 0, ILFSHT = 0, ILFSOT = 0, ILFSTT = 0;

	int ILKDFT = 0, ILKNFT = 0, ILSLAT = 0;

	int ILREDF = 0, ILSSGA = 0, ILTMCT = 0;



	std::vector<double> DRLVT;

	std::vector<double> EFFTB;

	std::vector<double> SLATB;

	std::vector<double> FLVTB;

	std::vector<double> FSHTB;

	std::vector<double> FSOTB;

	std::vector<double> FSTTB;

	std::vector<double> KDFTB;

	std::vector<double> KNFTB;

	std::vector<double> REDFTT;

	std::vector<double> SSGATB;

	std::vector<double> TMCTB;



	double TMAXC = 0, TMINC = 0;



	double ALAI = 0, AMAX = 0;

	double CBCHK = 0, CKCIN = 0, CKCFL = 0;

	double CO2 = 340, CO2EFF = 0, CO2LV = 0, CO2ST = 0, CO2STR = 0, CO2SO = 0, CO2REF = 340, CO2RT = 0;

	double CRGCR = 0, CRGLV = 0, CRGRT = 0, CRGSO = 0, CRGST = 0, CRGSTR = 0, CTRANS = 0;

	double DAYL = 0, DAYLP = 0, DLDR = 0, DLDRT = 0, DPAR = 0, DPARI = 0, DTGA = 0, DTR = 0;

	double DVEW = 0, DVR = 0, DVRI = 0, DVRJ = 0, DVRP = 0, DVRR = 0, DVSI = 0;

	double EFF = 0, FCLV = 0, FCRT = 0, FCSO = 0, FCST = 0, FCSTR = 0, FLV = 0, FSH = 0;

	double FRPAR = 0.5, FSO = 0, FRT = 0, FST = 0, FSTR = 0;

	double GCR = 0, GGR = 0, GLAI = 0, GLV = 0, GNGR = 0, GNSP = 0, GRT = 0, GRT1 = 0;

	double GSO = 0, GST = 0, GST1 = 0, GSTR = 0, GZRT = 0.01;

	double HU = 0, HULV = 0;

	double KEEP = 0, KDF = 0, KNF = 0, LAI = 0;

	double LAPE = 0.0001, LLV = 0, LRSTR = 0, LSTR = 0;

	double MAINLV = 0, MAINRT = 0, MAINSO = 0, MAINST = 0, MNDVS = 0, MOPP = 0;

	double NCOLD = 0, NFLV = 0, NGCR = 0, NGR = 0, NGRM2 = 0;

	double NH = 25, NPLDS = 200, NPLH = 5, NPLSB = 1000, NSP = 0, NSPM2 = 0;

	double PARCM1 = 0, PARCUM = 0, PARI1 = 0, PLTR = 0, PPSE = 0, PWRR = 0, Q10 = 2;

	double RAPCDT = 0, RDAE = 0, REDFT = 0, RGCR = 0, RGRL = 0, RMCR = 0, RTNASS = 0, RWLVG = 0;

	double RWLVG1 = 0, RWSTR = 0, RWSTR1 = 0;

	double SAI = 0, SCP = 0.2, SF1 = 1, SF2 = 1, SHCKD = 0, SHCKL = 0, SLA = 0;

	double SPFERT = 1, SPGF = 0, SSGA = 0;

	double TAV = 0, TAVD = 0, TBD = 8, TBLV = 8, TCLSTR = 10, TCOR = 0, TDRW = 0;

	double TEFF = 0, TMAX = 0, TMIN = 0, TMPCOV = 0, TMPSB = 0, TMD = 42, TNASS = 0;

	double TOD = 30, TS = 0, TSHCKD = 0, TSHCKL = 0, TSLV = 0, TREF = 25;

	double WAG = 0, WAGT = 0, WGRMX = 0, WLV = 0, WLVG = 0, WLVGI = 0, WLVGIT = 0;

	double WLVD = 0, WRR = 0, WRR14 = 0, WRT = 0, WRTI = 0, WST = 0, WSTI = 0;

	double WSO = 0, WSOI = 0, WSTS = 0, WSTR = 0;

	double ZRTI = 0.0001, ZRTM = 0, ZRTTR = 0.05, ZRTMCW = 0.25, ZRTMCD = 0.40, RGRLMX = 0, RGRLMN = 0;

	double ASLA = 0, BSLA = 0, CSLA = 0, DSLA = 0, SLAMAX = 0;

	double COLDMIN = 12, COLDEAD = 3;



	std::string SWISLA = "FUNCTION";



	double CHECKTB = 0, TSTCHK = 0;



	// Drought stress thresholds (crop.dat)

	double ULLS = 74.13, LLLS = 794.33;

	double ULDL = 630.95, LLDL = 1584.89;

	double ULLE = 1.45, LLLE = 1404.;

	double ULRT = 74.13, LLRT = 1584.89;

	std::string SWIRTR = "DATA";

	double SWIRTRF = 0.003297;



	int CROPSTA = 0, DTFSECMP = 0, EMYR = 0, EMD = 0, IDATE = 0, SBDUR = 0;

	double ANGA = 0.29, ANGB = 0.45, DAE = 0, DVS = 0, ETD = 0, EVSC = 0;

	double FAOF = 1, IR = 0, LDSTRS = 1, LESTRS = 1, LRSTRS = 1;

	double PCEW = 1, RAINCU = 0, SWR = 0, TKLT = 100, TMDA = 0, TRC = 0;

	double TRW = 0, WL0 = 0, ZRT = 0, ZRTMS = 100, LAIROL = 0, CPEW = 1;

	int NLXM = 10;

	std::vector<double> MSKPA = std::vector<double>(10, 0);

	std::vector<double> TKL = std::vector<double>(10, 0);

	std::vector<double> TRWL = std::vector<double>(10, 0);

	std::vector<double> WCAD = std::vector<double>(10, 0);

	std::vector<double> WCWP = std::vector<double>(10, 0);

	std::vector<double> WCFC = std::vector<double>(10, 0);

	std::vector<double> WCST = std::vector<double>(10, 0.3);

	std::vector<double> WCLQT = std::vector<double>(10, 0.3);



	double TNSOIL = 0, NACR = 0, NSLLV = 1, RNSTRS = 1;

	double NFLVI = 0.5;
	double FNLVI = 0.025;
	double NMAXUP = 8.;
	double NMAXSO = 0.0175;
	double RFNLV = 0.004;
	double RFNST = 0.0015;
	double TCNTRF = 10.;
	double FNTRT = 0.15;

	std::vector<double> NFLVTB, NMAXLT;
	std::vector<double> NMINLT = {0.0, 0.025, 1.0, 0.012, 2.1, 0.007, 2.5, 0.007};
	std::vector<double> NMINSOT = {0., 0.006, 50., 0.0008, 150., 0.0125, 250., 0.015, 400., 0.017, 1000., 0.017};
	std::vector<double> NSLLVT = {0., 1.0, 1.1, 1.0, 1.5, 1.4, 2.0, 1.5, 2.5, 1.5};

};





struct oryza_soil {

	std::string SCODE = "PADDY";

	int SWITPD = 0;

	int SWITGW = 1;

	int SWITPF = 1;

	int SWITVP = -1;

	int SWITKH = 1;

	int NL = 10;

	int NLPUD = 3;



	std::vector<double> TKL;       // layer thickness (m), for ORYZA1/WSTRESS

	std::vector<double> TKL_mm;    // layer thickness (mm), internal PADDY

	double ZRTMS = 1.0;

	double WL0MX = 250.;

	double WL0I = 0.;

	double WL0 = 0.;



	std::vector<double> KST;

	std::vector<double> WCST;

	std::vector<double> VGA, VGL, VGN, VGR, PN;

	std::vector<double> WCFC, WCWP, WCAD, WCLI;

	std::vector<double> WCSTRP;

	double FIXPERC = 10000.;

	std::vector<double> PERTB;

	std::vector<double> PTABLE;

	std::vector<double> ZWTB = {1., 200., 366., 200.};

	double ZWTBI = 100., MINGW = 100., MAXGW = 100.;

	double ZWA = 1.0, ZWB = 0.5;

	double PFCR = 6.0;

	std::string RIWCLI = "NO";

	std::vector<int> WCLINT;



	// Runtime water balance state

	std::vector<double> WCL, WCLQT;

	std::vector<double> MSKPA;

	std::vector<double> WL, WLFL;

	std::vector<double> WLFC, WLAD, WLST, KSAT;

	std::vector<double> CAPRI, GWFILL, WLCH, ZL, MS;

	std::vector<double> MSUC, TOTPOR, VL;



	bool PUDDLD = false;

	bool GRWAT = false;

	bool RWCLI = false;

	bool CRACKS = false;

	int IGW = 1;

	double ZW = 0., ZWPREV = 0.;

	double DSPW = 1.;

	double PERC = 0., PERCOL = 0., RUNOF = 0.;

	double EVSW = 0., EVSWS = 0.;

	double WL0CH = 0., WCUMCH = 0.;

	double CAPTOT = 0., GWTOT = 0., WL0FILL = 0.;

	double DRAIN = 0.;

	double TKLT = 0.;

	double WCUM = 0., WCUMI = 0.;

	double TIME = 0.;



	// Water balance check accumulators

	double WCUMCO = 0., WL0CO = 0.;

	double IRCUM1 = 0., RAINCUM1 = 0., RUNOFCUM1 = 0.;

	double EVSWCUM1 = 0., TRWCUM1 = 0., DRAICUM1 = 0.;

	double UPRICUM1 = 0., GWCUM1 = 0., WL0FCUM1 = 0.;

	double PERCCUM1 = 0., CAPTOTCUM1 = 0.;

};





struct Atmosphere {

	double latitude = 0;

	double TMMN = 0, TMMX = 0, VP = 0, WN = 0, RAIN = 0;

	double TMDA = 0;

};





struct Weather {

	std::vector<long> date;

	std::vector<double> tmin;

	std::vector<double> tmax;

	std::vector<double> srad;

	std::vector<double> wind;

	std::vector<double> vapr;

	std::vector<double> prec;

};





struct OryzaOutput {

	std::vector<std::string> names;

	std::vector<double> values;

};



struct OryzaSoilCollection {
	std::vector<oryza_soil> soils;
	size_t size() const { return soils.size(); }
	void push_back(oryza_soil s) { soils.push_back(std::move(s)); }
};



struct oryza_model {

	double DOY = 0;

	int IDOY = 0;

	int YEAR = 0;

	int IYEAR = 0;

	unsigned time = 0;

	unsigned step = 0;

	double DELT = 1;

	bool fatalError = false;

	int NL = 10;

	bool TERMINAL = false;

	std::vector<std::string> messages;

	OryzaOutput output;



	oryza_control control;

	oryza_crop crop;

	oryza_soil soil;

	Atmosphere atm;

	Weather wth;

	// NSOIL SAVE state (nitrogen balance)
	std::vector<double> FERTIL_TB; // AFGEN table (normalized from control.FERTIL)
	double NFERTP = 0.;
	double XFERT = 0.;

	// NCROP2 SAVE state
	double ANLV = 0., ANSO = 0., ANST = 0., ANLD = 0., ANCR = 0.;
	double ANLVA = 0., ANSTA = 0., ANCRF = 0.;
	double NALVS = 0., NASTS = 0., NASOS = 0., NACRS = 0., NTRTS = 0.;
	double NALV = 0., NAST = 0., NASO = 0.;
	double NLV = 0., NST = 0., NSO = 0., NLDLV = 0.;
	double NLVAN = 0., NSTAN = 0., NTRT = 0.;
	double FNLV = 0., FNST = 0., FNSO = 0.;
	double NMAXL = 0., NMINL = 0., NMINSO = 0.;

	void model_initialize();

	void model_rate();

	void model_state();

	void update_cropsta();



	void oryza_initialize();

	void oryza_rate();

	void oryza_state();



	void model_output();

	bool weather_step();

	void run();

	std::vector<double> run_batch(
		std::vector<double> tmin, std::vector<double> tmax, std::vector<double> srad,
		std::vector<double> prec, std::vector<double> vapr, std::vector<double> wind,
		std::vector<long> date, std::vector<long> mstart, std::vector<int> soilindex,
		OryzaSoilCollection soils, std::vector<double> depth,
		std::vector<double> elevation, std::vector<double> latitude);

	double cCO2 = 0, cKNF = 0, cNFLV = 0, cREDFT = 0;

};



#endif

