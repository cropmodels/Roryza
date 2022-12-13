#include <algorithm>
#include <vector>
#include <math.h>
#include <string>

using  namespace std;

void ET(int ITASK, double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, double IDOY, std::string ETMOD,
                  double CROPSTA, double NL, double FAOF, double WL0, std::vector<double> WCLQT, std::vector<double> WCST, double LAI,
                  double &EVSC, double &ETD, double &TRC);
void ET2_initialize(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
                    int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC);
void ET2_rate(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
              int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC);
void ET2_state(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
               int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC);
void GPPARGET(double xGAI, double xGAID, double xAmaxIn, double xEffIn, double &xAmaxOut, double &xEffOut);
void GPPARSET( double xCO2, double xKNF, double xNFLV, double xREDFT );

void NCROP2_initialization( double DELT, double TIME, bool TERMINAL, double DVS, double LLV, double DLDR, double WLVG, double WST,
                           double WSO, double GSO, double GST, double GLV, double PLTR, double LAI, double SLA, int CROPSTA, double TNSOIL,
                           double &NACR, double &NFLV, double &NSLLV, double &RNSTRS);
void SUBNBC( double CHKIN, double CKCFL, double TIME, double &NBCHK, bool &TERMINAL );
void NNOSTRESS2_initialization( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS );
void NNOSTRESS2_rate( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS );
double NSOIL( int ITASK, int IUNITD, int IUNITL, std::string FILEIT, double OUTPUT, double DELT, double DAE, double DVS, double NACR );
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
double SUBCD2(int COLDMIN, int CROPSTA, double TAV);
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
void WSTRESS_initialization(double DELT, double TRC, double ZRT, std::vector<double> TKL, int NL, int CROPSTA, std::vector<double> WCLQT,
                            std::vector<double> WCWP, std::vector<double> WCAD, std::vector<double> MSKPA, double &TRW, std::vector<double> &TRWL,
                            double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW);
void WSTRESS_rate(double DELT, double TRC, double ZRT, std::vector<double> TKL, int NL, int CROPSTA, std::vector<double> WCLQT,
                  std::vector<double> WCWP, std::vector<double> WCAD, std::vector<double> MSKPA, double &TRW, std::vector<double> &TRWL,
                  double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, std::string SWIRTR, double &CPEW);


//std::vector<double> SUBCBC( double CKCIN, double CKCFL, double TIME );
void SUBCBC( double CKCIN, double CKCFL, double TIME, double &CBCHK, bool &TERMNL );


struct oryza_control{

    std::string RICETYPE;
    std::string RUNMODE;
    std::string ESTAB;
    std::string ETMOD;
    std::string PRODENV;
    std::string NITROENV;

};

struct oryza_crop{
    double RDD;

//-----Local variables
    bool DLEAF, DROUT, INQOBS  , GRAINS;

    int IMX  = 40;
    int ILDRLV, ILEFFT, ILFLVT, ILFSHT, ILFSOT, ILFSTT;
    int ILKDFT, ILKNFT, ILSLAT;
    int ILREDF, ILSSGA, ILTMCT;
    
    
    
    std::vector<double> DRLVT = std::vector<double> (40, 0);
    std::vector<double> EFFTB = std::vector<double> (40, 0);
    std::vector<double> SLATB = std::vector<double> (40, 0);
    std::vector<double> FLVTB = std::vector<double> (40, 0);
    std::vector<double> FSHTB = std::vector<double> (40, 0);
    std::vector<double> FSOTB = std::vector<double> (40, 0);
    std::vector<double> FSTTB = std::vector<double> (40, 0);
    std::vector<double> KDFTB = std::vector<double> (40, 0);
    std::vector<double> KNFTB = std::vector<double> (40, 0);
    std::vector<double> REDFTT = std::vector<double> (40, 0);
    std::vector<double> SSGATB = std::vector<double> (40, 0);
    std::vector<double> TMCTB = std::vector<double> (40, 0);

    double TMAXC, TMINC;

    double    ALAI  , AMAX;
    double    CBCHK , CKCIN , CKCFL;
    double    CO2   , CO2EFF, CO2LV, CO2ST, CO2STR, CO2SO , CO2REF, CO2RT;
    double    CRGCR , CRGLV , CRGRT, CRGSO, CRGST , CRGSTR, CTRANS ;
    double    DAYL  , DAYLP , DLDR , DLDRT, DPAR  , DPARI , DTGA  , DTR;
    double    DVEW  , DVR   , DVRI , DVRJ , DVRP  , DVRR  , DVSI ;
    double    EFF   , FCLV  , FCRT , FCSO , FCST  , FCSTR , FLV   , FSH;
    double    FRPAR , FSO   , FRT  , FST  , FSTR;
    double    GCR   , GGR   , GLAI , GLV  , GNGR  , GNSP  , GRT   , GRT1;
    double    GSO   , GST   , GST1 , GSTR , GZRT;
    double    HU    , HULV ;
    double    KEEP  , KDF   , KNF     , LAI;
    double    LAPE  , LLV   , LRSTR   , LSTR;
    double    MAINLV, MAINRT, MAINSO  , MAINST, MNDVS , MOPP;
    double    NCOLD , NFLV  , NGCR  , NGR   , NGRM2;
    double    NH    , NPLDS , NPLH    , NPLSB , NSP   , NSPM2;
    double    PARCM1, PARCUM, PARI1   , PLTR  , PPSE  , PWRR  , Q10;
    double    RAPCDT, RDAE  , REDFT , RGCR    , RGRL  , RMCR  , RTNASS, RWLVG ;
    double    RWLVG1, RWSTR , RWSTR1;
    double    SAI   , SCP   , SF1     , SF2   , SHCKD , SHCKL , SLA  ;
    double    SPFERT, SPGF  , SSGA ;
    double    TAV   , TAVD  , TBD     , TBLV  , TCLSTR, TCOR  , TDRW;
    double    TEFF  , TMAX  , TMIN    , TMPCOV, TMPSB , TMD   , TNASS ;
    double    TOD   , TS    , TSHCKD  , TSHCKL, TSLV  , TREF  ;
    double    WAG   , WAGT  , WGRMX   , WLV   , WLVG  , WLVGI , WLVGIT;
    double    WLVD  , WRR   , WRR14   , WRT   , WRTI  , WST   , WSTI  ;
    double    WSO   , WSOI  , WSTS    , WSTR   ;
    double    ZRTI  , ZRTM  , ZRTTR   , ZRTMCW, ZRTMCD, RGRLMX, RGRLMN;
    double    ASLA  , BSLA  , CSLA    , DSLA  , SLAMAX;
    double    COLDMIN,COLDEAD;

    std::string SWISLA;

    double CHECKTB, TSTCHK ;


    //in model
    int CROPSTA, DTFSECMP, EMYR   , EMD    , IDATE , SBDUR;
    double ANGA   , ANGB    , DAE    , DVS    , ETD   , EVSC;
    double FAOF   , IR    , LDSTRS , LESTRS, LRSTRS;
    double PCEW   , RAINCU  , SWR    , TKLT   , TMDA  , TRC;
    double TRW    , WL0     , ZRT    , ZRTMS, LAIROL , CPEW;
    int NLXM = 10;
    std::vector<double> MSKPA = std::vector<double>(10, 0);
    std::vector<double> TKL = std::vector<double>(10, 0);
    std::vector<double> TRWL = std::vector<double>(10, 0);
    std::vector<double> WCAD = std::vector<double>(10, 0);
    std::vector<double> WCWP = std::vector<double>(10, 0);
    std::vector<double> WCFC = std::vector<double>(10, 0);
    std::vector<double> WCST = std::vector<double>(10, 0);
    std::vector<double> WCLQT = std::vector<double>(10, 0);

    //double LLV, DLDR, WLVG, WST, WSO, GSO, GGR, GST, GLV, PLTR;
    double TNSOIL, NACR, NSLLV, RNSTRS;

    //NNOSTRESS
    double       NFLVI;
    std::vector<double> NFLVTB, NMAXLT;




};

struct oryza_soil{


};

struct Atmosphere{
    double latitude;
    double TMMN, TMMX, VP, WN, RAIN;
    double TMDA;


};

struct Weather{
    std::vector<double> tmin;
    std::vector<double> tmax;
    std::vector<double> srad;
    std::vector<double> wind;
    std::vector<double> vapr;
    std::vector<double> prec;
    

};


struct oryza_model{
    double DOY;
    int IDOY;
    int YEAR;
    int IYEAR;
    int time;
    int step;
    double DELT;
    //string RUNMODE;
    //string ESTAB  , ETMOD   , PRODENV, NITROENV, WATBAL , RICETYPE;
    bool fatal_error;
    int NL;
    bool TERMINAL;
    std::vector<std::vector<double>> out;


    oryza_control control;
    oryza_crop crop;
    oryza_soil soil;
    Atmosphere atm;
    Weather wth;

    void model_initialize();
    void model_rate();
    void model_state();

    void oryza_initialize();
    void oryza_rate();
    void oryza_state();


    void model_output();
    void weather_step();
    void model_run();



//COMMON_GP.inc
 double cCO2, cKNF, cNFLV, cREDFT;
//COMMON_GWT.inc
 double ZWA,ZWB,MAXGW,MINGW,ZWTBI;
 int IZWTB;
 std::vector<double> ZWTB = std::vector<double>(400, 0);
//COMMON_HYDCON.INC
 std::vector<double> KST = std::vector<double>(10, 0);
 std::vector<double> WCAD = std::vector<double>(10, 0);
 std::vector<double> WCSTRP = std::vector<double>(10, 0);
//COMMON_NUCHT.INC
 std::vector<double> VGA = std::vector<double>(10, 0);
 std::vector<double> VGL = std::vector<double>(10, 0);
 std::vector<double> VGN = std::vector<double>(10, 0);
 std::vector<double> VGR = std::vector<double>(10, 0);
//COMMON_POWER.INC
 std::vector<double> PN = std::vector<double>(10, 0);
//COMMON_SWIT.INC
 int SWITKH;



};


