//----------------------------------------------------------------------*
// SUBROUTINE ORYZA1                                                    *
// Rice crop growth module of ORYZA2000 model                           *
//                                                                      *
// Date      : November 2002                                            *
//          Version august, 2003                                        *
// History   : Adapted from ORYZA1 (1995), and ORYZA_W (1996) models    *
//             This version for release under FSEWin                    *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// ITASK   I4  Task that subroutine should perform (-)               I  *
// IUNITD  I4  Unit that can be used for input files (-)             I  *
// IUNITL  I4  Unit number for log file messages (-)                 I  *
// FILEI1  C*  Name of file with input model data (-)                I  *
// FILEIT  C*  Name of experimental data file (-)                    I  *
// OUTPUT  L4  Flag to indicate if output should be done (-)         I  *
// TERMNL  L4  Flag to indicate if simulation is to stop (-)        I/O *
// IDOY    I4  Day number within year of simulation (d)              I  *
// DOY     R4  Day number (January 1 = 1) (d)                        I  *
// TIME    R4  Time of simulation (d)                                T  *
// DELT    R4  Time step of integration (d)                          I  *
// LAT     R4  Latitude of site (dec.degr.)                          I  *
// RDD     R4  Daily shortwave radiation (J m-2 d-1)                 I  *
// TMMN    R4  Daily minimum temperature (degrees C)                 I  *
// TMMX    R4  Daily maximum temperature (degrees C)                 I  *
// NFLV    R4  N fraction in the leaves (g N m-2)                    I  *
// NSLLV   R4  N stress factor that accelerates leaf death (-)       I  *
// RNSTRS  R4  N stress reduction factor for RGRL (-)                I  *
// ESTAB   C*  Mode of establishment (-)                             I  *
// TKLT    R4  Thickness of combined soil layers (m)                 I  *
// ZRTMS   R4  Maximum rooting depth of soil (m)                     I  *
// CROPSTA I4  Crop stage (-)                                        I  *
// LRSTRS  R4  Leaf rolling stress factor (-)                        I  *
// LDSTRS  R4  Leaf death stress factor (-)                          I  *
// LESTRS  R4  Leaf expansion stress factor (-)                      I  *
// PCEW    R4  Reduction in potential transpiration rate (-)         I  *
// DAE     R4  Days after emergence (d)                              O  *
// LAIROL  R4  Leaf area index rolled (ha ha-1)                      O  *
// ZRT     R4  Rooting depth (m)                                     O  *
// DVS     R4  Development stage of the crop (-)                     O  *
// LLV     R4  Loss rate of leaves (kg ha-1 d-1)                     O  *
// DLDR    R4  Death rate of leaves caused by drought (kg ha-1 d-1)  O  *
// WLVG    R4  Dry weight of green leaves (kg ha-1)                  O  *
// WST     R4  Dry weight of stems (kg ha-1)                         O  *
// WSO     R4  Dry weight of storage organs (kg ha-1)                O  *
// GSO     R4  Growth rate of storage organs (kg ha-1 d-1)           O  *
// GGR     R4  Rate of increase in grain weight (kg ha-1 d-1)        O  *
// GST     R4  Growth rate of stems (kg ha-1 d-1)                    O  *
// GLV     R4  Growth rate of leaves (kg ha-1 d-1)                   O  *
// PLTR    R4  Intermediate variable for planting density (-)        O  *
//                                                                      *
// SUBROUTINES called: PHENOL, SUBLAI3, SUBDD, SUBCD, SUBCBC, SUBGRN    *
//                     GPPARSET, SGPCDT                                 *
//                                                                      *
// Files included:      -                                               *
//                                                                      *
//----------------------------------------------------------------------*


#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"
#include "oryzaUtil.h"


// using namespace std;

void oryza_model::oryza_initialize() {
    crop.DVS = 0.;
    crop.PARCUM = 0.;
    crop.PARCM1 = 0.;
    crop.WLVG = 0.01;
    crop.WLVD = 0.;
    crop.WSTS = 0.01;
    crop.WSTR = 0.;
    crop.WSO = 0.;
    crop.WRT = 0.;

    crop.WST = crop.WSTS + crop.WSTR;
    crop.WLV = crop.WLVG + crop.WLVD;
    crop.WAG = crop.WLVG + crop.WST + crop.WSO;
    crop.WAGT = crop.WLV + crop.WST + crop.WSO;
    crop.TDRW = crop.WLV + crop.WST + crop.WSO + crop.WRT;
    crop.WRR    = 0.;
    crop.WRR14  = 0.;
    crop.PWRR   = 0.;
    crop.NGR    = 0.;
    crop.NSP    = 0.;
    crop.DAE    = 0.;
    crop.TS     = 0.;
    crop.TMAXC  = 0.;
    crop.TMINC  = 0.;
    crop.TSLV   = 0.;
    crop.TNASS  = 0.;
    crop.WLVGIT = 0.;
    crop.LAI    = 0.;
    crop.LAIROL = 0.;
    crop.DVR    = 0.;
    crop.TSHCKD = 0.;
    crop.TSHCKL = 0.;
    crop.NSPM2  = 0.;
    crop.NGRM2  = 0.;
    crop.HU     = 0.;
    crop.HULV   = 0.;
    crop.NCOLD  = 0.;
    crop.ZRT    = 0.;
    crop.DLDRT  = 0.;
    crop.DLDR   = 0.;
    crop.KEEP   = 0.;
    crop.DLEAF    = false;
    crop.DROUT    = false;
    crop.GRAINS   = false;

}


void oryza_model::oryza_rate() {
//--------Re-initialize weights and LAI at day of emergence
    if(crop.CROPSTA == 1){
        crop.DVS = crop.DVSI;
        crop.WLVG = crop.WLVGI;
        crop.WLVD = 0.;
        crop.WSTS = crop.WSTI;
        crop.WSTR = 0.;
        crop.WST = crop.WSTS + crop.WSTR;
        crop.WSO = crop.WSOI;
        crop.WRT = crop.WRTI;
        crop.ZRT = crop.ZRTI;
        if(control.ESTAB == "TRANSPLANT") crop.LAI = crop.LAPE * crop.NPLSB;
        if(control.ESTAB == "DIRECT-SEED") crop.LAI = crop.LAPE * crop.NPLDS;
    }

    //--------Re-initialize rooting depth at day of transplanting
    if(crop.CROPSTA == 3){
        crop.ZRT = crop.ZRTTR;
    }

    //=======SKIP ALL RATE CALCULATIONS BEFORE EMERGENCE
    if(crop.CROPSTA >= 1) {
        //----------Set DROUT when leaf expansion is reduced in the
        //          vegetative growth phase (=> extra root growth)
        if ((crop.DVS < 1) && (crop.LESTRS < 1.)) {
            crop.DROUT = true;
        } else {
            crop.DROUT = false;
        }

        //----------Computation of weather variables
        if (crop.CROPSTA <= 2) {
            crop.TMPCOV = crop.TMPSB;
        } else {
            crop.TMPCOV = 0.;
        }

        //problem
        crop.TCOR = AFGEN(crop.TMCTB, DOY);
        crop.TMAX = atm.TMMX + crop.TCOR + crop.TMPCOV;
        crop.TMIN = atm.TMMN + crop.TCOR;
        crop.TAV = (crop.TMIN + crop.TMAX) / 2.;
        crop.TAVD = (crop.TMAX + crop.TAV) / 2.;
        crop.DTR = crop.RDD;
        // RJH begin
        if (crop.TMIN < 5.) TERMINAL = true;
        // RJH end

        //----------Counter for days after emergence
        crop.RDAE = 1.;
        //----------Phenological development
        crop.HU = SUBDD(crop.TMAX, crop.TMIN, crop.TBD, crop.TOD, crop.TMD);
        crop.NCOLD = SUBCD2(crop.COLDMIN, crop.CROPSTA, crop.TAV);

        std::vector<double> phenol = PHENOL(crop.DVS, crop.DVRJ, crop.DVRI, crop.DVRP, crop.DVRR, crop.HU, crop.DAYL,
                                       crop.MOPP, crop.PPSE, crop.TS, crop.SHCKD, crop.CROPSTA);
        crop.DVR = phenol[0];
        crop.TSHCKD = phenol[1];
        //----------Effect of drought stress on development rate
        if (crop.DVS < 1.0) {
            // BB: REMOVE THIS; IT CAN TAKE MORE THAN 1 YEAR TO COMPLETE A CROP CYCLE////
            //              DVEW = LESTRS + (DVS*(1.-LESTRS))
            crop.DVEW = 1.;
        } else{
            if (crop.DVS >= 1.)  {
                crop.DVEW = 1.;
            }
            
            crop.DVR = crop.DVR*crop.DVEW;
            //----------CO2 concentration
            crop.CO2EFF = (1.-exp(-0.00305*crop.CO2-0.222)) / (1.-exp(-0.00305*crop.CO2REF-0.222));
            crop.EFF = AFGEN(crop.EFFTB, crop.TAVD)*crop.CO2EFF;
            //----------Leaf rolling under drought stress (only for photosynthesis)
            crop.LAIROL = crop.LAI*(0.5*crop.LRSTRS+0.5);

            //--------- Add specific stem area to leaf area
            crop.SSGA = AFGEN(crop.SSGATB,crop.DVS);
            crop.SAI  = crop.SSGA*crop.WST;
            crop.ALAI   = crop.LAIROL+0.5*crop.SAI;

            //----------Intercepted solar radiation
            crop.KDF   = AFGEN(crop.KDFTB,crop.DVS);
            crop.REDFT = AFGEN(crop.REDFTT,crop.TAVD);
            crop.KNF   = AFGEN(crop.KNFTB,crop.DVS);

            //----------Daily gross canopy CO2 assimilation (DTGA)
            GPPARSET (crop.CO2, crop.KNF, crop.NFLV, crop.REDFT);
            //----------The value 2 in next argument list: accuracy for Gauss integration over canopy.
            //          If value=1 => 3-points Gauss over canopy (as in TOTASP); of value = 2 =>
            //          enhanced accuracy if required as detrmined within the subroutine TRY//
            SGPCDT(1, IDOY, atm.latitude, crop.DTR, crop.FRPAR, crop.SCP, crop.AMAX, crop.EFF, crop.KDF, crop.ALAI, crop.DAYL, crop.DAYLP,
                    crop.DTGA, crop.RAPCDT);
            crop.PARI1 = crop.RAPCDT/1.e6;
            crop.DPARI = crop.RAPCDT/1.e6;
            crop.DPAR  = crop.FRPAR*crop.DTR/1.E6;

            //----------Unrolling of ALAI again
            crop.ALAI  = crop.LAI+0.5*crop.SAI;

            //----------Effect of drought stress on DTGA
            crop.DTGA  = crop.DTGA*crop.PCEW;

            //----------Relative growth rates of shoots and roots
            //          Effect of drought stress on shoot-root partitioning
            //BB: Changed according to SUCROS2
            crop.FSH = AFGEN(crop.FSHTB,crop.DVS);
            if (crop.DVS < 1.)  {
                crop.FSH  = (crop.FSH*crop.CPEW)/NOTNUL((1.+(crop.CPEW-1.)*crop.FSH));
            }
            crop.FRT = 1.- crop.FSH;

            //----------Relative growth rates of shoot organs
            crop.FLV = AFGEN(crop.FLVTB,crop.DVS);
            crop.FST = AFGEN(crop.FSTTB,crop.DVS);
            crop.FSO = AFGEN(crop.FSOTB,crop.DVS);

            //----------Check sink limitation based on yesterday's growth rates
            //          and adapt partitioning stem-storage organ accordingly
            if (crop.GRAINS)  {
                if (crop.GGR >= (crop.PWRR-crop.WRR))  {
            //RJH
            // WE HAD A CRASH WHEN NOTNUL(GCR*FSH) WAS VERY SMALL
            // WHILE std::max() WAS ZERO. THIS SHOULD RETURN ZERO BUT RETURNED INFINITY INSTEAD
                    if ( std::max(0.,(crop.PWRR-crop.WRR)) > 0)  {
                        crop.FSO = std::max(0.,(crop.PWRR-crop.WRR)/NOTNUL((crop.GCR*crop.FSH)));
                    } else { ;
                        crop.FSO = 0;
                    }
            // END RJH
                    crop.FST = 1.-crop.FSO-crop.FLV;
                }
            }

            //----------Loss rates of green leaves and stem reserves
            crop.LLV  = crop.NSLLV*crop.WLVG*AFGEN(crop.DRLVT,crop.DVS);
            crop.LSTR = INSW(crop.DVS-1.,0.,crop.WSTR/crop.TCLSTR);

            //----------Maintenance requirements
            crop.TEFF = pow(crop.Q10, ((crop.TAV-crop.TREF)/10.));
            crop.MNDVS = crop.WLVG/NOTNUL(crop.WLVG+crop.WLVD);
            crop.RMCR  = (crop.WLVG*crop.MAINLV+crop.WST*crop.MAINST+crop.WSO*crop.MAINSO+crop.WRT*crop.MAINRT)*crop.TEFF*crop.MNDVS;

            //----------Carbohydrate requirement for dry matter production (growth respiration)
            crop.CRGCR = crop.FSH*(crop.CRGLV*crop.FLV+crop.CRGST*crop.FST*(1.-crop.FSTR)+crop.CRGSTR*crop.FSTR*crop.FST+crop.CRGSO*crop.FSO)+crop.CRGRT*crop.FRT ;

            //----------Gross and net growth rate of crop (GCR, NGCR)
            crop.GCR   =((crop.DTGA*30./44.)-crop.RMCR+(crop.LSTR*crop.LRSTR*crop.FCSTR*30./12.))/NOTNUL(crop.CRGCR);
            crop.NGCR  = std::max(0.,crop.GCR-crop.LSTR*crop.LRSTR*crop.FCSTR*30./12.);

            //----------Set transplanting effect
            if (crop.CROPSTA  ==  3)  {
                crop.PLTR = crop.NPLH*crop.NH/crop.NPLSB;
            } else {
                crop.PLTR = 1.;
            }

            //----------Growth rates of crop organs at transplanting
            crop.RWLVG1 = (crop.WLVG*(1.-crop.PLTR))/DELT;
            crop.GST1   = (crop.WSTS*(1.-crop.PLTR))/DELT;
            crop.RWSTR1 = (crop.WSTR*(1.-crop.PLTR))/DELT;
            crop.GRT1   = (crop.WRT *(1.-crop.PLTR))/DELT;

            //----------Growth rates of crop organs
            crop.GRT    = crop.GCR*crop.FRT-crop.GRT1;
            crop.GLV    = crop.GCR*crop.FSH*crop.FLV-crop.RWLVG1;
            crop.RWLVG  = crop.GLV-crop.LLV;
            crop.GST    = crop.GCR*crop.FSH*crop.FST*(1.-crop.FSTR)-crop.GST1;
            crop.GSTR   = crop.GCR*crop.FSH*crop.FST*crop.FSTR-crop.RWSTR1;
            crop.RWSTR  = crop.GSTR-crop.LSTR;
            crop.GSO    = crop.GCR*crop.FSH*crop.FSO;
            if (crop.DVS > 0.95)  {
                crop.GGR = crop.GSO;
            } else { ;
                crop.GGR = 0.;
            }

//            void SUBGRN( double GCR, double CROPSTA, double LRSTRS, double DVS, double SF1, double SF2, double SPGF,
//                        double TAV, double TMAX, double NSP, double TIME, double &GNSP, double &GNGR, double &SPFERT, bool &GRAINS);

            //----------Growth rate of number of spikelets and grains
            //problem
            SUBGRN (crop.GCR,crop.CROPSTA,crop.LRSTRS,crop.DVS,crop.SF2,crop.SF1,crop.SPGF,crop.TAV,crop.TMAX,
                    crop.NSP, time, crop.GNSP,crop.GNGR,crop.SPFERT,crop.GRAINS);

            //--------- Leaf area growth (after calculation on leaf growth and loss rates//)

            //----------Temperature sum for leaf development
            crop.HULV = SUBDD(crop.TMAX,crop.TMIN,crop.TBLV,30.,42.);

            //----------Specific leaf area
            if (crop.SWISLA == "TABLE")  {
                crop.SLA  = AFGEN(crop.SLATB,crop.DVS);
            } else {
                crop.SLA = crop.ASLA + crop.BSLA*exp(crop.CSLA*(crop.DVS-crop.DSLA));
                crop.SLA = std::min(crop.SLAMAX, crop.SLA);
            }

            // BB: NEW LAI ROUTINE
//----------Leaf area index growth
            SUBLAI3(crop.CROPSTA,crop.RGRLMX,crop.RGRLMN,crop.TSLV,crop.HULV,
                    crop.SHCKL,crop.LESTRS,crop.RNSTRS,crop.SLA,crop.NH,crop.NPLH,crop.NPLSB,crop.DVS,crop.LAI,
                    control.ESTAB,crop.RWLVG,crop.DLDR,crop.WLVG,crop.GLAI,crop.RGRL);

//----------Leaf death as caused by drought stress
            crop.DLDR = 0.;
            if( crop.LDSTRS == 1. ){
                crop.DLEAF = false;
                crop.DLDRT = 0.;
            }
            if( crop.LDSTRS < 1. && !crop.DLEAF ){
                crop.WLVGIT = crop.WLVG;
                crop.DLEAF  = true;
                crop.KEEP   = crop.LDSTRS;
            }
            if(crop.DLEAF){
                if (crop.LDSTRS <= crop.KEEP)  {
                    crop.DLDR  = (crop.WLVGIT/DELT)*(1.-crop.LDSTRS)-crop.DLDRT/DELT;
                    crop.KEEP  = crop.LDSTRS;
                    crop.DLDRT = crop.DLDR*DELT+crop.DLDRT;
                }
            }

            //----------Growth respiration of the crop (RGCR)
            crop.CO2RT  = 44./12.*(crop.CRGRT *12./30.-crop.FCRT );
            crop.CO2LV  = 44./12.*(crop.CRGLV *12./30.-crop.FCLV );
            crop.CO2ST  = 44./12.*(crop.CRGST *12./30.-crop.FCST );
            crop.CO2STR = 44./12.*(crop.CRGSTR*12./30.-crop.FCSTR);
            crop.CO2SO  = 44./12.*(crop.CRGSO *12./30.-crop.FCSO );

            crop.RGCR = (crop.GRT+crop.GRT1)*crop.CO2RT + (crop.GLV+crop.RWLVG1)*crop.CO2LV +
            (crop.GST+crop.GST1)*crop.CO2ST + crop.GSO*crop.CO2SO+(crop.GSTR+crop.RWSTR1)*crop.CO2STR+
            (1.-crop.LRSTR)*crop.LSTR*crop.FCSTR*44./12.;

            crop.CTRANS = crop.RWLVG1*crop.FCLV+crop.GST1*crop.FCST+crop.RWSTR1*crop.FCSTR+crop.GRT1*crop.FCRT;
            crop.RTNASS = ((crop.DTGA*30./44.-crop.RMCR)*44./30.)-crop.RGCR-(crop.CTRANS*44./12.);

            //----------Carbon balance check
            crop.CKCIN  = (crop.WLVG+crop.WLVD-crop.WLVGI)*crop.FCLV+(crop.WSTS-crop.WSTI)*crop.FCST+crop.WSTR*crop.FCSTR
                        +(crop.WRT-crop.WRTI)*crop.FCRT+crop.WSO*crop.FCSO;
            crop.CKCFL  = crop.TNASS*(12./44.);
            //SUBCBC(<#double CKCIN#>, <#double CKCFL#>, <#double TIME#>, <#double TERMNL#>)
            
            SUBCBC(crop.CKCIN,crop.CKCFL, time ,crop.CBCHK, TERMINAL);
            
            

        }

    }
    else if(crop.CROPSTA == 0){
        crop.LAI    = 0.;
        crop.ALAI   = 0.;
        crop.LAIROL = 0;
    }

    //=======END OF SKIP WHOLE RATE CALCULATIONS BEFORE EMERGENCE

//===========Checks on simulation run
//-----------If biomass is negative: set at 0 and abort simulation

    //check error message
    if (crop.WSO < -5.  || crop.WLVG < -5.  || crop.WST < -5.  || crop.WRR < -5){
        if (crop.WRR < 0.) crop.WRR = 0.;
        TERMINAL = true;
    }

    //-----------If LAI is negative: set at 0 and abort simulation
    if(crop.LAI < -0.01){
        //error message
        TERMINAL = true;
    }

    //-----------The following only in main field
    if (crop.CROPSTA  >=  4) {
        if (crop.LDSTRS <= 0.) {
            //check error message
            TERMINAL = true;
        }
    }

}


void oryza_model::oryza_state() {
    //=======SKIP WHOLE STATE UPDATE BEFORE EMERGENCE
    if(crop.CROPSTA >= 1){
        //-----------Integrate rate variables
        crop.PARCUM = crop.PARCUM + crop.DPARI*DELT;
        crop.PARCM1 = crop.PARCM1 + crop.PARI1*DELT;
        crop.TS = crop.TS + crop.HU*DELT;
        crop.TSLV = crop.TSLV + crop.HULV*DELT;
        crop.TMAXC = crop.TMAXC + crop.TMAX*DELT;
        crop.TMINC = crop.TMINC + crop.TMIN*DELT;
        crop.DVS = crop.DVS + crop.DVR*DELT;

        crop.WLVG = crop.WLVG + (crop.RWLVG - crop.DLDR)*DELT;
        crop.WLVD = crop.WLVD + (crop.LLV + crop.DLDR)*DELT;
        crop.WSTS = crop.WSTS + crop.GST*DELT;
        crop.WSTR = crop.WSTR + crop.RWSTR*DELT;
        crop.WSO = crop.WSO + crop.GSO*DELT;
        crop.WRT = crop.WRT + crop.GRT*DELT;
        //problem
        //crop.WRR = crop.WRR + crop.GRR*DELT;
        if (crop.DVS > 0.95) {
			crop.WRR +=  crop.GSO * DELT;
		}
		
		
        crop.NGR = crop.NGR + crop.GNGR*DELT;
        crop.NSP = crop.NSP + crop.GNSP*DELT;
        crop.DAE = crop.DAE + crop.RDAE*DELT;
        crop.TNASS = crop.TNASS + crop.RTNASS*DELT;

        //-----------Calculate sums of states
        crop.WST    = crop.WSTS + crop.WSTR;
        crop.WLV    = crop.WLVG + crop.WLVD;
        crop.WAG    = crop.WLVG + crop.WST  + crop.WSO;
        crop.WAGT   = crop.WLV  + crop.WST  + crop.WSO;
        crop.TDRW   = crop.WLV  + crop.WST  + crop.WSO + crop.WRT;

        crop.PWRR   = crop.NGR*crop.WGRMX;
        crop.NGRM2  = crop.NGR/10000.;
        crop.NSPM2  = crop.NSP/10000.;

        //-----------Weight rough rice with 14% moisture
        crop.WRR14  = (crop.WRR/0.86);
        //-----------Leaf area index and total area index (leaves + stems)
        //problem
        crop.LAI = crop.LAI + crop.GLAI*DELT;

        crop.ALAI   = crop.LAI+0.5*crop.SAI;

        //-----------Root length
        if( !crop.DROUT && crop.ZRT <= crop.ZRTMCW ){
            crop.ZRTM = std::min( std::min(crop.ZRTMCW,crop.ZRTMS),crop.TKLT);
        }
        else if( !crop.DROUT && crop.ZRT > crop.ZRTMCW ){
            crop.ZRTM = std::min(std::min(crop.ZRT,crop.ZRTMS),crop.TKLT);
        }
        else if(crop.DROUT){
            crop.ZRTM = std::min(std::min(crop.ZRTMCD,crop.ZRTMS),crop.TKLT);
        }
        crop.ZRT = crop.ZRT + crop.GZRT*DELT;
        crop.ZRT = std::min(crop.ZRT, crop.ZRTM);

        //-----------Terminate simulation settings
        if(crop.DVS > 2.){
            TERMINAL = true;
            //error message
        }

        // BB 2006: Crop growth stops below certain lower threshold T days
        if(crop.NCOLD > crop.COLDEAD){
            TERMINAL = true;
        }
    }

}



























