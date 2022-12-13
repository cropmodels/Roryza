
#include <vector>
#include <math.h>
#include <algorithm>
#include "oryzaUtil.h"
#include "model.h"

// using namespace std;

//----------------------------------------------------------------------*
// SUBROUTINE GWTAB                                                     *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// ITASK   I4  Task that subroutine should perform (-)              C,I *
// SWITGW  I4  Groundwater switch (-)                                C  *
// NL      I4  Number of soil layers (-)                             I  *
// DOY     R4  Daynumber (January 1 = 1) (-)                         I  *
// DELT    R4  Time step of integration (d)                          T  *
// WLFL    R4  Array of fluxes out off soil compartments (mm d-1)    I  *
// TKL     R4  Array of thicknesses soil layers (m)                  I  *
// ZWPREV  R4  Groundwater tabel depth of previous day (cm)          O  *
// IGW     I4  Number of shallowest soil compartment in groundwater     *
//             (-)                                                   O  *
// ZW      R4  Depth of groundwater table below soil surface (cm)    O  *
//                                                                      *
// FUNCTIONS used : LINT2                                               *
//----------------------------------------------------------------------*
/*
void GWTAB( int ITASK, int SWITGW, int NL, double DOY, double DELT, std::vector<double> WLFL, std::vector<double> TKL,
            double &ZWPREV, int &IGW, double &ZW){
    //     Local variables
    double    ZLL, ZWL;

    if (ITASK == 1)  {
        if (SWITGW == 1)  {
            ZW = AFGEN(ZWTB, DOY);
        } else {
            ZW = ZWTBI;
        }
        ZWPREV = ZW;
    }
    IGW = 1;
    ZLL = 0.;
    ZWL = 0.;
    while( (IGW <= NL)  && (ZWL <= 0.) ){
        ZLL = ZLL+TKL[IGW - 1]/10.;
        ZWL = ZLL-ZWPREV;
        IGW = IGW+1;
    }

    if (ZWL > 0.) {IGW = IGW-1;}
    if (ITASK == 3)  {
        if (SWITGW == 1)  {
            ZW = AFGEN(ZWTB,DOY);
        } else {
            ZW = ZW+ZWA*DELT-ZWB*10.*WLFL[IGW - 1]*DELT;
            if (ZW < MINGW) {ZW = MINGW;}
            if (ZW > MAXGW) {ZW = MAXGW;}
        }
    }



}
*/

//----------------------------------------------------------------------*
// SUBROUTINE BACKFL                                                    *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// I       I4  Compartment index (-)                                 I  *
// WL      R4  actual water content (mm)                             I  *
// FLIN    R4  flux into soil compartment (mm d-1)                   I  *
// FLOUT   R4  flux out off soil compartment (mm d-1)                I  *
// EVSWS   R4  Actual evaporation rate soil compartment 1 (m d-1)    I  *
// TRWL    R4  Array of actual transpiration rate/layer (mm d-1)     I  *
// WLST    R4  Array amount of water per soil compartment at            *
//             saturation (mm)                                       I  *
// DELT    R4  Time step of integration (d)                          T  *
// FLNEW   R4  Boundary flow between soil compartments recalculated     *
//             via subroutine BACKFL (mm d-1)                        O  *
// HLP     R4  Help variable (mm)                                    O  *
//                                                                      *
// SUBROUTINES and FUNCTIONS used : none                                *
//----------------------------------------------------------------------*

/*
void BACKFL( int I, double WL, double FLIN, double FLOUT, double EVSWS, double TRWL, double WLST, double DELT,
            double &FLNEW, double &HLP){

    HLP = 0.;
    if (I == 1)  {
        HLP = WL+(FLIN-FLOUT-EVSWS-TRWL)*DELT;
    } else {
        HLP = WL+(FLIN-FLOUT-TRWL)*DELT;
    }

    if (HLP > WLST)  {
        FLNEW = FLIN-(HLP-WLST)/DELT;
    } else {
        FLNEW = FLIN;
    }

}

*/

//----------------------------------------------------------------------*
// SUBROUTINE DOWNFL                                                    *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// I       I4  Compartment index (-)                                 I  *
// KSAT    R4  Saturated hydraulic conductivity (cm d-1)             I  *
// FLIN    R4  Flux into soil compartment (mm d-1)                   I  *
// TRWL    R4  Array of actual transpiration rate/layer (mm d-1)     I  *
// EVSWS   R4  Actual evaporation rate soil compartment 1 (m d-1)    I  *
// WL      R4  Actual water content (mm)                             I  *
// WLFC    R4  Array amount of water per soil compartment at 'field     *
//             capacity' (mm)                                        I  *
// DELT    R4  Time step of integration (d)                          T  *
// FLOUT   R4  Flux out off soil compartment (mm d-1)                O  *
//                                                                      *
// SUBROUTINES and FUNCTIONS used : none                                *
//----------------------------------------------------------------------*

/*
double DOWNFL( int I, double KSAT, double FLIN, double TRWL, double EVSWS, double WL, double WLFC, double DELT ){
    double FLOUT;
    if (I == 1)  {
        FLOUT = std::min(10.*KSAT,std::max(0.,FLIN-EVSWS-TRWL+(WL-WLFC)/DELT));
    } else {
        FLOUT = std::min(10.*KSAT,std::max(0.,FLIN-TRWL+(WL-WLFC)/DELT));
    }
    return FLOUT;
}
*/

//----------------------------------------------------------------------*
// SUBROUTINE SATFLX                                                    *
//                                                                      *
// Date     : March 1993, Version: 1.0                                  *
// Purpose  : SATFLX determines percolation rate as a                   *
//            function of ponded water depth, hydraulic conductivity    *
//            compacted layer and hydraulic conductivity subsoil        *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// TKL     R4  Array of thicknesses soil layers (m)                  I  *
// NLPUD   I4  Number of puddled soil compartments including plow       *
//             sole (-)                                              I  *
// WL0     R4  Amount of ponded water (mm)                           I  *
// PERC    R4  Percolation rate (cm d-1)                             O  *
//----------------------------------------------------------------------*

/*
double SATFLX( std::vector<double> TKL, int NLPUD, double WL0 ){
    //-----Local variables
    int i;
    double    A, DF, DFX, DGX , DZL , F , F1, F2, FX, GX;
    double    L, N , M  , HLP1, HLP2, HS, HSPREV = 0.0, HT, TKLTOT, TINY;
    TKLTOT = 0.;
//-----Arbitrary initial value for pressure head subsoil
    TINY   = 1.E-5;
//---- Pressure head top compartments (WL0 and TKL are in mm//)
    HS     = -10.0;
    i = 1;
    while( i <= NLPUD-1){
        TKLTOT = TKLTOT+TKL[i - 1]/10.;
        i++;
    }
    HT  = WL0/10.+TKLTOT;
//------- Thickness plow sole in cm
    DZL = TKL[NLPUD - 1]/10.;

    //problem common
    if (SWITKH == 1){
        //------- Van Genuchten parameters of unsaturated compartment 3
        N = VGN[NLPUD];
        A = VGA[NLPUD];
        L = VGL[NLPUD];
        M = 1.-1./N;

        //------- Assign arbitrary value to F, the difference function
        F = 10.*TINY;
        i = 1;
        while(fabs(F) > TINY){
            if (HS > 0)  {
//-------- Estimated pressure head out of range, reset to previous
//-------- value, divided by 2
//               WRITE (*,*) ' HS > 0, RESET'
                HS = (HSPREV/2);
            }
            //-------- Calculate Van Genuchten parameters
            HLP1 = pow( (1.+ pow( A*fabs(HS) ,N) ) , M)- pow(A*fabs(HS), N - 1.) ;
            HLP2 =  1.+pow( (A*fabs(HS)) , N);
            FX   = HLP1*HLP1;
            GX   = pow(HLP2, (M*(L+2.)));

            //---------Estimated flux through unsaturated compartment
            F2   = KST[NLPUD]*(FX/GX);

            //---------Calculate derivative of Van Genuchten equation
            //problem
            DFX  = 2.*HLP1*(M*( pow(HLP2, M-1.) )*(-N)*A*( pow((A*fabs(HS)), N-1.) ) - A*(-N+1.)*( pow(A*fabs(HS), N-2.) ));
            DGX  = M*(L+2.)*( pow(HLP2, (M*(L+2.)-1.)) )*A*(-N) * ( pow(A*fabs(HS), N-1.) );

            //---------Estimated flux through saturated compartment
            F1   = KST[NLPUD - 1]*((HT-HS)/DZL+1.);

            //---------Continue iteration until difference f is approximately zero
            F    = F2 - F1;
            DF   = KST[NLPUD]*((GX*DFX-FX*DGX)/(GX*GX))+KST[NLPUD - 1]/DZL;

            //---------Keep current value of HS in case next is positive
            HSPREV = HS;

            //---------Estimate new value for pressure head unsaturated compartment
            HS   = HS-F/DF;
            i++;
        }

    }
    else if(SWITKH == 2){
        //------- Power function for hydraulic conductivity
        N = PN[NLPUD];
        F = 10.*TINY;

        i = 1;
        while( (fabs(F) > TINY) && (i <= 50) ){
            if (HS > 0) HS = HSPREV/2;
            //-----------Estimated flux through saturated compartment
            F1 = KST[NLPUD - 1]*((HT-HS)/DZL+1.);

            //-----------Estimated flux through unsaturated compartment
            F2 = KST[NLPUD]*(pow(fabs(HS), N));
            //-----------Continue iteration until difference f is approximately zero
            F  = F2-F1;
            DF = -N*( pow( (fabs(HS)), N-1.) )*KST[NLPUD]+KST[NLPUD - 1]/DZL;
            //-----------Keep current value of HS in case next is positive
            HSPREV = HS;

            //-----------Estimated new value for pressure head unsaturated compartment
            HS = HS-F/DF;
            i++;
        }
    }
    //-----Save flux (mm/d)
    
    double PERC = 10.*(F1+F2)/2.;
    return PERC;

}

*/

//----------------------------------------------------------------------*
//  SUBROUTINE SHRINK                                                   *
//                                                                      *
//  Version: 1.0                                                        *
//  Date:    April 1994                                                 *
//  Purpose: SHRINK calculates volumetric water content of              *
//           puddled soil compartments                                  *
//           using a shrinkage factor equal to WCSTRP/WCST              *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// ITASK   I4  Task that subroutine should perform (-)               I  *
// I       I4  Compartment index (-)                                 I  *
// MNL     I4  Maximum compartment numbers (-)                       I  *
// WL      R4  Actual water content / layer (mm)                     I  *
// TKL     R4  Array of thicknesses soil layers (mm)                 I  *
// WCST    R4  Array of water content saturation / layer (cm3 cm-3)  I  *
// WCSTRP  R4  Array saturated volumetric water content ripened      I  *
//             soil per soil compartment (cm3 cm-3)                  I  *
// WCL     R4  Actual water content / layer (cm3 cm-3)               O  *
// TOTPOR  R4  Total porosity / layer (cm3 cm-3)                     O  *
// VL      R4  Thickness soil compartment after shrinkage (mm)       O  *
//                                                                      *
// SUBROUTINES and FUNCTIONS called: none                               *
//----------------------------------------------------------------------*

//void SHRINK( int ITASK, int I, int MNL, double WL, std::vector<double> TKL, std::vector<double> WCST, std::vector<double> WCSTRP,
//            double &WCL, double &TOTPOR, double &VL){
//    //-----Local parameters
//    std::vector<double> WLLOW(MNL, 0);
//    double WLST, WLSTRP;
//
//    WLST   = WCST  *TKL;
//    WLSTRP = WCSTRP*TKL;
//
//    if (ITASK == 1)  {
//        WLLOW[I - 1] = WL;
//        if (WL >= WLSTRP)  {
//            VL     = TKL;
//            TOTPOR = WL/TKL;
//            WCL    = WL/TKL;
//        } else {
//            VL     = TKL;
//            TOTPOR = WLSTRP/TKL;
//            WCL    = WL/TKL;
//        }
//    }
//    else{
//        if (WL >= WLLOW[I - 1])  {
//            VL     = TKL;
//            WCL    = WL/TKL;
//        }
//        if ((WL < WLLOW[I - 1])  && (WL >= WLSTRP))  {
//            WLLOW[I - 1] = WL;
//            VL       = TKL;
//            TOTPOR   = WL/TKL;
//            WCL      = WL/TKL;
//        }
//        if ((WL < WLLOW[I - 1])  && (WL < WLSTRP))  {
//            WLLOW[I - 1] = WL;
//            VL       = TKL;
//            TOTPOR   = WLSTRP/TKL;
//            WCL      = WL/TKL;
//        }
//    }
//}

//----------------------------------------------------------------------*
// SUBROUTINE SUMSKM2                                                   *
//                                                                      *
// Purpose: SUMSKM2 calculates the hydraulic conductivity at            *
//          given suction for compartment I on the basis of chosen      *
//          option                                                      *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// I       I4  Compartment index (-)                                 I  *
// MS      R4  Soil water suction (cm)                               I  *
// WCST    R4  Array of water content saturation / layer (cm3 cm-3)  O  *
// KMS     R4  Hydraulic conductivity (cm d-1)                       O  *
//                                                                      *
// SUBROUTINES called : SUERR, SUWCMS2                                  *
//                                                                      *
// FUNCTIONS called   : none                                            *
//                                                                      *
// FILE usage : none                                                    *
//----------------------------------------------------------------------*

//void SUMSKM2( int I, double MS, std::vector<double> WCST, double KMS ){
//    //-----Local variables
//    double    HLP1, HLP2, HLP3, MSAD, TINY, VGM, WCL, WREL;
//    TINY = 1.e-10;
//    MSAD = 1.e7;
//
//    //-----Check input value MS
//    //problem
//    //if (MS < -TINY  || MS > 1.E8)  SUERR(1,MS,0.,1.E8);
//    if (MS >= MSAD-TINY)  {
////--------Air dry
//        KMS = 0.;
//    } else{
//        if (SWITKH == 1){
//            //-----------Van Genuchten conductivity
//            WCL = 0.;
//            //           Dummy value; WCL is returned by SUWCMS2//
//            SUWCMS2(I,2,WCST,WCL,MS);
//            VGM  = 1.0-1.0/VGN[I - 1];
//            WREL = (WCL-VGR[I - 1])/(WCSTRP[I - 1]-VGR[I - 1]);
//            HLP1 = pow(WREL, VGL[I - 1]);
//            HLP2 = 1.0-pow( WREL, (1./VGM) );
//            HLP3 = 1.0-pow(HLP2, VGM);
//            KMS  = KST[I - 1]*HLP1*HLP3*HLP3;
//        }
//        else if( SWITKH == 2 ){
//            //-----------Power function conductivity
//            if (MS <= 1.) KMS = KST[I - 1];
//            if (MS > 1.) KMS = KST[I - 1]*( pow(MS, PN[I-1]) );
//        }
//        else if(SWITKH == 5){
//            //problem
//        }
//        if (KMS < TINY) KMS = 0.;
//    }
//}


//----------------------------------------------------------------------*
// SUBROUTINE SUWCHK                                                    *
//                                                                      *
// Purpose: SUWCHK checks the soil water balance by comparing           *
//          time-integrated boundary fluxes versus change in            *
//          total amount of water contained in the system.              *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// CKWFL   R4  Sum of time-integrated boundary fluxes (mm)           I  *
// CKWIN   R4  Change in water storage since start (mm)              I  *
// TIME    R4  Time of simulation (d)                                I  *
//                                                                      *
// SUBROUTINES called : none                                            *
//                                                                      *
// FUNCTIONS called   : none                                            *
//                                                                      *
// FILE usage :       (screen)                                          *
//----------------------------------------------------------------------*

//problem check



//----------------------------------------------------------------------*
//  SUBROUTINE SUWCMS2                                                  *
//                                                                      *
//  Purpose: SUWCMS2 calculates volumetric soil water content from      *
//           soil water suction, and vice versa. Various options are    *
//           offered. See SWIT8 in input file or SAWAH manual.          *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// I       I4  Compartment index (-)                                 I  *
// SWIT4   I4  Switch to set request MS(WCL) or WCL(MS) (-)          I  *
// WCST    R4  Array of water content saturation / layer (cm3 cm-3)  I  *
// WCL     R4  Array of actual water content / layer (cm3 cm-3)     I/O *
// MS      R4  Soil water suction (cm)                              I/O *
//                                                                      *
//  SUBROUTINES called:  SUERR                                          *
//                                                                      *
//  FUNCTION called:     none                                           *
//                                                                      *
//  FILE usage:          none                                           *
//                                                                      *
//***********************************************************************

//void SUWCMS2( int I, int SWIT4, std::vector<double> WCST, std::vector<double> &WCL, double &MS ){
//    //-----Local variables
//    double    HLP1, HLP2, HLP3, HLP4, TINY, VGM, WREL;
//    TINY = 0.001;
//    if (SWIT4 == 1)  {
////--------Suction calculated from water content
//        if (WCL < WCAD[I - 1]  || WCL > WCST) SUERR(3,WCL,WCAD[I - 1],WCST);
//        if (WCL > WCSTRP[I - 1])  {
////           It is assumed that MS remains zero during shrinkage
//            MS = 0.;
//        } else {
////-----------Van Genuchten option
//            HLP1 = AMAX1(WCAD[I - 1],WCL);
//            WREL = (WCL-VGR[I - 1])/(WCSTRP[I -1]-VGR[I -1]);
//            VGM  = 1.-1./VGN[I - 1];
//            HLP2 = 1./VGA[I - 1];
//            HLP3 = -1./VGM;
//            HLP4 = 1./VGN[I - 1];
//            MS   = HLP2*pow((pow(WREL, HLP3)-1.), HLP4);
//        }
//    }
//    else if(SWIT4 == 2){
//        //--------Water content calculated from suction
//        if (MS < -TINY  || MS > 1.E8)  SUERR(4,MS,0.,1.E8);
//        //--------Van Genuchten option
//        VGM  = 1.-1./VGN[I - 1];
//        HLP1 = pow((MS*VGA[I - 1]), VGN[I - 1]);
//        WREL = pow((1.+HLP1), (-VGM));
//        WCL  = WREL*(WCSTRP[I - 1]-VGR[I -1])+VGR[I - 1];
//    }
//
//}





//----------------------------------------------------------------------*
// SUBROUTINE SUBSL2                                                    *
//                                                                      *
// Author   : C. Rappoldt                                               *
//            M. Wopereis (revision March 1993)                         *
// Date     : January 1986, revised June 1990                           *
//            Slightly changed to work with VG parameters, March 1993   *
// Purpose  :                                                           *
// Chapter 15 in documentation WOFOST Version 4.1 (1988) This routine   *
// calculates the rate of capillary flow or percolation between         *
// groundwater table and root zone. The stationary flow is found by     *
// integration of dZL = K.d(MH)/(K + FLW), where Z= height above        *
// groundwater, MH= matric head, K= conductivity and FLW= chosen flow.  *
// In an iteration loop the correct flow is found.  The integration     *
// goes at most over four intervals: [0,45],[45,170], [170,330] and     *
// [330,MH-rootzone] (last one on logarithmic scale).                   *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// PF      R4  pF value soil compartment (-)                         I  *
// D       R4  Distance to grounwater table (cm)                     I  *
// I       I4  Compartment index (-)                                 I  *
// WCST    R4  Array of water content saturation / layer (cm3 cm-3)  I  *
// FLOW    R4  Capillary rise calculated by subroutine SUBSL2 (mm       *
//             d-1)                                                  O  *
//                                                                      *
// SUBROUTINES called: SUMSKM2                                          *
//                                                                      *
//----------------------------------------------------------------------*


//double SUBSL2( double PF, double D, int I, std::vector<double> WCST ){
//    //output
//    double FLOW;
//
//    //-----Local variables
//    int I1,I2,I3,IINT,IMAX;
//    double    D1, DF, ELOG10, FL, FLW, FU, KMS, LOGST4, MH, PF1, Z ;
//    ELOG10 = 2.302585;
//    LOGST4 = 2.518514;
//
//    double START[4] = {0.,45.,170.,330.};
//    double PFSTAN[9] = {0.705143,1.352183,1.601282,1.771497,2.031409,2.192880, 2.274233,2.397940,2.494110};
//    double PGAU[3] = {.1127016654,.5,.8872983346};
//    double WGAU[3] = {.2777778,.4444444,.2777778};
//
//    double DEL[4];
//    double PFGAU[12];
//    double HULP[12];
//    double CONDUC[12];
//
//    double K0;
////-----Calculation of matric head and check on small pF
//    PF1  = PF;
//    D1   = D;
//    MH   = exp(ELOG10*PF1);
//    if (PF1 <= 0.) {
//        K0 = KST[i - 1];
//        //    END CHANGES
//        //    Flow in mm/d
//        FLOW = 10*K0*(MH/D-1.);
//        return FLOW;
//    }
//    IINT  = 0;
//
//    for(int i = 0; i < 4; i++){
//        if(i <= 2) DEL[i] = std::min( START[I1+1],MH ) - START[i];
//        if (i == 3) DEL[i] = PF1-LOGST4;
//        if (DEL[i] <= 0.){
//            break;
//        }
//        IINT = IINT+1;
//    }
//    for(int i = 1; i <= IINT; i++){
//        for(int j = 1; j <= 2; j++){
//            I3 = 3*(i - 1) + j;
//            if(i == IINT) continue;
//
//
//
//        }
//
//
//    }
//
//
//
//
//}










