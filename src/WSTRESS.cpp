//----------------------------------------------------------------------//
//  SUBROUTINE WSTRESS                                                  //
//  Used in ORYZA2000 model version 1.0                                 //
//  Date   : December 2001                                              //
//  Author : B.A.M. Bouman                                              //
//  Version: 1.0 (Based on earlier versions of DSTRES)                  //
//  Version: Redistribution of transpiration uptake per soil layer.     //
//           Avarage stress factors over all soil layers                //
//           Calculations based on Wopereis et al (1996a)               //
//                                                                      //
//  Purpose: Calculate actual transpiration of a crop, and the  effects //
//          of water stress on growth and development of rice.          //
//                                                                      //
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      //
// name   type meaning (unit)                                     class //
// ----   ---- ---------------                                    ----- //
// ITASK   I4  Task that subroutine should perform (-)               I  //
// DELT    R4  Time step of integration (d)                          T  //
// OUTPUT  R4  Flag to indicate if output should be done (-)         I  //
// IUNITD  I4  Unit that can be used for input files (-)             I  //
// IUNITL  I4  Unit number for log file messages (-)                 I  //
// FILEI1  C*  Name of file with input model data (-)                I  //
// TRC     R4  Potential transpiration rate (mm d-1)                 I  //
// ZRT     R4  Rooting depth (m)                                     I  //
// TKL     R4  Array of thicknesses soil layers (m)                  I  //
// NL      I4  Number of soil layers (-)                             I  //
// CROPSTA I4  Crop stage (-)                                        I  //
// WCLQT   R4  Array of actual soil water contents/layer (m3 m-3)    I  //
// WCWP    R4  Array of water content at wilting point/layer (m3 m-3)I  //
// WCAD    R4  Array of water content air dry/ layer (m3 m-3)        I  //
// MSKPA   R4  Array with soil water potential/layer (KPa)           I  //
// TRW     R4  Actual transpiration rate (mm)                        O  //
// TRWL    R4  Array of actual transpiration rate/layer (mm d-1)     O  //
// LRSTRS  R4  Stress factor for rolling of leaves (-)               O  //
// LDSTRS  R4  Stress factor for dead leaves (-)                     O  //
// LESTRS  R4  Stress factor for expansion of leaves (-)             O  //
// PCEW    R4  Stress factor for CO2 assimilation (-)                O  //
//                                                                      //
// Files included: -                                                    //
//                                                                      //
//----------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <math.h>
#include "oryzaUtil.h"
#include "model.h"
#include <string>
// using namespace std;

void WSTRESS_initialization(double DELT, double TRC, double ZRT, std::vector<double> TKL, int NL, int CROPSTA, std::vector<double> WCLQT,
                            std::vector<double> WCWP, std::vector<double> WCAD, std::vector<double> MSKPA, double &TRW, std::vector<double> &TRWL,
                            double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, double &CPEW){



    double LRAV   = 1.;
    double LDAV   = 1.;
    double LEAV   = 1.;
    LESTRS = 1.;
    PCEW   = 1.;
    CPEW   = 1.;
    LRSTRS = 1.;
    LDSTRS = 1.;
    for(int i = 0; i < NL; i++){
        TRWL[i] = 0.;
    }
    TRW = 0;
}


void WSTRESS_rate(double DELT, double TRC, double ZRT, std::vector<double> TKL, int NL, int CROPSTA, std::vector<double> WCLQT,
                  std::vector<double> WCWP, std::vector<double> WCAD, std::vector<double> MSKPA, double &TRW, std::vector<double> &TRWL,
                  double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW, std::string SWIRTR, double &CPEW){

    //local variables
    double       TRRM, ZRTL, ZLL , LRAV, LEAV, LDAV;
    double       ULLS, LLLS, ULDL, LLDL, LLLE, ULLE, LLRT, ULRT;

    int INL = 10;
    //std::vector<double> TRWL(10, 0.);
    TRWL = std::vector<double> (10, 0);
    std::vector<double> LR(10, 0.);
    std::vector<double> LE(10, 0.);
    std::vector<double> LD(10, 0.);
    std::vector<double> TRR(10, 0.);
    std::vector<double> WLA(10, 0.);
    double TINY = 0.0000000001;


    //--- Bas, 8 Sept, 2006. We calculate stress also in seedbed//
    //-------Only stress in main field after day of emergence
    if( CROPSTA > 1 ){
        TRRM = TRC/(ZRT+1.0E-10);

        TRW  = 0.;
        ZLL  = 0.;
        LRAV = 0.;
        LEAV = 0.;
        LDAV = 0.;

        for(int i = 0; i < NL; i++){
            //-------------Root length in each soil layer
            ZRTL  = std::min(TKL[i],std::max((ZRT-ZLL),0.0));

            //-------------Leaf-rolling factor
            LR[i] = (log10(MSKPA[i] + TINY) - log10(LLLS)) / (log10(ULLS)-log10(LLLS));
            LR[i] = LIMIT(0.,1.,LR[i]);
            LRAV  = LRAV+(ZRTL/(ZRT+TINY))*LR[i];

            //-------------Relative leaf expansion rate factor
            LE[i] = (log10(MSKPA[i]+TINY)-log10(LLLE)) / (log10(ULLE)-log10(LLLE));
            LE[i] = LIMIT(0.,1.,LE[i]);
            LEAV  = LEAV+(ZRTL/(ZRT+TINY))*LE[i];

            //-------------Relative death rate factor
            LD[i] = (log10(MSKPA[i]+TINY)-log10(LLDL)) / (log10(ULDL)-log10(LLDL));
            LD[i] = LIMIT(0.,1.,LD[i]);
            LDAV  = LDAV+(ZRTL/(ZRT+TINY))*LD[i];

            //-------------Relative transpiration ratio (actual/potential)
            if( MSKPA[i] > 10000. ) {
                TRR[i] = 0.;
            } else{
                if(SWIRTR == "DATA"){
                    TRR[i] = (log10(MSKPA[i]+TINY)-log10(LLRT)) / (log10(ULRT)-log10(LLRT));
                    TRR[i] = LIMIT(0.,1.,TRR[i]);
                } else{
                    TRR[i]  = 2./(1.+exp(0.003297*MSKPA[i]));
                }
            }
            TRR[i]  = LIMIT(0.,1.,TRR[i]);
            WLA[i]  = std::max(0.0,(WCLQT[i]-WCWP[i])*ZRTL*1000.);
            TRWL[i] = std::min(TRR[i]*ZRTL*TRRM,WLA[i]/DELT);
            TRW     = TRW + TRWL[i];
            ZLL     = ZLL+TKL[i];
        }

        //-----Compensation of water extraction from soil layers if drought stress occurs
//     Take water from soil layer that has a surplus, starting from top.
        for(int i = 0; i < NL; i++){
            if(TRW < TRC) {
                if (TRR[i] >= 1 && TRWL[i] < WLA[i] / DELT) {
                    TRWL[i] = std::min(WLA[i] / DELT, (TRWL[i] + (TRC - TRW) / DELT));
                }
                TRW = 0.;
                for (int j = 0; j < NL; j++) {
                    TRW = TRW + TRWL[j];
                }
            }
        }

        PCEW   = NOTNUL(TRW/TRC);
        //-----Set stress factors as average over all layers: DEFAULT

        LRSTRS = LRAV;
        LDSTRS = LDAV;
        LESTRS = LEAV;

    } else{
        //-------If crop is not in the main field, set all stres factors at 1.
        PCEW   = 1.;
        LRSTRS = 1.;
        LDSTRS = 1.;
        LESTRS = 1.;

    }
    CPEW = LESTRS;

}









