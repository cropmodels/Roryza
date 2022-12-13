//----------------------------------------------------------------------*
//  SUBROUTINE SUBLAI2                                                  *
//  Version 2: January 2001                                             *
//                                                                      *
//  Purpose: This subroutine calculates the rate of growth of LAI of    *
//    of the crop in the seedbed and after transplanting in the field.  *
//    Reductions by N-stress and water-stress are taken into account.   *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// CROPSTA I4  Crop stage (-)                                        I  *
// RWLVG   R4  Green leaves growth rate (kg d-1 ha-1)                I  *
// DLDR    R4  Death rate green leaves (kg d-1 ha-1)                 I  *
// TSLV    R4  Temperature sum for leaf development (oC)             I  *
// HULV    R4  Daily temperature for leaf development (oC)           I  *
// SHCKL   R4  Delay parameter in development ((oCd)(oCd)-1)         I  *
// LESTRS  R4  Reduction factor for leaf elongation (-)              I  *
// SLA     R4  Specific leaf area (ha kg-1)                          I  *
// NH      R4  Number of hills (hills m-2)                           I  *
// NPLH    I4  Number of plants per hill (pl/hill)                   I  *
// NPLSB   R4  Number of plants in seedbed (pl/m2)                   I  *
// DVS     R4  Development stage of the crop (-)                     I  *
// LAI     R4  Leaf area index (ha ha-1)                             I  *
// RGRLMX  R4  Maximum relative growth rate leaves ((oCd)-1)         I  *
// RGRLMN  R4  Minimum relative growth rate leaves ((oCd)-1)         I  *
// ESTAB   C*  Establishment method (-)                              I  *
// GLAI    R4  Growth rate leaf area index (ha d-1 ha-1)             O  *
// RGRL    R4  Actual relative growth rate leaves ((oCd)-1)          O  *
//----------------------------------------------------------------------*

#include <algorithm>
#include <math.h>
#include <string>
#include "oryzaUtil.h"
#include "model.h"

// using namespace std;

void SUBLAI2(int CROPSTA, double RGRLMX, double RGRLMN, double TSLV, double HULV, double SHCKL, double LESTRS, double RNSTRS, double SLA,
            double NH, double NPLH, double NPLSB, double DVS, double LAI, std::string ESTAB, double RWLVG, double DLDR, double WLVG, double &GLAI, double &RGRL){

    //-----Local parameters
    double          TSLVTR, TSHCKL, GLAI1,GLAI2, X, TESTSET;
    double          WLVGEXP, LAIEXP, WLVGEXS, LAIEXS, TEST;
    bool       TESTL ;


    //SAVE;
    if (CROPSTA  <=  1)  {
        X       = 1.;
        TESTL   = false;
        TESTSET = 0.00001;
    }

    //===================================================================*
    //------Transplanted rice                                            *
    //===================================================================*
    //      Calculate RGRL as function of N stress limitation
    RGRL = RGRLMX - (1.-RNSTRS)*(RGRLMX-RGRLMN);

    if (ESTAB  ==  "TRANSPLANT"){
        //------- 1. Seed-bed; no drought stress effects in seed-bed//
        if (CROPSTA  <  3){
            if (LAI < 1.){
                GLAI    = LAI*RGRL*HULV;
                WLVGEXS = WLVG;
                LAIEXS  = LAI;
            } else{
                TEST = fabs((LAI/NOTNUL(WLVG))-SLA)/SLA;
                if(!TESTL){
                    if (TEST  <  TESTSET) {TESTL = true;}
                }
                if(TESTL){
                    GLAI = ((WLVG+RWLVG)*SLA)-LAI;
                } else {
                    GLAI1 = ((WLVG+RWLVG-WLVGEXS)*SLA+LAIEXS)-LAI;
                    GLAI2 = ((WLVG+RWLVG)*SLA)-LAI;
                    GLAI  = (GLAI1+X*GLAI2)/(X+1.);
                    X     = X+1.;
                }
            }
//------- 2. Transplanting effects: dilution and shock-setting
        }
        else if(CROPSTA == 3){
            TSLVTR = TSLV;
            TSHCKL = SHCKL*TSLVTR;
            GLAI   = (LAI*NH*NPLH/NPLSB) - LAI;
            TESTL  = false;
            X      = 1.;
//--------3. After transplanting: main crop growth
        }
        else if(CROPSTA == 4){
            if (TSLV < (TSLVTR+TSHCKL))  {
                GLAI = 0.;
//--------3.2. After transplanting shock; drought stress effects
            } else{
                if ((LAI < 1.0)  && (DVS < 1.0))  {
                    GLAI = LESTRS * LAI*RGRL*HULV;
                    WLVGEXP = WLVG;
                    LAIEXP  = LAI;
                } else{
                    //                 There is a transition from RGRL to SLA determined growth
                    //                 when difference between simulated and imposed SLA is less than 10%
                    TEST = fabs((LAI/NOTNUL(WLVG))-SLA)/SLA;
                    if(!TESTL){
                        if (TEST  <  TESTSET) {TESTL = true;}
                    }
                    if (TESTL)  {
                        GLAI = ((WLVG+RWLVG-DLDR)*SLA)-LAI;
                    } else{
                        GLAI1 = ((WLVG+RWLVG-DLDR-WLVGEXP)*SLA+LAIEXP)-LAI;
                        GLAI2 = ((WLVG+RWLVG-DLDR)*SLA)-LAI;
                        GLAI  = (GLAI1+X*GLAI2)/(X+1.);
                        X     = X+1.;
                    }
                }
            }
        }
//===================================================================*
//------Direct-seeded rice                                           *
//===================================================================*
    }
    else if(ESTAB  ==  "DIRECT-SEED"){
        if ((LAI < 1.0)  && (DVS < 1.0))  {
            GLAI    = LAI*RGRL*HULV * LESTRS;
            WLVGEXP = WLVG;
            LAIEXP  = LAI;
        } else{
            //           There is a transition from RGRL to SLA determined growth
            //           when difference between simulated and imposed SLA is less than 10%
            TEST = fabs((LAI/NOTNUL(WLVG))-SLA)/SLA;
            if(!TESTL){
                if (TEST  <  TESTSET) {TESTL = true;}
            }
            if(TESTL){
                GLAI = ((WLVG+RWLVG-DLDR)*SLA)-LAI;
            } else{
                GLAI1 = ((WLVG+RWLVG-DLDR-WLVGEXP)*SLA+LAIEXP)-LAI;
                GLAI2 = ((WLVG+RWLVG-DLDR)*SLA)-LAI;
                GLAI  = (GLAI1+X*GLAI2)/(X+1.);
                X     = X+1.;
            }
        }

    }

}











