//----------------------------------------------------------------------*
//  SUBROUTINE PHENOL                                                   *
// Modified version 1; April, 2004. Bouman                              *
//  Purpose: This subroutine calculates the rate of phenological        *
//           development of the crop based on photoperiod and           *
//           temperature.                                               *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// DVS     R4  Development stage of the crop (-)                     I  *
// DVRJ    R4  Development rate juvenile ((oCd)-1)                   I  *
// DVRI    R4  Development rate, photoperiod-sensitive phase         I  *
//             ((oCd)-1)                                                *
// DVRP    R4  Development rate, PI phase ((oCd)-1)                  I  *
// DVRR    R4  Development rate, reproductive phase ((oCd)-1)        I  *
// HU      R4  Heat units (oCd d-1)                                  I  *
// DAYL    R4  Astronomic daylength (base = 0 degrees) (h)           I  *
// MOPP    R4  Maximum optimum photoperiod (h)                       I  *
// PPSE    R4  Photoperiod sensitivity (h-1)                         I  *
// TS      R4  Temperature sum (oCd)                                 I  *
// SHCKD   R4  Delay parameter in phenology ((oCd)(oCd)-1)           I  *
// CROPSTA I4  Crop stage (-)                                        I  *
// DVR     R4  Development rate of the crop (d-1)                    O  *
// TSHCKD  R4  Transpl. shock for phenol. development (oCd)          O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

std::vector<double> PHENOL( double DVS, double DVRJ, double DVRI, double DVRP, double DVRR, double HU, double DAYL, double MOPP,
                        double PPSE, double TS, double SHCKD, int CROPSTA){

    //-----Local parameters
    double    DL, TSTR, PPFAC;
    //output
    double DVR = 0.0, TSHCKD;

    //bb patch to re-initialize transplanting shock:
    if (DVS < 0.01) {
        TSHCKD = 0.;
    }

    if (DVS >= 0. && DVS < 0.40) {
        DVR = DVRJ*HU;
    }

    if (DVS >= 0.40  && DVS < 0.65)  {
        DL = DAYL+0.9;
        if (DL < MOPP)  {
            PPFAC = 1.;
        } else {
            PPFAC = 1.-(DL-MOPP)*PPSE;
        }
        PPFAC = std::min(1.,std::max(0.,PPFAC));
        DVR   = DVRI*HU*PPFAC;

    }

    if (DVS >= 0.65  && DVS < 1.00) {
        DVR = DVRP*HU;
    }
    if (DVS >= 1.00) {
        DVR = DVRR*HU;
    }

    if (CROPSTA  ==  3) {
        TSTR = TS;
    }
    TSHCKD = SHCKD*TSTR;

    if (CROPSTA > 3 && TS < (TSTR+TSHCKD)) {
        DVR = 0.;
    }
    return {DVR, TSHCKD};


}

