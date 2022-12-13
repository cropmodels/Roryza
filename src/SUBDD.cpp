//----------------------------------------------------------------------*
//  SUBROUTINE SUBDD                                                    *
//  Purpose: This subroutine calculates the daily amount of heat units  *
//           for calculation of the phenological development rate and   *
//           early leaf area growth.                                    *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// TMAX    R4  Daily maximum temperature (oC)                        I  *
// TMIN    R4  Daily minimum temperature (oC)                        I  *
// TBD     R4  Base temperature for development (oC)                 I  *
// TOD     R4  Optimum temperature for development (oC)              I  *
// TMD     R4  Maximum temperature for development (oC)              I  *
// HU      R4  Heat units (oCd d-1)                                  O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <algorithm>
#include <math.h>
#include "model.h"

// using namespace std;

double SUBDD( double TMAX, double TMIN, double TBD, double TOD, double TMD ){

    //-----Local parameters
    double TD, TM, TT;
    //output
    double HU;

    TM = (TMAX+TMIN)/2.;
    TT = 0.;

    for(int i = 1; i <= 24; i++){
        TD = TM + 0.5 * fabs(TMAX - TMIN) * cos(0.2618 * double(i-14));

        if ((TD > TBD)  && (TD < TMD)){
            if (TD > TOD) TD = TOD-(TD-TOD)*(TOD-TBD)/(TMD-TOD);
            TT = TT+(TD-TBD)/24.;
        }
    }
    HU = TT;
    return HU;
}
















