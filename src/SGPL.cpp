//----------------------------------------------------------------------*
// SUBROUTINE SGPL                                                      *
// Authors: Daniel van Kraalingen                                       *
// Date   : 15-Sept-1994, Version: 1.0                                  *
// Purpose: This subroutine calculates assimilation at a single depth   *
//          in the canopy                                               *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// CSLV    R4  Scattering coefficient of leaves for visible   -      I  *
//             radiation (PAR)                                          *
// AMAX1   R4  Assimilation rate at light saturation       kg CO2/   I  *
//                                                        ha leaf/h     *
// EFF1    R4  Initial light use efficiency               kg CO2/J/  I  *
//                                                       (ha/h/m2/s)    *
// ECPDF   R4  Extinction coefficient for diffuse light              I  *
// GAI     R4  Total green area                             ha/ha    I  *
// GAID    R4  Green area index above selected height       ha/ha    I  *
// SINB    R4  Sine of solar height                           -      I  *
// RDPDR   R4  Instantaneous flux of direct photo-          W/m2     O  *
//             synthetically active radiation (PAR)                     *
// RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2     O  *
//             synthetically active irradiation (PAR)                   *
// GPL     R4  Instantaneous assimilation rate of          kg CO2/   O  *
//             leaves at depth GAI                        ha leaf/h     *
// RAPL    R4  Absorbed radiation at depth GAI            W/m2 leaf  O  *
//                                                                      *
// Fatal error checks: none                                             *
// Warnings          : none                                             *
// Subprograms called: SRDPRF                                           *
// File usage        : none                                             *
//----------------------------------------------------------------------*


#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

void SGPL( double CSLV, double AMAX1, double EFF1, double ECPDF, double GAI, double GAID, double SINB, double &RDPDR, double &RDPDF,
            double &GPL, double &RAPL){
    //     Local parameters
    double AMAX2, EFF2;

    //     Miscellaneous
    double RAPSHL, RAPPPL, RAPSLL, FSLLA;
    double GPSHL , GPSLL , TMPR1;

    //     Gauss weights for three point Gauss
    double GSX[3] = {0.112702, 0.500000, 0.887298};
    double GSW[3] = {0.277778, 0.444444, 0.277778};

    std::vector<double> srdprf = SRDPRF( GAID,CSLV,SINB,ECPDF,RDPDR,RDPDF );
    RAPSHL = srdprf[0];
    RAPPPL = srdprf[1];
    FSLLA = srdprf[2];

    //     Get photosynthesis parameters from user routine
    GPPARGET(GAI,GAID,AMAX1,EFF1,AMAX2,EFF2);

    
    //     Assimilation of shaded leaf area
    // here
    if (AMAX2 > 0.)  {
        GPSHL = AMAX2 * (1.- exp(-RAPSHL*EFF2/AMAX2));
    } else {
        GPSHL = 0.;
    }

    //     Assimilation of sunlit leaf area
    GPSLL  = 0.;
    RAPSLL = 0.;
    for(int i = 0; i < 3; i++){
        TMPR1 = RAPSHL + RAPPPL * GSX[i];
        if (AMAX2 > 0.)  {
            GPSLL = GPSLL + AMAX2 * (1.-exp(-TMPR1*EFF2/AMAX2)) * GSW[i];
        } else {
            GPSLL = 0.;
        }
        RAPSLL = RAPSLL + TMPR1 * GSW[i];
    }
    //     Local assimilation rate (GPL) and rate of
    //     absorption of PAR by canopy (RAPL)
    GPL  = FSLLA * GPSLL  + (1.-FSLLA) * GPSHL;
    RAPL = FSLLA * RAPSLL + (1.-FSLLA) * RAPSHL;


}