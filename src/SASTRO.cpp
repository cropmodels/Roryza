//----------------------------------------------------------------------*
// SUBROUTINE SASTRO                                                    *
// Authors: Daniel van Kraalingen                                       *
// Date   : 12-June-1996, Version: 1.1                                  *
// Purpose: This subroutine calculates solar constant, daily            *
//          extraterrestrial radiation, daylength and some intermediate *
//          variables required by other routines. The routine has been  *
//          written such that latitudes from pole to pole can be used.  *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// IDOY    I4  Day of year (Jan 1st = 1)                      d      I  *
// LAT     R4  Latitude of the site                       degrees    I  *
// SOLCON  R4  Solar constant at day=IDOY                   W/m2     O  *
// ANGOT   R4  Daily extraterrestrial radiation            J/m2/d    O  *
// DAYL    R4  Astronomical daylength (base = 0 degrees)      h      O  *
// DAYLP   R4  Photoperiodic daylength (base = -4 degrees)    h      O  *
// DSINB   R4  Daily total of sine of solar height            s      O  *
// DSINBE  R4  Daily integral of sine of solar height         s      O  *
//             corrected for lower transmission at low                  *
//             elevation                                                *
// SINLD   R4  Intermediate variable for subroutine SSKYC     -      O  *
// COSLD   R4  Intermediate variable for subroutine SSKYC     -      O  *
//                                                                      *
// Fatal error checks: LAT > 90, LAT < -90                              *
// Warnings          : LAT above polar circle, LAT within polar circle  *
// Subprograms called: FATALERR                                         *
// File usage        : none                                             *
//----------------------------------------------------------------------*

#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

void SASTRO(double IDOY, double LAT, double &SOLCON, double &ANGOT, double &DAYL, double &DAYLP, double &DSINB, double &DSINBE,
            double &SINLD, double &COSLD){

    double AOB, DOY, DEC, ZZCOS, ZZSIN, ZZA;
    double PI = 3.1415927;
    double DEGTRAD=0.017453292;


    //     Error check and conversion of day number
    if ( fabs(LAT) > 90.)  //problem messgae

    DOY = double(IDOY);

//     Declination of the sun as a function of daynumber,
//     calculation of daylength from intermediate variables
//     SINLD, COSLD and AOB

    DEC = - asin( sin(23.45*DEGTRAD)*cos(2.*PI*(DOY+10.)/365.));
    SINLD = sin(DEGTRAD*LAT)*sin(DEC);
    COSLD = cos(DEGTRAD*LAT)*cos(DEC);
    AOB = SINLD/COSLD;

    if (AOB < -1.)  {
        //error message
        DAYL  = 0.;
        ZZCOS = 0.;
        ZZSIN = 1.;
    }
    else if(AOB > 1.){
        //error message
        DAYL  = 24.;
        ZZCOS =  0.;
        ZZSIN = -1.;
    }
    else{
        DAYL  = 12.*(1.+2.* asin(AOB)/PI);
        DAYLP = 12.0*(1.+2.* asin((-sin(-4.*DEGTRAD)+SINLD)/COSLD)/PI);
        ZZA   = PI*(12.+DAYL)/24.;
        ZZCOS = cos(ZZA);
        ZZSIN = sin(ZZA);
    }

    //     Daily integral of sine of solar height (DSINB) with a
    //     correction for lower atmospheric transmission at lower solar
    //     elevations (DSINBE)
    DSINB  = 2.*3600.*(DAYL*0.5*SINLD-12.*COSLD*ZZCOS/PI);
    DSINBE = 2.*3600.*(DAYL*(0.5*SINLD+0.2*SINLD*SINLD+0.1*COSLD*COSLD) - (12.*COSLD*ZZCOS+9.6*SINLD*COSLD*ZZCOS + 2.4*COSLD*COSLD*ZZCOS*ZZSIN)/PI);

    //     Solar constant and daily extraterrestrial radiation
    SOLCON = 1370.*(1.+0.033* cos(2.*PI*DOY/365.));
    ANGOT  = SOLCON*DSINB;
}