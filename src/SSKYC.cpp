//----------------------------------------------------------------------*
// SUBROUTINE SSKYC                                                     *
// Authors: Daniel van Kraalingen                                       *
// Date   : 2-Jun-1993, Version 1.0                                     *
// Purpose: This subroutine estimates solar inclination and fluxes of   *
//          diffuse and direct irradiation at a particular time of      *
//          the day.                                                    *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// HOUR    R4  Hour for which calculations should be done     h      I  *
// SOLCON  R4  Solar constant                               W/m2     I  *
// FRPAR   R4  Fraction of total shortwave irradiation        -      I  *
//             that is photosynthetically active (PAR)                  *
// DSINBE  R4  Daily integral of sine of solar height         s      I  *
//             corrected for lower transmission at low                  *
//             elevation                                                *
// SINLD   R4  Intermediate variable from SASTRO              -      I  *
// COSLD   R4  Intermediate variable from SASTRO              -      I  *
// RDD     R4  Daily shortwave radiation                   J/m2/d    I  *
// SINB    R4  Sine of solar inclination at HOUR              -      O  *
// RDPDR   R4  Instantaneous flux of direct photo-          W/m2     O  *
//             synthetically active radiation (PAR)                     *
// RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2     O  *
//             synthetically active irradiation (PAR)                   *
//                                                                      *
// Fatal error checks: RDD <= 0                                         *
// Warnings          : ATMTR > 0.9                                      *
// Subprograms called: FATALERR                                         *
// File usage        : none                                             *
//----------------------------------------------------------------------*


#include <math.h>
#include <algorithm>
#include "model.h"

void SSKYC( double HOUR, double SOLCON, double FRPAR, double DSINBE, double SINLD, double COSLD, double RDD, double &SINB, double &RDPDR, double &RDPDF ){
    //     Local parameters
    double ATMTR, FRDif, SOLHM, TMPR1;
    SOLHM=12.;

    if (RDD <= 0.) {}//error message

    ATMTR  = 0.;
    FRDif  = 0.;
    RDPDF = 0.;
    RDPDR = 0.;

    //     Sine of solar inclination, 0.2617993 is 15 degrees in radians
    SINB = SINLD+COSLD* cos((HOUR-SOLHM)*0.2617993);
    if (SINB > 0.){
        //        Sun is above the horizon
        TMPR1 = RDD*SINB*(1.+0.4*SINB)/DSINBE;
        ATMTR = TMPR1/(SOLCON*SINB);

        if (ATMTR > 0.9) {} //error message

        if (ATMTR <= 0.22) {
            FRDif = 1.;
        }
        else if(ATMTR > 0.22   &&  ATMTR <= 0.35){
            FRDif = 1.-6.4*pow(ATMTR-0.22, 2);
        }
        else{
            FRDif = 1.47-1.66*ATMTR;
        }
        //        Apply lower limit to fraction diffuse
        FRDif = std::max(FRDif, 0.15+0.85*(1.-exp(-0.1/SINB)));

        //        Diffuse and direct PAR
        RDPDF = TMPR1*FRPAR*FRDif;
        RDPDR = TMPR1*FRPAR*(1.-FRDif);
    }
}







