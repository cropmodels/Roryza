//----------------------------------------------------------------------*
// SUBROUTINE SGPC1                                                     *
// Authors: Daniel van Kraalingen                                       *
// Date   : 15-Dec-1993, Version: 1.0                                   *
// Purpose: This subroutine performs a Gaussian integration over        *
//          depth of canopy by selecting three different GAI's and      *
//          computing assimilation at these GAI levels. The             *
//          integrated variable is GPC. Also absorbed PAR is            *
//          integrated in RAPC                                          *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// CSLV    R4  Scattering coefficient of leaves for           -      I  *
//             visible radiation (PAR)                                  *
// AMAX    R4  Assimilation rate at light saturation       kg CO2/   I  *
//                                                        ha leaf/h     *
// EFF     R4  Initial light use efficiency               kg CO2/J/  I  *
//                                                        ha/h m2 s     *
// ECPDF   R4  Extinction coefficient for          ha ground/ha leaf I  *
//             diffuse light                                            *
// GAI     R4  Green area index                             ha/ha    I  *
// SINB    R4  Sine of solar height                           -      I  *
// RDPDR   R4  Instantaneous flux of direct photo-          W/m2     I  *
//             synthetically active radiation (PAR)                     *
// RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2     I  *
//             synthetically active irradiation (PAR)                   *
// GPC     R4  Instantaneous assimilation rate of          kg CO2/   O  *
//             whole canopy                               ha soil/h     *
// RAPC    R4  Absorbed PAR                                 W/m2     O  *
//                                                                      *
// Fatal error checks: none                                             *
// Warnings          : none                                             *
// Subprograms called: SGPL                                             *
// File usage        : none                                             *
//----------------------------------------------------------------------*


#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

std::vector<double> SGPC1( double CSLV, double AMAX, double EFF, double ECPDF, double GAI, double SINB, double RDPDR, double RDPDF ){

    //output
    double GPC, RAPC;
    //     Gauss weights for three point Gauss
    double GSX[3] = {0.112702, 0.500000, 0.887298};
    double GSW[3] = {0.277778, 0.444444, 0.277778};

    //     Miscellaneous
    double GAID, GPL, RAPL;
    int I1;

    for(int i = 0; i < 3; i++){
        GAID = GAI * GSX[i];
        SGPL(CSLV , AMAX , EFF , ECPDF, GAI, GAID, SINB, RDPDR, RDPDF, GPL , RAPL);
        //        Integration of local assimilation rate to canopy
        //        assimilation (GPC) and absorption of PAR by canopy (RAPC)
        GPC  = GPC  + GPL  * GSW[i];
        RAPC = RAPC + RAPL * GSW[i];
    }
    GPC  = GPC  * GAI;
    RAPC = RAPC * GAI;
    return {GPC, RAPC};

}