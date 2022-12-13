//----------------------------------------------------------------------*
// SUBROUTINE SGPC2                                                     *
// Authors: Daniel van Kraalingen                                       *
// Date   : 15-Dec-1993, Version: 1.1                                   *
// Purpose: Iterative procedure to calculate whole canopy assimilation  *
//          using internal error control. Routine can be used in place  *
//          of SGPC1 for greater accuracy. The method used is a combi-  *
//          nation of the TRAPZD and QSIMP routines from Numerical      *
//          Recipes (Press et al., 1990, page 111-113)                  *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// CSLV    R4  Scattering coefficient of leaves for          -       I  *
//             visible radiation (PAR)                                  *
// AMAX    R4  Assimilation rate at light saturation      kg CO2/    I  *
//                                                       ha leaf/h      *
// EFF     R4  Initial light use efficiency            kg CO2/J/ha/  I  *
//                                                       (h/m2/s)       *
// ECPDF   R4  Extinction coefficient for          ha ground/ha leaf I  *
//             diffuse light                                            *
// GAI     R4  Green area index                            ha/ha     I  *
// SINB    R4  Sine of solar height                          -       I  *
// RDPDR   R4  Instantaneous flux of direct photo-          W/m2     I  *
//             synthetically active radiation (PAR)                     *
// RDPDF   R4  Instantaneous flux of diffuse photo-         W/m2     I  *
//             synthetically active irradiation (PAR)                   *
// GPC     R4  Instantaneous assimilation rate of         kg CO2/    O  *
//             whole canopy                              ha soil/h      *
// RAPC    R4  Absorbed radiation                          W/m2      O  *
//                                                                      *
// Fatal error checks: N > NMAX (6)                                     *
// Warnings          : none                                             *
// Subprograms called: SGPL, FATALERR                                   *
// File usage        : none                                             *
//----------------------------------------------------------------------*

#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

std::vector<double> SGPC2( double CSLV, double AMAX, double EFF, double ECPDF, double GAI, double SINB, double RDPDR, double RDPDF ){


    //     Local variables

    double EPS, GPCTO, GPCT, RAPCT, RAPCTO, SUM1, SUM2, GPCO, RAPCO;
    double DEL, TNM, X, GPL1, GPL2, RAPL1, RAPL2;
    int NMAX, N, J, IT;
    EPS=0.10;
    NMAX=7;
    //output
    double GPC, RAPC;

    GPC   = 0.;
    GPCO  = 1.;
    GPCT  = 0.;

    RAPC  = 0.;
    RAPCO = 1.;
    RAPCT = 0.;

    N = 1;

    while( fabs(GPC-GPCO) > EPS*fabs(GPCO) || fabs(RAPC-RAPCO) > EPS* fabs(RAPCO) ){
        GPCO   = GPC;
        GPCTO  = GPCT;

        RAPCO  = RAPC;
        RAPCTO = RAPCT;
        if(N == 1){
            SGPL (CSLV, AMAX, EFF,ECPDF, GAI, 0., SINB, RDPDR, RDPDF, GPL1 , RAPL1);
            SGPL (CSLV, AMAX, EFF,ECPDF, GAI, GAI, SINB, RDPDR, RDPDF, GPL2 , RAPL2);
            //           First approximation of total gross photosynthesis of canopy
            GPCT  = 0.5*GAI*(GPL1 +GPL2);
            //           First approximation of absorption of PAR by canopy
            RAPCT = 0.5*GAI*(RAPL1+RAPL2);
        } else{
            IT = pow(2, N - 2);
            TNM  = IT;
            DEL  = GAI/TNM;
            X    = 0.5*DEL;
            SUM1 = 0.;
            SUM2 = 0.;

            for(int i = 0; i < IT; i++){
                SGPL (CSLV, AMAX, EFF,ECPDF, GAI, X, SINB, RDPDR, RDPDF, GPL2 , RAPL2);
                SUM1 = SUM1+GPL2;
                SUM2 = SUM2+RAPL2;
                X    = X+DEL;
            }

            //           Higher order approximation of total gross photosynthesis
            //           of canopy
            GPCT  = 0.5*(GPCT +GAI*SUM1/TNM);

            //           Higher order approximation of absorption of PAR by
            //           canopy
            RAPCT = 0.5*(RAPCT+GAI*SUM2/TNM);
        }
        GPC  = (4.*GPCT -GPCTO)/3.;
        RAPC = (4.*RAPCT-RAPCTO)/3.;

        N = N+1;
        if (N > NMAX)  {} //error message
    }

    return {GPC, RAPC};

}