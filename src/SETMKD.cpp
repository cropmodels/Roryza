//----------------------------------------------------------------------*
// SUBROUTINE SETMKD (Subroutine Evap. Trans. MaKkink Daily)            *
// Authors: Daniel van Kraalingen                                       *
// Date   : 7-March-1997                                                *
// Version: 1.0                                                         *
// Purpose: This subroutine calculates reference evapotranspiration     *
//          according to Makkink (1957). To obtain crop evapo-          *
//          transpiration, multiplication with a Makkink crop factor    *
//          should be done. The use of this formula is basically limited*
//          to areas with large amounts of radiation                    *
// Refs.  : Kraalingen, D.W.G. van, W. Stol, 1997. Evapotranspiration   *
//          modules for crop growth simulation. Quantitative Approaches *
//          in Systems Analysis No. 11. DLO Research Institute for      *
//          Agrobiology and Soil Fertility (AB-DLO), The C.T. de Wit    *
//          graduate school for Production Ecology (PE). Wageningen.    *
//          The Netherlands.                                            *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (units)                                    class *
// ----   ---- ---------------                                    ----- *
// RDD     R4  Daily short-wave radiation (J.m-2.d)                  I  *
// TMDA    R4  24 hour average temperature (degrees C)               I  *
// ETD     R4  Potential evapotranspiration (mm.d-1)                 O  *
//                                                                      *
// Fatal error checks : none                                            *
// Warnings           : RDD < 0.5E6                                     *
// Subprograms called : SVPS1                                           *
// Required libraries : none                                            *
// File usage         : none                                            *
//----------------------------------------------------------------------*

#include <math.h>
#include "model.h"
// using namespace std;

double SETMKD(double RDD, double TMDA){
    double LHVAP, PSCH, VPS, VPSL, MAKFAC;
    LHVAP = 2454.e3;
    PSCH = 0.067;
    MAKFAC = 0.63;
    //     Checks
    if(RDD < 0.5e6) {
        //problem message
    }

    std::vector<double> result = SVPS1(TMDA);
    VPS = result[0];
    VPSL = result[1];
    //     Calculate Makkink evaporation, MAKFAC factor is calibrated for the
    //     Netherlands
    double ETD = MAKFAC*(RDD*(VPSL/(VPSL+PSCH)))/LHVAP;
    return ETD;
}