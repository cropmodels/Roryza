//----------------------------------------------------------------------//
//  SUBROUTINE WNOSTRESS                                                //
//  Used in ORYZA2000 model version 1.0                                 //
//  Date  : December 2001                                               //
//  Author: B.A.M. Bouman                                               //
//                                                                      //
//  Purpose: Invoked when production environment is POTENTIAL.          //
//           Sets actual transpiration of a crop at zero, and sets      //
//           effects of water stress on growth and development of rice  //
//           at unity.                                                  //
//                                                                      //
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      //
// name   type meaning (unit)                                     class //
// ----   ---- ---------------                                    ----- //
// NL      I4  Number of soil layers (-)                             I  //
// TRW     R4  Actual transpiration rate (mm d-1)                    O  //
// TRWL    R4  Array of actual transpiration rate/layer (mm d-1)     O  //
// LRSTRS  R4  Stress factor for rolling of leaves (-)               O  //
// LDSTRS  R4  Stress factor for accelerating leaf death (-)         O  //
// LESTRS  R4  Stress factor for reducing expansion of leaves (-)    O  //
// PCEW    R4  Stress factor for CO2 assimilation (-)                O  //
//                                                                      //
//----------------------------------------------------------------------//

#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

void WNOSTRESS( int NL, double &TRW, std::vector<double> &TRWL, double &LRSTRS, double &LDSTRS, double &LESTRS, double &PCEW,double &CPEW ){

    TRWL = std::vector<double>(NL, 0);
    TRW    = 0.;
    LRSTRS = 1.;
    LDSTRS = 1.;
    LESTRS = 1.;
    CPEW   = 1.;
    PCEW   = 1.;
}






