//----------------------------------------------------------------------*
//  SUBROUTINE SUBCD2                                                    *
//  Purpose: This subroutine calculates number of days below a certain  *
//           average temperature (TAV), which is used to terminate the  *
//           simulation after a maximum number of cold days the crop    *
//           can survive.                                               *
// Adapted from version 1, March 2006, Bouman                           *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// COLDMIN I4 lower T threshold for growth                           I  *
// CROPSTA I4  Crop stage (-)                                        I  *
// TAV     R4  Average daily temperature (oC)                        I  *
// TIME    R4  Time of simulation (d)                                T  *
// NCOLD   R4  Number of cold days (-)                               O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <math.h>
#include "model.h"
// using namespace std;

double SUBCD2(int COLDMIN, int CROPSTA, double TAV){

    //output
    double NCOLD;

    if (CROPSTA  ==  3) NCOLD = 0.;
    if (TAV < COLDMIN)  {
// RJH  NCOLD = 1.
        NCOLD = 100.;
// end RJH
    } else {
        NCOLD = 0.;
    }
    return NCOLD;
}

