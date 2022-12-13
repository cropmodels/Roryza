//----------------------------------------------------------------------*
//  SUBROUTINE SUBCD                                                    *
//  Purpose: This subroutine calculates number of days below a certain  *
//           average temperature (TAV), which is used to terminate the  *
//           simulation after a maximum number of cold days the crop    *
//           can survive.                                               *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// CROPSTA I4  Crop stage (-)                                        I  *
// TAV     R4  Average daily temperature (oC)                        I  *
// TIME    R4  Time of simulation (d)                                T  *
// NCOLD   R4  Number of cold days (-)                               O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

double SUBCD(int CROPSTA, double TAV, double TIME){
    //output
    double NCOLD = 0.;
    if (CROPSTA  ==  3) {
        NCOLD = 0.;
    }

    if( TAV < 12. ){
        NCOLD = NCOLD+1.;
    } else{
        NCOLD = 0.;
    }

    if (NCOLD > 3.){
        //error messag
    }
    return NCOLD;
}








