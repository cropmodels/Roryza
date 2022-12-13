//----------------------------------------------------------------------*
//  SUBROUTINE SUBCBC                                                   *
//  Purpose: This subroutine checks the Crop Carbon Balance             *
//           and stops the simulation if the difference between         *
//           CKCIN and CKCFL exceeds 0.1 %                              *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// CKCIN   R4  Accumulated C in the crop (kg C ha-1)                 I  *
// CKCFL   R4  Sum of integrated C fluxes (kg C ha-1)                I  *
// TIME    R4  Time of simulation (d)                                T  *
// CBCHK   R4  Carbon balance check, relative value to the sums of      *
//             CKIN and CKCFL (-)                                    O  *
// TERMNL  R4  Flag to indicate if simulation is to stop (-)         O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <vector>
#include <math.h>
#include <algorithm>
#include "model.h"

// using namespace std;

void SUBCBC( double CKCIN, double CKCFL, double TIME, double &CBCHK, bool &TERMNL ){

    //output
    //double CBCHK;
    //double TERNML = 0.;

    CBCHK = 2.0*(CKCIN-CKCFL)/(CKCIN+CKCFL+1.E-10);

    if (fabs(CBCHK) > 0.001)  {
        //error message
        TERMNL = true;
    }
    //return {CBCHK, TERNML};
}