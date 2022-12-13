#include <math.h>
#include <vector>
#include <algorithm>
#include "model.h"
//----------------------------------------------------------------------//
// SUBROUTINE NCROP                                                     //
// Authors:                                                             //
//          Version august, 2003                                        //
// Date   : December 2001, Version: 1                                   //
// Purpose: This subroutine calculates the nitrogen dynamics in a rice  //
//          crop and the effects on growth and development              //
//                                                                      //
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      //
// name   type meaning                                     units  class //
// ----   ---- -------                                     -----  ----- //
// ITASK   I4  Task that subroutine should perform (-)               I  //
// IUNITD  I4  Unit that can be used for input files (-)             I  //
// IUNITL  I4  Unit number for log file messages (-)                 I  //
// FILEI1  C*  Name of file with crop data (-)                       I  //
// DELT    R4  Time step of integration (d)                          I  //
// TIME    R4  Time of simulation (d)                                I  //
// OUTPUT  R4  Flag to indicate if output should be done (-)         I  //
// TERMNL  L4  Flag to indicate if simulation is to stop (-)         I  //
// DVS     R4  Development stage of the crop (-)                     I  //
// LLV     R4  Loss rate of leaves caused by senescence (kg ha-1 d-1)I  //
// DLDR    R4  Loss rate of leaves caused by drought (kg ha-1 d-1)   I  //
// WLVG    R4  Dry weight of green leaves (kg ha-1)                  I  //
// WST     R4  Dry weight of stems (kg ha-1)                         I  //
// WSO     R4  Dry weight of storage organs (kg ha-1)                I  //
// GSO     R4  Growth rate of storage organs (kg ha-1 d-1)           I  //
// GST     R4  Growth rate of stems (kg ha-1 d-1)                    I  //
// GLV     R4  Growth rate of leaves (kg ha-1 d-1)                   I  //
// PLTR    R4  Intermediate variable for planting density (-)        I  //
// LAI     R4  Leaf area index (ha ha-1)                             I  //
// CROPSTA I4  Crop stage (-)                                        I  //
// TNSOIL  R4  Soil-N available for crop uptake (kg N ha-1 d-1)      I  //
// NACR    R4  Actual N uptake rate by the crop (kg ha-1 d-1)        O  //
// NFLV    R4  Nitrogen fraction in the leaves (g N m-2 leaf)        O  //
// NSLLV   R4  Stress factor for leaf death caused by N stress  (-)  O  //
// RNSTRS  R4  Decrease factor for RGRL caused by N stress (-)       O  //
//                                                                      //
// Subroutine called: SUBNBC                                            //
//                                                                      //
//----------------------------------------------------------------------*


void NCROP2_initialization( double DELT, double TIME, bool TERMINAL, double DVS, double LLV, double DLDR, double WLVG, double WST,
                            double WSO, double GSO, double GST, double GLV, double PLTR, double LAI, double SLA, int CROPSTA, double TNSOIL,
                            double &NACR, double &NFLV, double &NSLLV, double &RNSTRS){

    //       Initialize variables
//    NUPP  = 0.;
//    ANLV   = 0.;
//    ANSO   = 0.;
//    ANST   = 0.;
//    ANLD   = 0.;
//    ANCR   = 0.;
//    ANLVA  = 0.;
//    ANSTA  = 0.;
//    ANCRF  = 0.;
//    NALVS  = 0.;
//    NASTS  = 0.;
//    NASOS  = 0.;
//    NACRS  = 0.;
//    NTRTS  = 0.;
//    NALV   = 0.;
//    NAST   = 0.;
//    NASO   = 0.;
//    NACR   = 0.;
//    NLV    = 0.;
//    NST    = 0.;
//    NSO    = 0.;
//    NLDLV  = 0.;
//    NLVAN  = 0.;
//    NSTAN  = 0.;
//    NTRT   = 0.;
//    ANCRPT = 0.;
//    NBCHK  = 0.;
//    NCHCK  = 0.;
//
//    FNLV   = FNLVI;
//    FNST   = 0.5*FNLVI;
//    FNSO   = 0.;
//    NFLV   = NFLVI;
//
//    NSLLV  = 1.;
//    RNSTRS = 1.;

}












//----------------------------------------------------------------------*
//  SUBROUTINE SUBNBC                                                   *
//  Purpose: This subroutine checks the Crop Nitrogen Balance           *
//           and stops the simulation if the difference between         *
//           CHKIN and CHKFL exceeds 0.1 %                              *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// CHKIN   R4  Accumulated N in the crop (kg N ha-1)                 I  *
// CHKFL   R4  Sum of N supplied by soil and roots (kg N ha-1)       I  *
// TIME    R4  Time of simulation (d)                                T  *
// NBCHK   R4  Nitrogen balance check, relative value to the sums of    *
//             CHKIN and CHKFL (-)                                   O  *
// TERMNL  R4  Flag to indicate if simulation is to stop (-)         O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*
void SUBNBC( double CHKIN, double CKCFL, double TIME, double &NBCHK, bool &TERMINAL ){
    //IMPLICIT NONE;
//-----Formal parameters
    NBCHK = 2.0*(CHKIN-CKCFL)/(CHKIN+CKCFL+1.E-10);

    if( fabs(NBCHK) > 0.001 ){
        //error message
        TERMINAL = true;
    }
}














