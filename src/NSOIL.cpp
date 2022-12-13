//----------------------------------------------------------------------//
//  SUBROUTINE NSOIL                                                    //
//  Used in ORYZA2000 model version 1.0                                 //
//  Date   : December 2001                                              //
//  Author : B.A.M. Bouman                                              //
//  Purpose: This module calculates the supply of N from the soil.      //
//                                                                      //
//                                                                      //
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      //
// name   type meaning (unit)                                     class //
// ----   ---- ---------------                                    ----- //
// ITASK   I4  Task that subroutine should perform (-)               I  //
// IUNITD  I4  Unit that can be used for input files (-)             I  //
// IUNITL  I4  Unit number for log file messages (-)                 I  //
// FILEIT  C*  Name of file with experimental data (-)               I  //
// OUTPUT  R4  Flag to indicate if output should be done (-)         I  //
// DELT    R4  Time step of integration (d)                          I  //
// DAE     R4  Days after emergence (d)                              I  //
// DVS     R4  Development stage of the crop (-)                     I  //
// NACR    R4  Actual N uptake rate  by the crop (kg N ha-1 d-1)     I  //
// TNSOIL  R4  Total N in the soil for uptake by crop (kg ha-1 d-1)  O  //
//                                                                      //
// SUBROUTINES called: -                                                //
//                                                                      //
//----------------------------------------------------------------------//

#include <algorithm>
#include <math.h>
#include <vector>
#include <string>

// using namespace std;

double NSOIL( int ITASK, int IUNITD, int IUNITL, std::string FILEIT, double OUTPUT, double DELT, double DAE, double DVS, double NACR ){
    //     OUTPUT
    double TNSOIL;

//-----Local variables
    int INX, ILFERT, ILREC;
    INX=100;
    double FERT, RECOV, SOILSP, XFERT, NFERTP;
    double FERTIL[100];
    double RECNIT[100];



    return 0.;
}








