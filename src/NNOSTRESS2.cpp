//----------------------------------------------------------------------*
// SUBROUTINE NNOSTRESS                                                 *
// Authors: Bas Bouman                                                  *
// Date   : Jan-2001, Version: 1.0                                      *
//          Version august, 2003                                        *
// Purpose: This subroutine sets nitrogen stress factors on crop growth *
//          to unity, and gets leaf N content as function of DVS or     *
//          from observed values in the experiment data file.           *
//          A programming 'trick' using the function INTGR2 is used to  *
//          get DVS-interpolated values for leaf N content before and   *
//          after the first and last day of observation, respectively.  *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// DELT    R4  Time step of integration (d)                          I  *
// IUNITD  I4  Unit that can be used for input files (-)             I  *
// IUNITL  I4  Unit number for log file messages (-)                 I  *
// ITASK   I4  Task that subroutine should perform (-)               I  *
// FILEI1  C*  Name of file with crop input data (-)                 I  *
// FILEIT  C*  Name of file with experimental data (-)               I  *
// CROPSTA I4  Crop stage (-)                                        I  *
// DVS     R4  Development stage (-)                                 I  *
// NFLV    R4  Nitrogen fraction in the leaves (g N m-2 leaf)        O  *
// NSLLV   R4  N stress factor on leaf death (-)                     O  *
// RNSTRS  R4  N stress factor on relative leaf growth (-)           O  *
//                                                                      *
//----------------------------------------------------------------------*
#include <math.h>
#include <string>
#include <algorithm>
#include "model.h"
#include "oryzaUtil.h"
#include <vector>
// using namespace std;

/*
double AFGEN(std::vector<double> xy, double x) {
    int n = xy.size();
    double y = -1;
    if (x < xy[0]) {
        y = xy[1];
    } else if (x > xy[n-2]) {
        y = xy[n-1];
    } else {
        for(int i=2; i<n; i=i+2) {
            if (xy[i] > x) {
                double slope = (xy[i+1] - xy[i-1]) / (xy[i] - xy[i-2]);
                y = xy[i-1] + (x - xy[i-2]) * slope;
                break;
            }
        }
    }
    return(y);
}

*/
void NNOSTRESS2_initialization( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS ){
    NFLV     = NFLVI;
    NSLLV    = 1.;
    RNSTRS   = 1.;

}

void NNOSTRESS2_rate( double NFLVI, std::vector<double> NMAXLT, std::vector<double> NFLVTB, double DELT, int CROPSTA, double DVS, double WLVG, double LAI, double SLA, double &NFLV, double &NSLLV, double &RNSTRS ){
    if(CROPSTA < 4){
        NFLV   = NFLVI;
    } else{
        //          Read NFLV as function of development state
        //            NFLV1 = LINT2('NFLVTB',NFLVTB,ILNFLV,DVS)
        double NMAXL = AFGEN (NMAXLT, DVS);
        double NFLV1 = NMAXL/(10.*SLA);
        
        //          This is a programming 'trick' to enable forcing of observed
        //          NFLV as function of day-of-observation. Below the first observation
        //          day, NFLV1 is used;l between first and last observation day, interpolated
        //          observed values are used; after last observation day, NFLV1 is used again
        //          Forcing is determined by the variable NFLV_FRC in the experiment data file.
        NFLV = NFLV1 * DELT;
    }
}

