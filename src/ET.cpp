//----------------------------------------------------------------------//
//  SUBROUTINE ET                                                       //
//  Used in ORYZA model version 4.0                                     //
//  Date  : December 2001                                              //
//  Author: B.A.M. Bouman                                               //
//                                                                      //
//  Purpose: Calculates potential evaporation of soil/water layer and   //
//           potential transpiration of a crop. Calculations done with: //
//           Penman, Priestley-Taylor or Makkink subroutines.           //
//           All calculations pertain to the main field, and not to     //
//           the seedbed.                                              //
//                                                                      //
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      //
// name   type meaning (unit)                                     class //
// ----   ---- ---------------                                    ----- //
// ITASK   I4  Task that subroutine should perform (-)               I  //
// ANGA    R4  Angstrom parameter A                                  I  //
// ANGB    R4  Angstrom parameter B                                  I  //
// RDD     R4  Daily shortwave radiation (J.m-2.d-1)                 I  //
// TMDA    R4  Daily average temperature (degrees C)                 I  //
// VP      R4  Early morning vapour pressure (kPa)                   I  //
// WN      R4  Average wind speed (m.s-1)                            I  //
// LAT     R4  Latitude of site (dec.degr.)                          I  //
// IDOY    I4  Day number within year of simulation (d)              I  //
// ETMOD   C*  Name of subroutine to calculate E and T (-)           I  //
// CROPSTA I4  Crop stage (-)                                        I  //
// NL      I4  Number of soil layers (-)                             I  //
// FAOF    R4  Correction factor for E and T (FAO factor) (-)        I  //
// WL0     R4  Depth f ponded water layer (mm)                       I  //
// WCLQT   R4  Array of actual soil water contents/layer (m3 m-3)    I  //
// WCST    R4  Array of water content saturation / layer (m3 m-3)    I  //
// LAI     R4  Leaf Area Index (-)                                   I  //
// EVSC    R4  Potential soil evaporation (mm d-1)                   O  //
// ETD     R4  Reference evapotranspiration (mm d-1)                 O  //
// TRC     R4  Potential transpiration of crop at given LAI (mm d-1) O  //
//                                                                      //
// SUBROUTINES called: SETPMD, SETMKD, SETPTD                           //
//                                                                      //
// Files included: -                                                    //
//----------------------------------------------------------------------//

#include <vector>
#include <algorithm>
#include <math.h>
#include <string>
#include "model.h"


// using namespace std;

void ET(int ITASK, double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, double IDOY, std::string ETMOD,
                    double CROPSTA, double NL, double FAOF, double WL0, std::vector<double> WCLQT, std::vector<double> WCST, double LAI,
                    double &EVSC, double &ETD, double &TRC){

    
    //     Local variables
    int       ISURF;
    double          ALB, DT, ETAE, ETRD;
    double          RF , RFS;
    //SAVE;

    //output
    //double EVSC, ETD, TRC;

    if(ITASK == 2){
        //---- Set value for reflection coefficient of soil or water background
        //     If there is standing water:
        if(WL0 > 5.){
            ALB = 0.05;
            RFS = ALB;
            //     If there is moist or dry soil
        } else{
            ALB = 0.25;
            RFS = ALB*(1. - 0.5*WCLQT[0]/WCST[0] );
        }

        //---- The soil or water background is shielded by the crop
        RF = RFS * exp(-0.5*LAI) + 0.25*(1. - exp(-0.5*LAI));

        //-----Penman evapotranspiration
        if(ETMOD == "PENMAN"){
            //-----Set ISURF value (soil or water background) for wind function in main field
            //     Before transplanting: ISURF equals 1=open water, or 2 = bare soil)
            //     After transplanting : ISURF equals 3
            if(CROPSTA < 3.){
                if(WL0 > 5.){
                    ISURF = 1;
                } else{
                    ISURF = 2;
                }
            } else{
                ISURF = 3;
            }
            std::vector<double> setpmd = SETPMD(IDOY,LAT,ISURF,RF,ANGA,ANGB,0.,RDD,TMDA,WN,VP);
            ETD = setpmd[0];
            ETRD = setpmd[1];
            ETAE = setpmd[2];
            DT = setpmd[3];
        }
        else if(ETMOD == "MAKKINK"){
            ETD = SETMKD(RDD, TMDA);
            //        Estimate radiation-driven and wind- and humidity-driven part
            ETRD = 0.75*ETD;
            ETAE = ETD-ETRD;
        }
        else if(ETMOD == "PRIESTLEY TAYLOR"){
            ETD = SETPTD(IDOY,LAT,RF,RDD,TMDA);
            //        Estimate radiation-driven and wind- and humidity-driven part
            ETRD = 0.75*ETD;
            ETAE = ETD-ETRD;
        }
        //-----Multiplied by a factor according to FAO (1998)
        ETD  = ETD  * FAOF;
        ETRD = ETRD * FAOF;
        ETAE = ETAE * FAOF;
        //---- Calculate potential soil evaporation taking into account the standing crop
        EVSC = exp(-0.5*LAI)*(ETRD+ETAE);
        EVSC = std::max(EVSC, 0.);

        //---- Calculate potential transpiration of rice in main field
        if (CROPSTA  >=  4)  {
            TRC = ETRD*(1.-exp(-0.5*LAI))+ETAE*std::min(2.0,LAI);
        //     There is no transpiration from main field before transplanting
        } else {
            TRC = 0.;
        }
        //     End ITASK=2:
    }

}










