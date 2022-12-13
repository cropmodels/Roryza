//----------------------------------------------------------------------//
//  SUBROUTINE ET2                                                      //
//  Used in ORYZA model version 4.0                                     //
//  Date  : December 2001; modified May 30, 2006                        //
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

#include <math.h>
#include <algorithm>
#include "model.h"
#include <string>

// using namespace std;

void ET2_initialize(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
        int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC){
    //     Local variables
    int       ISURF;
    double          ALB, DT, ETAE, ETRD;
    double          RF , RFS;

    std::string ESTAB;
    double ETDCUM1, EVSCCUM1, TRCCUM1;
    double ETDCUM2, EVSCCUM2, TRCCUM2;
    double ETDCUM3, EVSCCUM3, TRCCUM3;


    ETD = 0.;
    EVSC = 0.;
    TRC = 0.;

    ETDCUM1  = 0.;
    EVSCCUM1 = 0.;
    TRCCUM1  = 0.;

    ETDCUM2  = 0.;
    EVSCCUM2 = 0.;
    TRCCUM2  = 0.;

    ETDCUM3  = 0.;
    EVSCCUM3 = 0.;
    TRCCUM3  = 0.;




}


void ET2_rate(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
              int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC){
    //     Local variables
    int       ISURF;
    double          ALB, DT, ETAE, ETRD;
    double          RF , RFS;

    std::string ESTAB;
    double ETDCUM1, EVSCCUM1, TRCCUM1;
    double ETDCUM2, EVSCCUM2, TRCCUM2;
    double ETDCUM3, EVSCCUM3, TRCCUM3;


}


void ET2_state(double ANGA, double ANGB, double RDD, double TMDA, double VP, double WN, double LAT, int IDOY, std::string ETMOD, int CROPSTA,
               int NL, double FAOF, double WL0, double WCLOT, double WCST, double LAI, double &EVSC, double &ETD, double &TRC){
    //     Local variables
    int       ISURF;
    double          ALB, DT, ETAE, ETRD;
    double          RF , RFS;

    std::string ESTAB;
    double ETDCUM1, EVSCCUM1, TRCCUM1;
    double ETDCUM2, EVSCCUM2, TRCCUM2;
    double ETDCUM3, EVSCCUM3, TRCCUM3;



}















