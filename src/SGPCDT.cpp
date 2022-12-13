//----------------------------------------------------------------------*
// SUBROUTINE SGPCDT                                                    *
// Authors: Daniel van Kraalingen                                       *
// Date   : 28-Apr-1995, Version: 1.2                                   *
// Purpose: This subroutine calculates daily total gross                *
//          assimilation (GPCDT) by performing a Gaussian integration   *
//          over time. At three different times of the day,             *
//          radiation is computed and used to determine assimilation    *
//          whereafter integration takes place. Accuracy of computation *
//          can be controlled the switch IACC (IACC=1 calls three point *
//          GAUSS procedure, IACC=2 calls combination routine with      *
//          error control                                               *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning                                     units  class *
// ----   ---- -------                                     -----  ----- *
// IACC    I4  Switch used to determine accuracy of            -     I  *
//             assimilation calculations                                *
// IDOY    I4  Day of year                                    d      I  *
// LAT     R4  Latitude of the site                        degrees   I  *
// RDD     R4  Daily total of shortwave radiation          J/m2/d    I  *
// FRPAR   R4  Fraction of total shortwave irradiation        -      I  *
//             that is photosynthetically active                        *
// CSLV    R4  Scattering coefficient of leaves for PAR        -     I  *
// AMAX    R4  Assimilation rate at light saturation       kg CO2/   I  *
//                                                        ha leaf/h     *
// EFF     R4  Initial light use efficiency               kg CO2/J/  I  *
//                                                       (ha/h/m2/s)    *
// ECPDF   R4  Extinction coefficient for          ha ground/ha leaf I  *
//             diffuse light                                            *
// GAI     R4  Green area index                    ha leaf/ha ground I  *
// DAYL    R4  Astronomic daylength (base = 0 degrees)        h      O  *
// DAYLP   R4  Photoperiodic daylength (base = -4 degrees)    h      O  *
// GPCDT   R4  Daily total gross assimilation            kg CO2/ha/d O  *
// RAPCDT  R4  Daily rate of absorbed PAR                  J/m/d     O  *
//                                                                      *
// Fatal error checks: none                                             *
// Warnings          : none                                             *
// Subprograms called: SASTRO, SGPC1, SGPC2, SSKYC                      *
// File usage        : none                                             *
//----------------------------------------------------------------------*


#include <vector>
#include <algorithm>
#include <math.h>
#include "model.h"

// using namespace std;

void SGPCDT(int IACC, int IDOY, double LAT, double RDD, double FRPAR, double CSLV, double AMAX, double EFF, double ECPDF,
            double GAI, double &DAYL, double &DAYLP, double &GPCDT, double &RAPCDT ){

    //     Miscellaneous
    double SOLCON,DSINB,DSINBE,SINLD,COSLD,HOUR,SINB,RDPDR,RDPDF;
    double GPC = 0.0, ANGOT,RAPC = 0.0;

    //     Gauss parameters
    double GSX[3] = {0.112702, 0.500000, 0.887298};
    double GSW[3] = {0.277778, 0.444444, 0.277778};

    //     Compute daylength and related data
    SASTRO(IDOY,LAT,SOLCON,ANGOT,DAYL,DAYLP,DSINB,DSINBE,SINLD,COSLD);

    //     Assimilation set to zero and three different times of the day
    GPCDT  = 0.;
    RAPCDT = 0.;

    for(int i = 0; i < 3; i++){
        //        At the specified HOUR, external radiation
        //        conditions are computed
        HOUR = 12.0+DAYL*0.5*GSX[i];
        SSKYC(HOUR, SOLCON, FRPAR, DSINBE, SINLD, COSLD, RDD, SINB, RDPDR, RDPDF);

        //        Assimilation rate of canopy is computed
        if (IACC == 1){
            //           Choose three point Gauss method
            std::vector<double> result = SGPC1(CSLV, AMAX, EFF,ECPDF, GAI, SINB, RDPDR, RDPDF);
            GPC = result[0];
            RAPC = result[1];
        }
        else if(IACC == 2){
            //           Choose routine with internal error control
            std::vector<double> result2 = SGPC2(CSLV, AMAX, EFF,ECPDF, GAI, SINB, RDPDR, RDPDF);
            GPC = result2[0];
            RAPC = result2[1];
        }
        else{
            //error message
        }
        GPCDT  = GPCDT +GPC *GSW[i];
        RAPCDT = RAPCDT+RAPC*GSW[i];
    }
    //     Integration of assimilation rate to a daily total (GPCDT)
    GPCDT  = GPCDT *DAYL;

    //     Integration of absorption rate to a daily total (RAPCDT)
    RAPCDT = RAPCDT*DAYL*3600.;

}











