//----------------------------------------------------------------------*
//  SUBROUTINE SUBGRN                                                   *
//  Adapted: Bouman, July 1999                                          *
//  Purpose: This subroutine calculates spikelet formation rate and     *
//           spikelet fertility as affected by low and high temperature *
//           and the grain growth rate.  Spikelet sterility component   *
//           is according to Horie et al., 1992.                        *
//                                                                      *
// FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
// name   type meaning (unit)                                     class *
// ----   ---- ---------------                                    ----- *
// GCR     R4  Gross growth rate of the crop (kg DM/ha/d)            I  *
// CROPSTA R4  Crop stage (-)                                        I  *
// LRSTRS  R4  Leaf rolling stress factor (-)                        I  *
// DVS     R4  Development stage of the crop (-)                     I  *
// SF1     R4  Spikelet fertility factor due to low temperatures (-) I  *
// SF2     R4  Spikelet fertility factor due to high temperatures (-)I  *
// SPGF    R4  Spikelet growth factor (no kg-1)                      I  *
// TAV     R4  Average daily temperature (oC)                        I  *
// TMAX    R4  Daily maximum temperature (oC)                        I  *
// NSP     R4  Number of spikelets (no)                              I  *
// TIME    R4  Time of simulation (d)                                T  *
// GNSP    R4  Rate of increase in spikelet number (no ha-1 d-1)     O  *
// GNGR    R4  Rate of increase in grain number (no ha-1 d-1)        O  *
// SPFERT  R4  Spikelet fertility (-)                                O  *
// GRAINS  L*  Fortran logical function whether grains are formed    O  *
//                                                                      *
//  FILE usage : none                                                   *
//----------------------------------------------------------------------*

#include <algorithm>
#include <math.h>
#include "model.h"

// using namespace std;

void SUBGRN( double GCR, double CROPSTA, double LRSTRS, double DVS, double SF1, double SF2, double SPGF,
            double TAV, double TMAX, double NSP, double TIME, double &GNSP, double &GNGR, double &SPFERT, bool &GRAINS){

    //-----Local parameters
    double    DVSPI, DVSF, CTT, COLDTT, TFERT, NTFERT, TINCR;

    //     Initialization
    if (CROPSTA  <=  1)  {
        GRAINS = false;
        COLDTT = 0.;
        TFERT  = 0.;
        NTFERT = 0.;
        SF1    = 1.;
        SF2    = 1.;
        SPFERT = 1.;
    }

    //-----Temperature increase due to leaf rolling (BAS): 1.6 degree
    //     per unit leaf rolling (Turner et al., 1986; p 269)
    TINCR = 5.*(1.-LRSTRS)*1.6;

    //-----Spikelet formation between PI and Flowering
    DVSPI = 0.65;
    DVSF  = 1.;
    if ((DVS >= DVSPI)  && (DVS <= DVSF))  {
        GNSP = GCR*SPGF;
    } else {
        GNSP = 0.;
    }

    //-----Grain formation from spikelets (GNGR)
    //-----Calculate GNGR reduction factors
    if ((DVS >= 0.75)  && (DVS <= 1.2))  {
        CTT    = std::max(0.,22.-(TAV-TINCR));
        COLDTT = COLDTT+CTT;
    }
    if ((DVS >= 0.96)  && (DVS <= 1.2))  {
        TFERT  = TFERT +(TMAX+TINCR);
        NTFERT = NTFERT+1.;
    }

    //-----Apply GNGR reduction factors when DVS is 1.2
    if ((DVS >= 1.2)  && !GRAINS )  {
        GRAINS = true;
        SF1    = 1.-(4.6+0.054* pow(COLDTT, 1.56))/100.;
        SF1    = std::min(1.,std::max(0.,SF1));
        TFERT  = TFERT/(NTFERT);
        SF2    = 1./(1.+exp(0.853*(TFERT-36.6)));
        SF2    = std::min(1.,std::max(0.,SF2));
        SPFERT = std::min(SF1,SF2);
        GNGR   = NSP*SPFERT;
    } else {
        GNGR   = 0.;
    }

}








