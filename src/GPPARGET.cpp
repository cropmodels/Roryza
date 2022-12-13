
// user subroutine ////////////////

// xGAI     - Total leaf area index
// xGAID    - Leaf area index above point of calculation
// xAmaxIn  - Uncorrected amax
// xEffIn   - Uncorrected efficiency
// xAmaxOut - Corrected amax
// xEffOut  - Corrected efficiency

#include <math.h>
#include <algorithm>
#include "oryzaUtil.h"
#include "model.h"
extern double cCO2, cKNF, cNFLV, cREDFT;

// using namespace std;

void GPPARGET(double xGAI, double xGAID, double xAmaxIn, double xEffIn, double &xAmaxOut, double &xEffOut){

    //     local variables
    double AmaxCO2,SLNI;
    double Amax;
    double tmpr1;
    
    //     avoid compiler warnings on unused variables
    tmpr1 = xAmaxIn;

    AmaxCO2 = 49.57/34.26*(1.-exp (-0.208*(cCO2-60.)/49.57));
    AmaxCO2 = max (0.,AmaxCO2);

    if (xGAI > 0.01  && cKNF > 0.) {
        SLNI = cNFLV*xGAI*cKNF*exp(-cKNF*xGAID)/(1.-exp(-cKNF*xGAI));
    } else{
        SLNI = cNFLV;
    }

    //-----Calculate actual photosynthesis from SLN, CO2 and temperature
    //     calculation of AMAX according to Van Keulen  Seligman (1987):
    //     AMAX = 32.4 * (SLNI-0.2) * REDFT * CO2AMX

    if (SLNI >= 0.5){
        //        According to Shaobing Peng (IRRI, unpublished data):
        Amax = 9.5+(22.*SLNI)*cREDFT*AmaxCO2;
    } else{
        Amax = std::max(0.,68.33*(SLNI-0.2)*cREDFT*AmaxCO2);
    }

    xAmaxOut = Amax;
    xEffOut  = xEffIn;

}