#include <math.h>
#include <algorithm>
#include "oryzaUtil.h"
#include "model.h"
double cCO2;
double cKNF;
double cNFLV;
double cREDFT;

// using namespace std;

void GPPARSET( double xCO2, double xKNF, double xNFLV, double xREDFT ){
    cCO2   = xCO2;
    cKNF   = xKNF;
    cNFLV  = xNFLV;
    cREDFT = xREDFT;
}