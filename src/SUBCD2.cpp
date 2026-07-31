//----------------------------------------------------------------------*
//  SUBROUTINE SUBCD2                                                    *
//  Consecutive cold days counter                                        *
//----------------------------------------------------------------------*

#include "model.h"

void SUBCD2(double COLDMIN, int CROPSTA, double TAV, double &NCOLD) {
	if (CROPSTA == 3) NCOLD = 0.;
	if (TAV < COLDMIN) {
		NCOLD = NCOLD + 1.;
	} else {
		NCOLD = 0.;
	}
}
