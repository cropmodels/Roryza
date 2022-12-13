#include <vector>
#include <math.h>
#include <algorithm>
#include <string>
#include "model.h"
// using namespace std;


void oryza_model::weather_step() {
    atm.TMMN = wth.tmin[time];
    atm.TMMX = wth.tmax[time];
    //atm.TEMP  = (atm.TMIN + atm.TMAX) / 2.;
    //atm.DTEMP = (atm.TMAX + atm.TEMP) / 2.;

    //atm.AVRAD = wth.srad[time] * 1000;
    //atm.WIND = wth.wind[time];
    //atm.VAP = wth.vapr[time] * 10;
    //atm.RAIN = wth.prec[time] / 10 ; // cm !
    
    //oroblem
    //DOY = wth.simdate[time].dayofyear();

}


void oryza_model::model_output(){
    out.push_back( { double(step), crop.DVS, crop.LAI});
}

//void oryza_model::model_initialize(){
//    DELT = 1.;
//    TERMINAL = false;
//    fatal_error = false;
//
//    oryza_initialize();
//
//
//}


void oryza_model::model_run(){
    int nruns = 1;
    int maxstep = 300;
    //std::vector<std::vector<double> > out;

    step = 0;
    //time = 110.;
    for (int run=0; run < nruns; run++) {
        //crop.alive =true;
        TERMINAL = false;
        model_initialize();

        while ((!TERMINAL) && (step < maxstep)) {
            weather_step();

            model_rate();
            model_state();
            model_output();

            //time++;
            step++;
            //model_date++;
            if (fatal_error) {
                break;
            }

        }
    }
    out.resize(step);

}



void oryza_model::model_initialize(){

    DELT = 1.;
    TERMINAL = false;
    fatal_error = false;

    if(control.RUNMODE == "EXPLORATION"){
        crop.SWR = IYEAR;
        //problem
        //crop.EMD = NINT(time);
    }
//--------Initialize variables
    crop.CROPSTA  = 0;
    crop.WL0      = 0.;
    crop.RAINCU   = 0.;
    for(int i = 0; i < 10; i++){
        crop.WCLQT[i] = 0.3;
        crop.WCST[i] = 0.3;
    }

    //problem
    //--------Choose rice type to be modelled
//    if (RICETYPE == 'LOWLAND')  {
//        WRITE (IUNITO,'(A,T7,A)') '*','Lowland rice modelled';
//    } else { if (RICETYPE == 'AEROBIC')  {
//            WRITE (IUNITO,'(A,T7,A)') '*','Upland or aerobic rice modelled';
//        } else {
//            FATALERR ('MODELS','unknown name for RICETYPE');
//        }


    //skip
    atm.TMDA = (atm.TMMN + atm.TMMX) / 2.;
    
//    ET2_initialize(crop.ANGA, crop.ANGB, crop.RDD, atm.TMDA, atm.VP, atm.WN, atm.latitude, IDOY, control.ETMOD, crop.CROPSTA, NL, crop.FAOF, crop.WL0, crop.WCLQT, crop.WCST, crop.LAI, crop.EVSC, crop.ETD, crop.TRC);

    //problem ET2
    //ET2_initialize(crop.ANGA, crop.ANGB, crop.RDD, atm.TMDA, atm.VP, atm.WN, atm.latitude, IDOY, control.ETMOD, crop.CROPSTA, NL, crop.FAOF, crop.WL0, crop.WCLQT, crop.WCST, crop.LAI, crop.EVSC, crop.ETD, crop.TRC);

    if(control.PRODENV == "POTENTIAL"){
        WNOSTRESS(NL, crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
    }

    //call oryza
    oryza_initialize();

    if(control.NITROENV == "POTENTIAL"){
        NNOSTRESS2_initialization( crop.NFLVI, crop.NMAXLT, crop.NFLVTB, DELT, crop.CROPSTA, crop.DVS, crop.WLVG, crop.LAI, crop.SLA, crop.NFLV, crop.NSLLV, crop.RNSTRS );
    }

    crop.NFLV = crop.NFLVI;
    crop.NSLLV = 1.;
    crop.RNSTRS = 1.;


    if(control.PRODENV == "POTENTIAL"){
        crop.TKLT = 100.;
        crop.ZRTMS = 100.;
        crop.WL0 = 0.;
        NL = crop.NLXM;
        for(int i = 0; i < NL; i++){
            crop.WCLQT[i] = 0.3;
            crop.WCST[i] = 0.3;
        }
    }


    if (crop.CROPSTA  ==  3) crop.CROPSTA = 4;

    if (crop.CROPSTA  ==  2)  {
        if (crop.DAE  ==  double(crop.SBDUR)) crop.CROPSTA = 3;
// 	 GIS. direct seeded. NO TRANSPLANTING SHOCK
        if (control.ESTAB == "DIRECT-SEED") crop.CROPSTA = 4;
    }
    if (crop.CROPSTA  ==  1) {
        if (control.ESTAB == "TRANSPLANT") {
            crop.CROPSTA = 2;
        } else if (control.ESTAB == "DIRECT-SEED") {
//               CROPSTA = 4
// GIS to allow initialization of the water balance.
            crop.CROPSTA = 2;
        }
    }

    if (crop.CROPSTA  ==  0)  {
        //problem
        //crop.IDATE = DTFSECMP(crop.EMYR, crop.EMD, IYEAR, IDOY);
        if (crop.IDATE  ==  0)  {
            crop.CROPSTA = 1;
        } else if (crop.IDATE  ==  1)  {
                //problem error message
            _exit(2);
        }
    }

}



void oryza_model::model_rate() {

    atm.TMDA = (atm.TMMX+atm.TMMN)/2.;
    //problem ET2
    //ET2_initialize(crop.ANGA, crop.ANGB, crop.RDD, atm.TMDA, atm.VP, atm.WN, atm.latitude, IDOY, control.ETMOD, crop.CROPSTA, NL, crop.FAOF, crop.WL0, crop.WCLQT, crop.WCST, crop.LAI, crop.EVSC, crop.ETD, crop.TRC);
    
    if(control.PRODENV == "POTENTIAL"){
        WNOSTRESS(NL, crop.TRW, crop.TRWL, crop.LRSTRS, crop.LDSTRS, crop.LESTRS, crop.PCEW, crop.CPEW);
    }

    //call oryza
    oryza_rate();

    if(control.NITROENV == "POTENTIAL"){
        NNOSTRESS2_initialization( crop.NFLVI, crop.NMAXLT, crop.NFLVTB, DELT, crop.CROPSTA, crop.DVS, crop.WLVG, crop.LAI, crop.SLA, crop.NFLV, crop.NSLLV, crop.RNSTRS );
    }

    crop.NFLV = crop.NFLVI;
    crop.NSLLV = 1.;
    crop.RNSTRS = 1.;


    if(control.PRODENV == "POTENTIAL"){
        crop.TKLT = 100.;
        crop.ZRTMS = 100.;
        crop.WL0 = 0.;
        NL = crop.NLXM;
        for(int i = 0; i < NL; i++){
            crop.WCLQT[i] = 0.3;
            crop.WCST[i] = 0.3;
        }
    }




}



void oryza_model::model_state(){
    //-------Summation of some state variables
    crop.RAINCU = crop.RAINCU + atm.RAIN;


    if (crop.CROPSTA  ==  3) crop.CROPSTA = 4;

    if (crop.CROPSTA  ==  2)  {
        if (crop.DAE  ==  double(crop.SBDUR)) crop.CROPSTA = 3;
// 	 GIS. direct seeded. NO TRANSPLANTING SHOCK
        if (control.ESTAB == "DIRECT-SEED") crop.CROPSTA = 4;
    }
    if (crop.CROPSTA  ==  1) {
        if (control.ESTAB == "TRANSPLANT") {
            crop.CROPSTA = 2;
        } else if (control.ESTAB == "DIRECT-SEED") {
//               CROPSTA = 4
// GIS to allow initialization of the water balance.
            crop.CROPSTA = 2;
        }
    }

    if (crop.CROPSTA  ==  0)  {
        //problem
        //crop.IDATE = DTFSECMP(crop.EMYR, crop.EMD, IYEAR, IDOY);
        if (crop.IDATE  ==  0)  {
            crop.CROPSTA = 1;
        } else if (crop.IDATE  ==  1)  {
            //problem error message
            _exit(2);
        }
    }
}


