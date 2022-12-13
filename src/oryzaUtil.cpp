#include "oryzaUtil.h"
#include <vector>
#include <algorithm>
#include <math.h>

// using namespace std;



double LIMIT(double min, double max, double v) {
    if (v < min) {
        v = min;
    } else if (v > max) {
        v = max;
    }
    return(v);
}

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

double NOTNUL(double x){
    if(x == 0.) return 1.;
    else return x;
}


double INSW(double x, double y1, double y2){
    if(x < 0.) return y1;
    else return y2;
}

double FCNSW(double x, double y1, double y2, double y3){
    if(x < 0.) return y1;
    else if(x == 0.) return y2;
    else return y3;
}

double REAAND(double x1, double x2){
    if(x1 > 0 && x2 > 0) return 1.;
    else return 0.;
}

double REANOR(double x1, double x2){
    if(x1 <= 0 && x2 <= 0) return 1.;
    else return 0.;
}