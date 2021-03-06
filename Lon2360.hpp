#ifndef ASU_LON2360
#define ASU_LON2360


/***********************************************************
 * This C++ template convert input longitude (in deg) to
 * [0,360] deg.
 *
 * input(s):
 * const double &lon  ----  Longitude.
 *
 * return(s):
 * double ans  ----  Longitude in 0 ~ 360.
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: coordinates, longitude.
***********************************************************/

inline double Lon2360(const double &lon){

    double ans=lon;

    if (0 <= ans && ans <= 360.0) {
        return ans;
    }
    else if (ans >= 0) {
        ans -= 360.0 * ((int)(ans / 360));
    }
    else {
        ans += 360.0 * (1 + (int)(-ans / 360));
    }
    return ans;
}

#endif
