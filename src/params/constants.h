//
// Created by Summer Gu on 8/16/22.
//

#ifndef RNA_DESIGN_CONSTANTS_H
#define RNA_DESIGN_CONSTANTS_H
#include <limits>

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** Infinity as used in minimization routines */
#define inf 10000000 /* (INT_MAX/10) */
#define INF 10000000

#define EMAX (INF/10)
/** forbidden */
#define FORBIDDEN 9999
/** bonus contribution */
#define BONUS 10000
/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30


#define UNIT 100

#define EPSILON 0.00001

const string PARAMSPATH = {}; //"../data/InputFiles/";


#endif //RNA_DESIGN_CONSTANTS_H
