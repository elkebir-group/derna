//
// Created by Summer Gu on 8/16/22.
//

#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "default.h"
#include "params/constants.h"
#include "utils.h"
#include "params/intl11.h"
#include "params/intl11dH.h"
#include "params/intl21.h"
#include "params/intl21dH.h"
#include "params/intl22.h"
#include "params/intl22dH.h"

using namespace std;

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */
PUBLIC double lxc37=107.856;
PUBLIC int ML_intern37=-90;
PUBLIC int ML_interndH=-220;
PUBLIC int ML_closing37=930;
PUBLIC int ML_closingdH=3000;
PUBLIC int ML_BASE37=0;
PUBLIC int ML_BASEdH=0;
PUBLIC int MAX_NINIO=300;
PUBLIC int ninio37=60;
PUBLIC int niniodH=320;
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=370;


// loop length penalty based on loop type
PUBLIC int hairpin37[31] = {   INF,   INF,   INF,   540,   560,   570,   540,   600,   550,   640,   650,   660,   670,   680,   690,   690,   700,   710,   710,   720,   720,   730,   730,   740,   740,   750,   750,   750,   760,   760,   770};
PUBLIC int hairpindH[31] = {   INF,   INF,   INF,   130,   480,   360,  -290,   130,  -290,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500};
PUBLIC int bulge37[31] = {   INF,   380,   280,   320,   360,   400,   440,   460,   470,   480,   490,   500,   510,   520,   530,   540,   540,   550,   550,   560,   570,   570,   580,   580,   580,   590,   590,   600,   600,   600,   610};
PUBLIC int bulgedH[31] = {   INF,  1060,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710};
PUBLIC int internal_loop37[31] = {   INF,   INF,   100,   100,   110,   200,   200,   210,   230,   240,   250,   260,   270,   280,   290,   290,   300,   310,   310,   320,   330,   330,   340,   340,   350,   350,   350,   360,   360,   370,   370};
PUBLIC int internal_loopdH[31] = {   INF,   INF,   -720,   -720,  -720,  -680,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130};

// special case hairpin loop energies
PUBLIC vector<string> TriloopSeq = {
        "CAACG",
        "GUUAC"
};
PUBLIC int Triloop37[2] = {   680,   690};
PUBLIC int TriloopdH[2] = {  2370,  1080};

PUBLIC vector<string> TetraloopSeq = {
        "CAACGG",
        "CCAAGG",
        "CCACGG",
        "CCCAGG",
        "CCGAGG",
        "CCGCGG",
        "CCUAGG",
        "CCUCGG",
        "CUAAGG",
        "CUACGG",
        "CUCAGG",
        "CUCCGG",
        "CUGCGG",
        "CUUAGG",
        "CUUCGG",
        "CUUUGG"
};
PUBLIC int Tetraloop37[16] = {   550,   330,   370,   340,   350,   360,   370,   250,   360,   280,   370,   270,   280,   350,   370,   370};
PUBLIC int TetraloopdH[16] = {   690, -1030,  -330,  -890,  -660,  -750,  -350, -1390,  -760, -1070,  -660, -1290, -1070,  -620, -1530,  -680};

PUBLIC vector<string> HexaloopSeq = {
        "ACAGUACU",
        "ACAGUGAU",
        "ACAGUGCU",
        "ACAGUGUU",
};
PUBLIC int Hexaloop37[4] = {   280,   360,   290,   180};
PUBLIC int HexaloopdH[4] = { -1680, -1140, -1280, -1540};

PUBLIC int BP_pair[5][5] =
/* _  A  C  G  U  */
        { { 0, 0, 0, 0, 0 },
          { 0, 0, 0, 0, 5 },
          { 0, 0, 0, 1, 0 },
          { 0, 0, 2, 0, 3 },
          {0, 6, 0, 4, 0 }
        };
PUBLIC int rtype[7] = { 0, 2, 1, 4, 3, 6, 5 };

PUBLIC vector<char> to_char = {
        'A','C','G','U'
};

PUBLIC int internal_loop_map[3][3] = {
        {1,2,0},
        {3,4,5},
        {0,6,0},
};

PUBLIC vector<int> n_codon {4, 6, 2, 2, 2, 2, 2, 4, 2, 1, 3, 6, 2, 2, 4, 6, 4, 1, 2, 4};

PUBLIC bool codonMatch[20][3][3][4][4] = {
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, true, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, true, true, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, true, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, true, false}, {false, false, true, false}, {false, false, true, false}, {false, false, true, false}} ,
                        {{false, true, false, false}, {false, true, false, false}, {false, true, false, false}, {false, true, false, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, true, false}, {false, false, true, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, true, false}, {true, true, true, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {true, true, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, true, true, true}, {false, false, false, false}} ,
                },
                {
                        {{true, true, false, false}, {false, true, false, false}, {true, true, false, false}, {false, true, false, false}} ,
                        {{false, false, true, false}, {false, false, true, false}, {false, false, true, false}, {false, false, true, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, true, false, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, true, false}, {false, false, false, false}, {false, false, true, false}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, true, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, true}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, true}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, true, false, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, true, false}, {false, false, false, false}, {false, false, true, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, false, true, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, true, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, true, false, false}, {false, false, false, false}, {false, true, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, false, true, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, true, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, true, true, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, true, true, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, true, false}, {false, false, true, false}, {false, false, true, false}, {false, false, true, false}} ,
                        {{false, false, true, false}, {false, false, true, false}, {false, false, true, false}, {false, false, true, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, true, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, true, false, false}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, true, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, true}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, true, false, true}} ,
                },
                {
                        {{true, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, true}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, false, false, false}, {true, false, true, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, true, true, true}} ,
                },
                {
                        {{false, true, false, true}, {false, true, false, false}, {false, true, false, true}, {false, true, false, false}} ,
                        {{false, false, false, true}, {false, false, false, true}, {false, false, false, true}, {false, false, false, true}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{true, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {true, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, true}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, true}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, true, false, false}, {false, true, false, false}, {false, true, false, false}, {false, true, false, false}} ,
                        {{false, true, false, false}, {false, true, false, false}, {false, true, false, false}, {false, true, false, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, true, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {true, true, true, true}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, true}, {true, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, true, false, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, true}, {true, false, false, true}, {false, false, false, true}, {true, false, false, true}} ,
                        {{false, true, false, false}, {false, true, true, false}, {false, true, false, false}, {false, true, true, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, true, true, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {true, true, true, true}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{true, false, false, false}, {true, false, false, false}, {true, false, false, false}, {true, false, false, false}} ,
                        {{false, true, false, false}, {false, true, false, false}, {false, true, false, false}, {false, true, false, false}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, true, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, true, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, true}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, true}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, true, false, true}} ,
                },
                {
                        {{false, false, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{true, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                        {{false, true, false, true}, {false, false, false, false}, {false, false, false, false}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, true}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {true, false, false, false}, {false, false, false, false}, {true, false, false, false}} ,
                        {{false, false, false, false}, {false, true, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                },
        },
        {
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, true, false}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, true}, {false, false, false, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {true, true, true, true}, {false, false, false, false}} ,
                },
                {
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, true, false}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {false, false, false, true}} ,
                        {{false, false, false, false}, {false, false, false, false}, {false, false, false, false}, {true, true, true, true}} ,
                },
                {
                        {{false, false, true, false}, {false, false, true, false}, {false, false, true, false}, {false, false, true, false}} ,
                        {{false, false, false, true}, {false, false, false, true}, {false, false, false, true}, {false, false, false, true}} ,
                        {{true, false, false, false}, {false, true, false, false}, {false, false, true, false}, {false, false, false, true}} ,
                },
        },

};

PUBLIC int nucleotides[20][6][3] = {
        // A: GCU,GCC,GCA,GCG
        {{2,1,3},{2,1,1},{2,1,0},{2,1,2},{-1,-1,-1},{-1,-1,-1}},
        // R: CGU,CGC,CGA,CGG,AGA,AGG
        {{1,2,3},{1,2,1},{1,2,0},{1,2,2},{0,2,0},{0,2,2}},
        // N: AAU,AAC
        {{0,0,3},{0,0,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // D: GAU,GAC
        {{2,0,3},{2,0,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // C: UGU,UGC
        {{3,2,3},{3,2,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // Q: CAA,CAG
        {{1,0,0},{1,0,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // E: GAA,GAG
        {{2,0,0},{2,0,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // G: GGU,GGC,GGA,GGG
        {{2,2,3},{2,2,1},{2,2,0},{2,2,2},{-1,-1,-1},{-1,-1,-1}},
        // H: CAU,CAC
        {{1,0,3},{1,0,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // M: AUG
        {{0,3,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // I: AUU,AUC,AUA
        {{0,3,3},{0,3,1},{0,3,0},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // L: CUU,CUC,CUA,CUG,UUA,UUG
        {{1,3,3},{1,3,1},{1,3,0},{1,3,2},{3,3,0},{3,3,2}},
        // K: AAA,AAG
        {{0,0,0},{0,0,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // F: UUU,UUC
        {{3,3,3},{3,3,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // P: CCU,CCC,CCA,CCG
        {{1,1,3},{1,1,1},{1,1,0},{1,1,2},{-1,-1,-1},{-1,-1,-1}},
        // S: UCU,UCC,UCA,UCG,AGU,AGC
        {{3,1,3},{3,1,1},{3,1,0},{3,1,2},{0,2,3},{0,2,1}},
        // T: ACU,ACC,ACA,ACG
        {{0,1,3},{0,1,1},{0,1,0},{0,1,2},{-1,-1,-1},{-1,-1,-1}},
        // W: UGG
        {{3,2,2},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // Y: UAU,UAC
        {{3,0,3},{3,0,1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
        // V: GUU,GUC,GUA,GUG
        {{2,3,3},{2,3,1},{2,3,0},{2,3,2},{-1,-1,-1},{-1,-1,-1}}
};


PUBLIC double codon_cai[20][6];


PUBLIC double codon_cai_s[20][6];



PUBLIC double max_codon_cai[20][6][16];


// amino acid A = 0, R = 1, N = 2, D = 3, C = 4, Q = 5, E = 6, G = 7, H = 8, M = 9, I = 10, L = 11, K = 12, F = 13, P = 14
// S = 15, T = 16, W = 17, Y = 18, V = 19


double codon_usage[20][6] = {
        {0.26,0.4,0.23,0.11,0,0},
// R: CGU,CGC,CGA,CGG,AGA,AGG
        {0.08,0.19,0.11,0.21,0.2,0.2},
// N: AAU,AAC
        {0.46,0.54,0,0,0,0},
// D: GAU,GAC
        {0.46,0.54,0,0,0,0},
// C: UGU,UGC
        {0.45,0.55,0,0,0,0},
// Q: CAA,CAG
        {0.25,0.75,0,0,0,0},
// E: GAA,GAG
        {0.42,0.58,0,0,0,0},
// G: GGU,GGC,GGA,GGG
        {0.16,0.34,0.25,0.25,0,0},
// H: CAU,CAC
        {0.41,0.59,0,0,0,0},
// M: AUG
        {1.0,0,0,0,0,0},
// I: AUU,AUC,AUA
        {0.36,0.48,0.16,0,0,0},
// L: CUU,CUC,CUA,CUG,UUA,UUG
        {0.13,0.2,0.07,0.41,0.07,0.13},
// K: AAA,AAG
        {0.42,0.58,0,0,0,0},
// F: UUU,UUC
        {0.45,0.55,0,0,0,0},
// P: CCU,CCC,CCA,CCG
        {0.28,0.33,0.27,0.11,0,0},
// S: UCU,UCC,UCA,UCG,AGU,AGC
        {0.18,0.22,0.15,0.06,0.15,0.24},
// T: ACU,ACC,ACA,ACG
        {0.24,0.36,0.28,0.12,0,0},
// W: UGG
        {1.0,0,0,0,0,0},
// Y: UAU,UAC
        {0.43,0.57,0,0,0,0},
// V: GUU,GUC,GUA,GUG
        {0.18,0.24,0.11,0.47,0,0}
};



PUBLIC int max_cai_pos[20] = {1,3,1,1,1,1,1,1,1,0,1,3,1,1,1,5,1,0,1,3};

// (a == 0 && b == 3) || (a == 3 && b == 0) || (a == 2 && b == 3) || (a == 3 && b == 2)
PUBLIC int AU[4][4] = {
        {0,0,0,1},
        {0,0,0,0},
        {0,0,0,1},
        {1,0,1,0}
};

PUBLIC unordered_map<string, int> pair2pos = {
        {"NP", 0},
        {"CG", 1},
        {"GC", 2},
        {"GU", 3},
        {"UG", 4},
        {"AU", 5},
        {"UA", 6},
        {"NN", 7}
};


// 0: "NP", 1: "CG", 2: "GC", 3: "GU",4: "UG",5: "AU",6: "UA", 7: "NN"
PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
        {{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
                ,{   INF,  -240,  -330,  -210,  -140,  -210,  -210,  -140}
                ,{   INF,  -330,  -340,  -250,  -150,  -220,  -240,  -150}
                ,{   INF,  -210,  -250,   130,   -50,  -140,  -130,   130}
                ,{   INF,  -140,  -150,   -50,    30,   -60,  -100,    30}
                ,{   INF,  -210,  -220,  -140,   -60,  -110,   -90,   -60}
                ,{   INF,  -210,  -240,  -130,  -100,   -90,  -130,   -90}
                ,{   INF,  -140,  -150,   130,    30,   -60,   -90,   130}};
PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] =
        {{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
                ,{   INF, -1060, -1340, -1210,  -560, -1050, -1040,  -560}
                ,{   INF, -1340, -1490, -1260,  -830, -1140, -1240,  -830}
                ,{   INF, -1210, -1260, -1460, -1350,  -880, -1280,  -880}
                ,{   INF,  -560,  -830, -1350,  -930,  -320,  -700,  -320}
                ,{   INF, -1050, -1140,  -880,  -320,  -940,  -680,  -320}
                ,{   INF, -1040, -1240, -1280,  -700,  -680,  -770,  -680}
                ,{   INF,  -560,  -830,  -880,  -320,  -320,  -680,  -320}};



// 0: "NP", 1: "CG", 2: "GC", 3: "GU",4: "UG",5: "AU",6: "UA", 7: "NN"
PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,   -80,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,  -100,     0,  -100,     0}
                 ,{     0,     0,     0,     0,   -60}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,   -80,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,  -100,     0,  -100,     0}
                 ,{     0,     0,     0,     0,   -60}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,   -10,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -30,    70,   -30,    70}
                 ,{    70,    70,    70,    70,    10}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,   -10,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -30,    70,   -30,    70}
                 ,{    70,    70,    70,    70,    10}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,   -10,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -30,    70,   -30,    70}
                 ,{    70,    70,    70,    70,    10}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,   -10,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -30,    70,   -30,    70}
                 ,{    70,    70,    70,    70,    10}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,   -10,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -30,    70,   -30,    70}
                 ,{    70,    70,    70,    70,    10}
         }};
PUBLIC int mismatchIdH[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{   280,     0,     0,   280,     0}
                 ,{     0,     0,     0,  -340,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{   280,  -760,     0,   280,     0}
                 ,{     0,     0,     0,     0,  -580}
         }
                ,{{   280,     0,     0,   280,     0}
                 ,{     0,     0,     0,  -340,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{   280,  -760,     0,   280,     0}
                 ,{     0,     0,     0,     0,  -580}
         }
                ,{{   790,   500,   500,   790,   500}
                 ,{   500,   500,   500,   170,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   790,  -260,   500,   790,   500}
                 ,{   500,   500,   500,   500,   -80}
         }
                ,{{   790,   500,   500,   790,   500}
                 ,{   500,   500,   500,   170,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   790,  -260,   500,   790,   500}
                 ,{   500,   500,   500,   500,   -80}
         }
                ,{{   790,   500,   500,   790,   500}
                 ,{   500,   500,   500,   170,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   790,  -260,   500,   790,   500}
                 ,{   500,   500,   500,   500,   -80}
         }
                ,{{   790,   500,   500,   790,   500}
                 ,{   500,   500,   500,   170,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   790,  -260,   500,   790,   500}
                 ,{   500,   500,   500,   500,   -80}
         }
                ,{{   790,   500,   500,   790,   500}
                 ,{   500,   500,   500,   170,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   790,  -260,   500,   790,   500}
                 ,{   500,   500,   500,   500,   -80}
         }};

PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{   -80,  -100,  -110,  -100,   -80}
                 ,{  -140,  -150,  -150,  -140,  -150}
                 ,{   -80,  -100,  -110,  -100,   -80}
                 ,{  -150,  -230,  -150,  -240,  -150}
                 ,{  -100,  -100,  -140,  -100,  -210}
         }
                ,{{   -50,  -110,   -70,  -110,   -50}
                 ,{  -110,  -110,  -150,  -130,  -150}
                 ,{   -50,  -110,   -70,  -110,   -50}
                 ,{  -150,  -250,  -150,  -220,  -150}
                 ,{  -100,  -110,  -100,  -110,  -160}
         }
                ,{{    20,    20,   -20,   -10,   -20}
                 ,{    20,    20,   -50,   -30,   -50}
                 ,{   -10,   -10,   -20,   -10,   -20}
                 ,{   -50,  -100,   -50,  -110,   -50}
                 ,{   -10,   -10,   -30,   -10,  -100}
         }
                ,{{     0,   -20,   -10,   -20,     0}
                 ,{   -30,   -50,   -30,   -60,   -30}
                 ,{     0,   -20,   -10,   -20,     0}
                 ,{   -30,   -90,   -30,  -110,   -30}
                 ,{   -10,   -20,   -10,   -20,   -90}
         }
                ,{{   -10,   -10,   -20,   -10,   -20}
                 ,{   -30,   -30,   -50,   -30,   -50}
                 ,{   -10,   -10,   -20,   -10,   -20}
                 ,{   -50,  -120,   -50,  -110,   -50}
                 ,{   -10,   -10,   -30,   -10,  -120}
         }
                ,{{     0,   -20,   -10,   -20,     0}
                 ,{   -30,   -50,   -30,   -50,   -30}
                 ,{     0,   -20,   -10,   -20,     0}
                 ,{   -30,  -150,   -30,  -150,   -30}
                 ,{   -10,   -20,   -10,   -20,   -90}
         }
                ,{{    20,    20,   -10,   -10,     0}
                 ,{    20,    20,   -30,   -30,   -30}
                 ,{     0,   -10,   -10,   -10,     0}
                 ,{   -30,   -90,   -30,  -110,   -30}
                 ,{   -10,   -10,   -10,   -10,   -90}
         }};
PUBLIC int mismatchHdH[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{   560,  -570,   560,  -560,  -270}
                 ,{  -560,  -910,  -560,  -560,  -560}
                 ,{  -270,  -570,  -340,  -570,  -270}
                 ,{   560, -1400,   560,  -920,  -560}
                 ,{  -530,  -570,  -530,  -570, -1440}
         }
                ,{{    50,  -520,    50,  -560,  -400}
                 ,{  -400,  -520,  -400,  -560,  -400}
                 ,{    50,  -720,    50,  -720,  -420}
                 ,{  -400, -1290,  -400,  -620,  -400}
                 ,{   -30,  -720,   -30,  -720, -1080}
         }
                ,{{   970,   140,   970,   140,   570}
                 ,{   570,    30,   570,    20,   570}
                 ,{   970,   140,   970,   140,   340}
                 ,{   570,  -270,   570,    20,   570}
                 ,{   830,   140,   830,   140,   -50}
         }
                ,{{   230,   100,   230,   220,   190}
                 ,{  -110,  -110,  -260,  -520,  -260}
                 ,{   190,   -60,  -140,   -60,   190}
                 ,{   220,   100,  -260,   220,  -260}
                 ,{   230,   -60,   230,   -60,   -70}
         }
                ,{{   970,   140,   970,   140,   570}
                 ,{   570,   -20,   570,    20,   570}
                 ,{   970,   140,   970,   140,   340}
                 ,{   570,  -520,   570,    20,   570}
                 ,{   830,   140,   830,   140,  -380}
         }
                ,{{   230,   -30,   230,   -60,   190}
                 ,{   -30,   -30,  -260,  -520,  -260}
                 ,{   190,   -60,  -140,   -60,   190}
                 ,{  -260,  -590,  -260,  -520,  -260}
                 ,{   230,   -60,   230,   -60,   -70}
         }
                ,{{   970,   140,   970,   220,   570}
                 ,{   570,    30,   570,    20,   570}
                 ,{   970,   140,   970,   140,   340}
                 ,{   570,   100,   570,   220,   570}
                 ,{   830,   140,   830,   140,   -50}
         }};

/* mismatch_multi */
// NP: 0, CG: 1, GC:2, GU:3,UG:4,AU:5,UA:6, NN:7
PUBLIC int mismatchM37[NBPAIRS+1][5][5] =
        {{ /* NP.. */
                 {   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         },
         { /* CG.. */
                 {   -50,  -110,   -50,  -140,   -70}
                 ,{  -110,  -110,  -110,  -160,  -110}
                 ,{   -70,  -150,   -70,  -150,  -100}
                 ,{  -110,  -130,  -110,  -140,  -110}
                 ,{   -50,  -150,   -50,  -150,   -70}
         },
         { /* GC.. */
                 {   -80,  -140,   -80,  -140,  -100}
                 ,{  -100,  -150,  -100,  -140,  -100}
                 ,{  -110,  -150,  -110,  -150,  -140}
                 ,{  -100,  -140,  -100,  -160,  -100}
                 ,{   -80,  -150,   -80,  -150,  -120}
         },
         { /* GU.. */
                 {   -50,   -80,   -50,   -50,   -50}
                 ,{   -50,  -100,   -70,   -50,   -70}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -110,   -70,   -80,   -70}
                 ,{   -50,   -80,   -50,   -80,   -50}
         },
         { /* UG.. */
                 {   -30,   -30,   -60,   -60,   -60}
                 ,{   -30,   -30,   -60,   -60,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,  -100,   -70,  -100,   -60}
         },
         { /* AU.. */
                 {   -50,   -80,   -50,   -80,   -50}
                 ,{   -70,  -100,   -70,  -110,   -70}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -110,   -70,  -120,   -70}
                 ,{   -50,   -80,   -50,   -80,   -50}
         },
         { /* UA.. */
                 {   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
         },
         { /* NN.. */
                 {   -30,   -30,   -50,   -50,   -50}
                 ,{   -30,   -30,   -60,   -50,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -50,   -80,   -50,   -80,   -50}
         }};

/* mismatch_multi_enthalpies */
PUBLIC int mismatchMdH[NBPAIRS+1][5][5] =
        {{ /* NP.. */
                 {   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         },
         { /* CG.. */
                 {    50,  -400,    50,  -400,   -30}
                 ,{  -520,  -520,  -720,  -710,  -720}
                 ,{    50,  -400,    50,  -400,   -30}
                 ,{  -560,  -560,  -720,  -620,  -720}
                 ,{  -400,  -400,  -420,  -400,  -500}
         },
         { /* GC.. */
                 {  -270,  -560,  -270,  -560,  -530}
                 ,{  -570,  -910,  -570,  -820,  -570}
                 ,{  -340,  -560,  -340,  -560,  -530}
                 ,{  -560,  -560,  -570,  -920,  -570}
                 ,{  -270,  -560,  -270,  -560,  -860}
         },
         { /* GU.. */
                 {   310,  -480,  -180,   310,   140}
                 ,{   310,  -480,  -430,   310,  -430}
                 ,{  -140,  -630,  -510,  -630,  -140}
                 ,{  -150,  -890,  -430,  -150,  -430}
                 ,{   140,  -630,  -180,  -630,   140}
         },
         { /* UG.. */
                 {   600,   200,   600,   200,   460}
                 ,{   -60,  -340,  -230,   -60,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -230,  -350,  -230,  -350,  -230}
                 ,{   200,   200,   -30,   200,   160}
         },
         { /* AU.. */
                 {   140,  -400,  -180,  -380,   140}
                 ,{  -380,  -400,  -430,  -380,  -430}
                 ,{  -140,  -630,  -510,  -630,  -140}
                 ,{  -430,  -890,  -430,  -890,  -430}
                 ,{   140,  -630,  -180,  -630,   140}
         },
         { /* UA.. */
                 {   600,   200,   600,   200,   460}
                 ,{  -230,  -390,  -230,  -310,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -230,  -350,  -230,  -350,  -230}
                 ,{   200,   200,   -30,   200,  -170}
         },
         { /* NN.. */
                 {   600,   200,   600,   310,   460}
                 ,{   310,  -340,  -230,   310,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -150,  -350,  -230,  -150,  -230}
                 ,{   200,   200,   -30,   200,   160}
         }};

PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
         }};
PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
         }};

PUBLIC int mismatch23I37[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,   -50,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,  -110,     0,   -70,     0}
                 ,{     0,     0,     0,     0,   -30}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,  -120,     0,   -70,     0}
                 ,{     0,     0,     0,     0,   -30}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -40,    70,     0,    70}
                 ,{    70,    70,    70,    70,    40}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    20,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -40,    70,     0,    70}
                 ,{    70,    70,    70,    70,    40}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -40,    70,     0,    70}
                 ,{    70,    70,    70,    70,    40}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    20,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -40,    70,     0,    70}
                 ,{    70,    70,    70,    70,    40}
         }
                ,{{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,    70,    70,    70,    70}
                 ,{    70,   -40,    70,     0,    70}
                 ,{    70,    70,    70,    70,    40}
         }};
PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] =
        {{{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,  -570,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,  -860,     0,  -900,     0}
                 ,{     0,     0,     0,     0,  -640}
         }
                ,{{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0,     0,     0,     0,     0}
                 ,{     0, -1090,     0,  -900,     0}
                 ,{     0,     0,     0,     0,  -640}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,  -580,   500,  -400,   500}
                 ,{   500,   500,   500,   500,  -140}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   -60,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,  -360,   500,  -400,   500}
                 ,{   500,   500,   500,   500,  -140}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,  -580,   500,  -400,   500}
                 ,{   500,   500,   500,   500,  -140}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   -60,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,  -360,   500,  -400,   500}
                 ,{   500,   500,   500,   500,  -140}
         }
                ,{{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,   500,   500,   500,   500}
                 ,{   500,  -360,   500,  -400,   500}
                 ,{   500,   500,   500,   500,  -140}
         }};

/* mismatch_exterior */
PUBLIC int mismatchExt37[NBPAIRS+1][5][5] =
        {{ /* NP.. */
                 {   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         },
         { /* CG.. */
                 {   -50,  -110,   -50,  -140,   -70}
                 ,{  -110,  -110,  -110,  -160,  -110}
                 ,{   -70,  -150,   -70,  -150,  -100}
                 ,{  -110,  -130,  -110,  -140,  -110}
                 ,{   -50,  -150,   -50,  -150,   -70}
         },
         { /* GC.. */
                 {   -80,  -140,   -80,  -140,  -100}
                 ,{  -100,  -150,  -100,  -140,  -100}
                 ,{  -110,  -150,  -110,  -150,  -140}
                 ,{  -100,  -140,  -100,  -160,  -100}
                 ,{   -80,  -150,   -80,  -150,  -120}
         },
         { /* GU.. */
                 {   -50,   -80,   -50,   -50,   -50}
                 ,{   -50,  -100,   -70,   -50,   -70}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -110,   -70,   -80,   -70}
                 ,{   -50,   -80,   -50,   -80,   -50}
         },
         { /* UG.. */
                 {   -30,   -30,   -60,   -60,   -60}
                 ,{   -30,   -30,   -60,   -60,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,  -100,   -70,  -100,   -60}
         },
         { /* AU.. */
                 {   -50,   -80,   -50,   -80,   -50}
                 ,{   -70,  -100,   -70,  -110,   -70}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -110,   -70,  -120,   -70}
                 ,{   -50,   -80,   -50,   -80,   -50}
         },
         { /* UA.. */
                 {   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -70,  -100,   -70,  -100,   -80}
         },
         { /* NN.. */
                 {   -30,   -30,   -50,   -50,   -50}
                 ,{   -30,   -30,   -60,   -50,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -60,   -80,   -60,   -80,   -60}
                 ,{   -50,   -80,   -50,   -80,   -50}
         }};

/* mismatch_exterior_enthalpies */
PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] =
        {{ /* NP.. */
                 {   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
                 ,{   INF,   INF,   INF,   INF,   INF}
         },
         { /* CG.. */
                 {    50,  -400,    50,  -400,   -30}
                 ,{  -520,  -520,  -720,  -710,  -720}
                 ,{    50,  -400,    50,  -400,   -30}
                 ,{  -560,  -560,  -720,  -620,  -720}
                 ,{  -400,  -400,  -420,  -400,  -500}
         },
         { /* GC.. */
                 {  -270,  -560,  -270,  -560,  -530}
                 ,{  -570,  -910,  -570,  -820,  -570}
                 ,{  -340,  -560,  -340,  -560,  -530}
                 ,{  -560,  -560,  -570,  -920,  -570}
                 ,{  -270,  -560,  -270,  -560,  -860}
         },
         { /* GU.. */
                 {   310,  -480,  -180,   310,   140}
                 ,{   310,  -480,  -430,   310,  -430}
                 ,{  -140,  -630,  -510,  -630,  -140}
                 ,{  -150,  -890,  -430,  -150,  -430}
                 ,{   140,  -630,  -180,  -630,   140}
         },
         { /* UG.. */
                 {   600,   200,   600,   200,   460}
                 ,{   -60,  -340,  -230,   -60,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -230,  -350,  -230,  -350,  -230}
                 ,{   200,   200,   -30,   200,   160}
         },
         { /* AU.. */
                 {   140,  -400,  -180,  -380,   140}
                 ,{  -380,  -400,  -430,  -380,  -430}
                 ,{  -140,  -630,  -510,  -630,  -140}
                 ,{  -430,  -890,  -430,  -890,  -430}
                 ,{   140,  -630,  -180,  -630,   140}
         },
         { /* UA.. */
                 {   600,   200,   600,   200,   460}
                 ,{  -230,  -390,  -230,  -310,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -230,  -350,  -230,  -350,  -230}
                 ,{   200,   200,   -30,   200,  -170}
         },
         { /* NN.. */
                 {   600,   200,   600,   310,   460}
                 ,{   310,  -340,  -230,   310,  -230}
                 ,{   600,   200,   600,   200,   460}
                 ,{  -150,  -350,  -230,  -150,  -230}
                 ,{   200,   200,   -30,   200,   160}
         }};

/* dangle5 */
PUBLIC int dangle5_37[NBPAIRS+1][5] =
        { /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -10,   -50,   -30,   -20,   -10},
/* GC */ {    -0,   -20,   -30,    -0,    -0},
/* GU */ {   -20,   -30,   -30,   -40,   -20},
/* UG */ {   -10,   -30,   -10,   -20,   -20},
/* AU */ {   -20,   -30,   -30,   -40,   -20},
/* UA */ {   -10,   -30,   -10,   -20,   -20},
/* NN */ {    -0,   -20,   -10,    -0,    -0}
        };

/* dangle3 */
PUBLIC int dangle3_37[NBPAIRS+1][5] =
        { /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -40,  -110,   -40,  -130,   -60},
/* GC */ {   -80,  -170,   -80,  -170,  -120},
/* GU */ {   -10,   -70,   -10,   -70,   -10},
/* UG */ {   -50,   -80,   -50,   -80,   -60},
/* AU */ {   -10,   -70,   -10,   -70,   -10},
/* UA */ {   -50,   -80,   -50,   -80,   -60},
/* NN */ {   -10,   -70,   -10,   -70,   -10}
        };

/* dangle5_enthalpies */
PUBLIC int dangle5_dH[NBPAIRS+1][5] =
        { /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   330,  -240,   330,    80,  -140},
/* GC */ {    70,  -160,    70,  -460,   -40},
/* GU */ {   310,   160,   220,    70,   310},
/* UG */ {   690,   -50,   690,    60,    60},
/* AU */ {   310,   160,   220,    70,   310},
/* UA */ {   690,   -50,   690,    60,    60},
/* NN */ {   690,   160,   690,    80,   310}
        };

/* dangle3_enthalpies */
PUBLIC int dangle3_dH[NBPAIRS+1][5] =
        { /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {  -280,  -740,  -280,  -640,  -360},
/* GC */ {  -410,  -900,  -410,  -860,  -750},
/* GU */ {   -70,  -570,   -70,  -580,  -220},
/* UG */ {   -90,  -490,   -90,  -550,  -230},
/* AU */ {   -70,  -570,   -70,  -580,  -220},
/* UA */ {   -90,  -490,   -90,  -550,  -230},
/* NN */ {   -70,  -490,   -70,  -550,  -220}
        };



PUBLIC double lxc;
PUBLIC int stackE[NBPAIRS+1][NBPAIRS+1];
PUBLIC int hairpins[31];
PUBLIC int bulge[31];
PUBLIC int internal_loop[31];
PUBLIC int mismatchI[NBPAIRS+1][5][5];
PUBLIC int mismatch1nI[NBPAIRS+1][5][5];
PUBLIC int mismatch23I[NBPAIRS+1][5][5];
PUBLIC int mismatchH[NBPAIRS+1][5][5];
PUBLIC int int11[NBPAIRS+1][NBPAIRS+1][5][5];
PUBLIC int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
PUBLIC int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
PUBLIC int ML_BASE;
PUBLIC int ML_closing;
PUBLIC int ML_intern;
PUBLIC int ninio;
PUBLIC int TerminalAU;
PUBLIC int Tetraloop[16];
PUBLIC int Triloop[2];
PUBLIC int Hexaloop[4];
PUBLIC unordered_map<string, int> hairpinE;
PUBLIC unordered_map<string, int> max_cai_map;

PUBLIC vector<double> Log(9000);

void fill_stack(const string &filename, int data, char delimeter) {
    ifstream file(filename);

//    cout << data << endl;

    string line, word;
    int idx = 0, cidx = 0;

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "Pairs") break;
            if (idx == 0) {
                cidx = pair2pos[word];
            } else {
                if (!data) {
//                    cout << idx << " " << word << endl;
                    stack37[cidx][idx-1] = stoi(word);
                } else {
                    stackdH[cidx][idx-1] = stoi(word);
                }

            }
            idx += 1;
        }

    }
}

void fill_mismatch(const string &filename, int data, char delimeter) {
    ifstream file(filename);
//    cout << data << endl;
    unordered_map<string, int> cnt = {
            {"NP", 0},
            {"CG", 0},
            {"GC", 0},
            {"GU", 0},
            {"UG", 0},
            {"AU", 0},
            {"UA", 0},
            {"NN", 0}
    };

    string line, word, pair;
    int idx = 0, cidx = 0;

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "Pairs") break;
            if (idx == 0) {
                cidx = pair2pos[word];
                pair = word;
            } else {
//                cout << idx << " " << word << endl;
                switch (data) {
                    case 0:
                        mismatchI37[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 1:
                        mismatchIdH[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 2:
                        mismatchH[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 3:
                        mismatchHdH[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 4:
                        mismatch1nI[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 5:
                        mismatch1nIdH[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 6:
                        mismatch23I[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    case 7:
                        mismatch23IdH[cidx][cnt[pair]][idx-1] = stoi(word);
                        break;
                    default:
                        break;

                }

//                codon_usage[cidx][idx-1] = stof(word);
            }
            idx += 1;
        }
        cnt[pair] += 1;

    }
}

void fill_intl11(const string &filename, int data, char delimeter) {
    ifstream file(filename);
//    cout << data << endl;
    unordered_map<string, int> cnt = {
            {"NP", 0},
            {"CG", 0},
            {"GC", 0},
            {"GU", 0},
            {"UG", 0},
            {"AU", 0},
            {"UA", 0},
            {"NN", 0}
    };

    string line, word, pair;
    int idx = 0, cidx = 0;

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "Pairs") break;
            if (idx == 0) {
                cidx = pair2pos[word];
                pair = word;
            } else {
                if (!data) {
                    int11_37[cidx][cnt[pair]/5][cnt[pair]%5][idx-1] = stoi(word);
                }
                else {
                    int11_dH[cidx][cnt[pair]/5][cnt[pair]%5][idx-1] = stoi(word);
                }

            }
            idx += 1;
        }
        cnt[pair] += 1;

    }
}

void fill_intl21(const string &filename, int data, char delimeter) {
    ifstream file(filename);
//    cout << data << endl;
    unordered_map<string, int> cnt = {
            {"NP", 0},
            {"CG", 0},
            {"GC", 0},
            {"GU", 0},
            {"UG", 0},
            {"AU", 0},
            {"UA", 0},
            {"NN", 0}
    };

    string line, word, pair;
    int idx = 0, cidx = 0;

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "Pairs") break;
            if (idx == 0) {
                cidx = pair2pos[word];
                pair = word;
            } else {
                if (!data) {
                    int21_37[cidx][cnt[pair]/25][(cnt[pair]/5)%5][cnt[pair]%5][idx-1] = stoi(word);
                }
                else {
                    int21_dH[cidx][cnt[pair]/25][(cnt[pair]/5)%5][cnt[pair]%5][idx-1] = stoi(word);
                }
            }
            idx += 1;
        }
        cnt[pair] += 1;
    }
}

void fill_intl22(const string &filename, int data, char delimeter) {
    ifstream file(filename);
//    cout << data << endl;
    unordered_map<string, int> cnt = {
            {"NP", 0},
            {"CG", 0},
            {"GC", 0},
            {"GU", 0},
            {"UG", 0},
            {"AU", 0},
            {"UA", 0},
            {"NN", 0}
    };

    string line, word, pair;
    int idx = 0, cidx = 0;

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "Pairs") break;
            if (idx == 0) {
                cidx = pair2pos[word];
                pair = word;
            } else {
                if (!data) {
                    int22_37[cidx][cnt[pair]/125][(cnt[pair]/25)%5][(cnt[pair]/5)%5][cnt[pair]%5][idx-1] = stoi(word);
                }
                else {
                    int22_dH[cidx][cnt[pair]/125][(cnt[pair]/25)%5][(cnt[pair]/5)%5][cnt[pair]%5][idx-1] = stoi(word);
                }
            }
            idx += 1;
        }
        cnt[pair] += 1;

    }
}

void fill_codon(const string &filename, char delimeter) {
    ifstream file(filename);

    string line, word;
    int idx = 0, cidx = 0;
    unordered_map<string, int> miscellaneous = {
            {"lxc37", 0},
            {"ML_intern37",1},
            {"ML_interndH", 2},
            {"ML_closing37",3},
            {"ML_closingdH",4},
            {"ML_BASE37", 5},
            {"ML_BASEdH", 6},
            {"MAX_NINIO", 7},
            {"ninio37", 8},
            {"niniodH", 9},
            {"TerminalAU37", 10},
            {"TerminalAUdH", 11}
    };

    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word == "0") break;
            if (idx == 0) {
                cidx = miscellaneous[word];
            } else {
                switch (cidx) {
                    case 0:
                        lxc37 = stod(word);
                        break;
                    case 1:
                        ML_intern37 = stoi(word);
                        break;
                    case 2:
                        ML_interndH = stoi(word);
                        break;
                    case 3:
                        ML_closing37 = stoi(word);
                        break;
                    case 4:
                        ML_closingdH = stoi(word);
                        break;
                    case 5:
                        ML_BASE37 = stoi(word);
                        break;
                    case 6:
                        ML_BASEdH = stoi(word);
                        break;
                    case 7:
                        MAX_NINIO = stoi(word);
                        break;
                    case 8:
                        ninio37 = stoi(word);
                        break;
                    case 9:
                        niniodH = stoi(word);
                        break;
                    case 10:
                        TerminalAU37 = stoi(word);
                        break;
                    case 11:
                        TerminalAUdH = stoi(word);
                        break;
                    default:
                        break;
                }
            }
            idx += 1;
        }
    }

}


void fill_miscellaneous(const string &filename, char delimeter) {
    ifstream file(filename);

    string line, word;
    int idx = 0, cidx = 0;
    while (getline(file, line)) {
        stringstream str(line);
        idx = 0;
        while(getline(str, word, delimeter)) {
            if (idx == 0 && word.empty()) break;
            if (idx == 0) {
                cidx = aa_index(word[0]);
            } else {
                codon_usage[cidx][idx-1] = stof(word);
            }
            idx += 1;
        }
    }

}

void scale_params(const string & file, const string & paramspath, double temp) {
    int  i, j, k, l;
    double tempf = (temp + K0) / Tmeasure;

    if (!paramspath.empty()) {
        string stack_file = PARAMSPATH + "stack.csv";
        string stackh_file = PARAMSPATH + "stack_H.csv";
        string mismatchI_file = PARAMSPATH + "mismatchI.csv";
        string mismatchIh_file = PARAMSPATH + "mismatchI_H.csv";
        string mismatchH_file = PARAMSPATH + "mismatchH.csv";
        string mismatchHh_file = PARAMSPATH + "mismatchH_H.csv";
        string mismatch1nI_file = PARAMSPATH + "mismatch1nI.csv";
        string mismatch1nIh_file = PARAMSPATH + "mismatch1nI_H.csv";
        string mismatch23I_file = PARAMSPATH + "mismatch23I.csv";
        string mismatch23Ih_file = PARAMSPATH + "mismatch23I_H.csv";
        string intl11_file = PARAMSPATH + "intl11.csv";
        string intl11h_file = PARAMSPATH + "intl11_H.csv";
        string intl21_file = PARAMSPATH + "intl21.csv";
        string intl21h_file = PARAMSPATH + "intl21_H.csv";
        string intl22_file = PARAMSPATH + "intl22.csv";
        string intl22h_file = PARAMSPATH + "intl22_H.csv";
        string miscellaneous_file = PARAMSPATH + "miscellaneous.csv";
        if (exists(stack_file)) fill_stack(stack_file, 0);
        if (exists(stackh_file)) fill_stack(stackh_file, 1);
        if (exists(mismatchI_file)) fill_mismatch(mismatchI_file, 0);
        if (exists(mismatchIh_file)) fill_mismatch(mismatchIh_file, 1);
        if (exists(mismatchH_file)) fill_mismatch(mismatchH_file, 2);
        if (exists(mismatchHh_file)) fill_mismatch(mismatchHh_file, 3);
        if (exists(mismatch1nI_file)) fill_mismatch(mismatch1nI_file, 4);
        if (exists(mismatch1nIh_file)) fill_mismatch(mismatch1nIh_file, 5);
        if (exists(mismatch23I_file)) fill_mismatch(mismatch23I_file, 6);
        if (exists(mismatch23Ih_file)) fill_mismatch(mismatch23Ih_file, 7);
        if (exists(intl11_file)) fill_intl11(intl11_file, 0);
        if (exists(intl11h_file)) fill_intl11(intl11h_file, 1);
        if (exists(intl21_file)) fill_intl21(intl21_file, 0);
        if (exists(intl21h_file)) fill_intl21(intl21h_file, 1);
        if (exists(intl22_file)) fill_intl22(intl22_file, 0);
        if (exists(intl22h_file)) fill_intl22(intl22h_file, 1);
        if (exists(miscellaneous_file)) fill_miscellaneous(miscellaneous_file);
    }


    if (!file.empty()) {
        fill_codon(file);
    }


    ninio              = RESCALE_dG(ninio37, niniodH, tempf);
    lxc                   = lxc37 * tempf;
    TerminalAU            = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
    ML_intern              = RESCALE_dG(ML_intern37, ML_interndH, tempf);
    ML_BASE                = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
    ML_closing             = RESCALE_dG(ML_closing37, ML_closingdH, tempf);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            if (AU[i][j]) AU[i][j] = TerminalAU;
        }
    }



    for (i = 0; i < 20; i++) {
        max_cai_pos[i] = max_element(codon_usage[i], codon_usage[i] + 6) - codon_usage[i];
        double max_c = codon_usage[i][max_cai_pos[i]];
        for (j = 0; j < 6; j++) {
            if (j < n_codon[i] && codon_usage[i][j] == 0) {
                codon_cai[i][j] = log((codon_usage[i][j] + EPSILON)/max_c); //*100.0
                codon_cai_s[i][j] = log((codon_usage[i][j] + EPSILON)/max_c);
            } else if (codon_usage[i][j] == 0) {
                codon_cai[i][j] = -INF; //*100.0
                codon_cai_s[i][j] = -INF;

            } else {
                codon_cai[i][j] = log(codon_usage[i][j]/max_c); //*100.0
                codon_cai_s[i][j] = log(codon_usage[i][j]/max_c);
            }
        }
    }

    int pairs[16][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3},{1,0},{2,0},{3,0},{2,1},{3,1},{3,2},{0,0},{1,1},{2,2},{3,3}};

    for (i = 0; i < 16; i++) {

        string s = {to_char[pairs[i][0]], to_char[pairs[i][1]]};
        max_cai_map.insert(make_pair(s, i));
    }


    for (i = 0; i < 20; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 16; k++) {
                double m = -INF;
                for (l = 0; l < 6; l++) {
                    if (nucleotides[i][l][j] == k) {
                        m = max(m, codon_cai[i][l]);
                    }
                }
                max_codon_cai[i][j][k] = m;
            }
        }
        for (j = 3; j < 6; j++) {
            for (k = 0; k < 16; k++) {
                double m = -INF;
                for (l = 0; l < 6; l++) {

                    switch (j) {

                        case 3:
                            if (nucleotides[i][l][0] == pairs[k][0] && nucleotides[i][l][1] == pairs[k][1]) {

                                m = max(m, codon_cai[i][l]);

                            }
                            break;
                        case 4:
                            if (nucleotides[i][l][0] == pairs[k][0] && nucleotides[i][l][2] == pairs[k][1]) {
                                m = max(m, codon_cai[i][l]);
                            }
                            break;
                        case 5:
                            if (nucleotides[i][l][1] == pairs[k][0] && nucleotides[i][l][2] == pairs[k][1]) {
                                m = max(m, codon_cai[i][l]);
                            }
                            break;
                        default:
                            break;

                    }
                }
                max_codon_cai[i][j][k] = m;
//                cout << i << " " << j << " " << k << " " << m << endl;
            }
        }
    }


    for (i = 0; i < (int)TriloopSeq.size(); i++)
        hairpinE.insert(make_pair(TriloopSeq[i], RESCALE_dG(Triloop37[i], TriloopdH[i], tempf)));

    for (i = 0; i < (int)TetraloopSeq.size(); i++)
        hairpinE.insert(make_pair(TetraloopSeq[i], RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf)));

    for (i = 0; i < (int)HexaloopSeq.size(); i++)
        hairpinE.insert(make_pair(HexaloopSeq[i], RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf)));

    for (i = 0; i < 9000; i++)
        Log[i] = log(i/30.);

    for (i = 0; i < 2; i++)
        Triloop[i] = RESCALE_dG(Triloop37[i], TriloopdH[i], tempf);

    for (i = 0; i < 16; i++)
        Tetraloop[i] = RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf);

    for (i = 0; i < 4; i++)
        Hexaloop[i] = RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf);

    for (i = 0; i < 31; i++)
        hairpins[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

    for (i = 0; i <= min(30, MAXLOOP); i++) {
        bulge[i]          = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
        internal_loop[i]  = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
    }

    for (; i <= MAXLOOP; i++) {
        bulge[i] = bulge[30] +
                           (int)(lxc * log((double)(i) / 30.));
        internal_loop[i] = internal_loop[30] +
                                   (int)(lxc * log((double)(i) / 30.));
    }

    for (i = 0; i <= NBPAIRS; i++)
        for (j = 0; j <= NBPAIRS; j++)
            stackE[i][j] = RESCALE_dG(stack37[i][j],
                                             stackdH[i][j],
                                             tempf);


    for (i = 0; i <= NBPAIRS; i++)
        for (j = 0; j < 5; j++)
            for (k = 0; k < 5; k++) {
                mismatchI[i][j][k] = RESCALE_dG(mismatchI37[i][j][k],
                                                        mismatchIdH[i][j][k],
                                                        tempf);
                mismatchH[i][j][k] = RESCALE_dG(mismatchH37[i][j][k],
                                                        mismatchHdH[i][j][k],
                                                        tempf);
                mismatch1nI[i][j][k] = RESCALE_dG(mismatch1nI37[i][j][k],
                                                          mismatch1nIdH[i][j][k],
                                                          tempf);
                mismatch23I[i][j][k] = RESCALE_dG(mismatch23I37[i][j][k],
                                                          mismatch23IdH[i][j][k],
                                                          tempf);
            }

    /* dangles */
//    for (i = 0; i <= NBPAIRS; i++)
//        for (j = 0; j < 5; j++) {
//            int dd;
//            dd = RESCALE_dG(dangle5_37[i][j],
//                            dangle5_dH[i][j],
//                            tempf);
//        }

    /* interior 1x1 loops */
    for (i = 0; i <= NBPAIRS; i++)
        for (j = 0; j <= NBPAIRS; j++)
            for (k = 0; k < 5; k++)
                for (l = 0; l < 5; l++)
                    int11[i][j][k][l] = RESCALE_dG(int11_37[i][j][k][l],
                                                           int11_dH[i][j][k][l],
                                                           tempf);

    /* interior 2x1 loops */
    for (i = 0; i <= NBPAIRS; i++)
        for (j = 0; j <= NBPAIRS; j++)
            for (k = 0; k < 5; k++)
                for (l = 0; l < 5; l++) {
                    int m;
                    for (m = 0; m < 5; m++)
                        int21[i][j][k][l][m] = RESCALE_dG(int21_37[i][j][k][l][m],
                                                                  int21_dH[i][j][k][l][m],
                                                                  tempf);
                }

    /* interior 2x2 loops */
    for (i = 0; i <= NBPAIRS; i++)
        for (j = 0; j <= NBPAIRS; j++)
            for (k = 0; k < 5; k++)
                for (l = 0; l < 5; l++) {
                    int m, n;
                    for (m = 0; m < 5; m++)
                        for (n = 0; n < 5; n++)
                            int22[i][j][k][l][m][n] = RESCALE_dG(int22_37[i][j][k][l][m][n],
                                                                         int22_dH[i][j][k][l][m][n],
                                                                         tempf);
                }

    

}