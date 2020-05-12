//
//  CubicFE.hpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/4.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#ifndef CubicFE_hpp
#define CubicFE_hpp
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <stdio.h>
#include <fstream>
#include "pp.h"
#include "LinFE.hpp"
using namespace NTL;
using namespace std;

struct Cubic_pp
{
    int n;
    int Q;
    long p0;
    long p1;
    double sig0;
    double sig1;
    ZZ_pX Ring_P;
    Lin_pp Lin_pp;
};

struct Cubic_MPK
{
    Lin_PK Lmpk;
    Vec<ZZ_pE> u;
};

void Cubic_pp_Print(Cubic_pp &para);
void Cubic_MPK_Print(Cubic_MPK &mpk);

void Cubic_SetUp(Cubic_pp &para, Cubic_MPK &mpk, int &n, int &Q);

void Cubic_Enc(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &x, Vec<ZZ_pE> &CT);

void Cubic_KeyGen(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &g, int &k, Vec<ZZ_pE> &SK);

void Cubic_Dec(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &g, int &k,  Vec<ZZ_pE> &SK,  Vec<ZZ_pE> &CT, ZZ_pX &ans);

void Cubic_test();

#include <stdio.h>

#endif /* CubicFE_hpp */
