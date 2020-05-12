//
//  LinFE.hpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/4/28.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#ifndef LinFE_hpp
#define LinFE_hpp
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <stdio.h>
#include "pp.h"
using namespace NTL;
using namespace std;

struct Lin_pp
{
    long Lin_p0;
    long Lin_p1;
    int L;
    double sig;
    ZZ_pX Ring_P;
};

struct Lin_PK
{
    Vec<ZZ_pE> u;
    Vec<ZZ_pE> ome;
};

void Lin_pp_Print(Lin_pp &para);
void Lin_PK_Print(Lin_PK &mpk);

void Lin_SetUp(Lin_pp &para, Lin_PK &mpk, int &L);

void Lin_Enc(Lin_pp &para, Lin_PK &mpk, Vec<ZZ_pE> &x, Vec<ZZ_pE> &CT);

void Lin_KeyGen(Lin_pp &para, Lin_PK &mpk, Vec<ZZ_pE> &y, Vec<ZZ_pE> &SK);

void Lin_Dec(Lin_pp &para, Vec<ZZ_pE> &y, Vec<ZZ_pE> &CT, Vec<ZZ_pE> &SK, ZZ_pX &ans);


void Lin_test();


#endif /* LinFE_hpp */
