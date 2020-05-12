//
//  NCFE.hpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/5.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#ifndef NCFE_hpp
#define NCFE_hpp
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <fstream>
#include <stdio.h>
#include "pp.h"
#include "CubicFE.hpp"
#include "AIK.hpp"
using namespace NTL;
using namespace std;

struct FE_pp
{
    int n;
    int l;
    Cubic_pp Cpp;
    AIK_PP App;
    ZZ_pX Ring_P;
};

struct FE_MPK
{
    Cubic_MPK Cmpk;
};

void FE_pp_Print(FE_pp &para);
void FE_MPK_Print(FE_MPK &mpk);

void FE_SetUp(FE_pp &para, FE_MPK &mpk);
void FE_Enc(FE_pp &para, FE_MPK &mpk, Vec<ZZ_p> &x, Vec<ZZ_pE> &cubic_CT);
void FE_KeyGen(FE_pp &para, FE_MPK &mpk, int *Ax, Vec<ZZ_pE> *cubic_SK, Vec<ZZ_p> *F);
void FE_Dec(FE_pp &para, FE_MPK &mpk, Vec<ZZ_pE> &cubic_CT,  Vec<ZZ_pE> *cubic_SK, Vec<ZZ_p> *F, ZZ_p &ans);

#endif /* NCFE_hpp */
