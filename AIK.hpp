//
//  AIK.hpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/5.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#ifndef AIK_hpp
#define AIK_hpp

#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "pp.h"
using namespace NTL;
using namespace std;

struct AIK_PP
{
    int n;
    int l;
    int num;
    int omega;
    ZZ p;
};

void AIK_pp_Print(AIK_PP &para);

void Matrix_Print(int &len, int *A);

void AIK_SetUp(AIK_PP &para);
void AIK_Rebuild(AIK_PP &para, int *Ax, Vec<ZZ_p> *F);

#endif /* AIK_hpp */
