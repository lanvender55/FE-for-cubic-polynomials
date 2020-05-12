//
//  main.cpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/4/28.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include "pp.h"
#include "NCFE.hpp"
#include <fstream>

using namespace std;
using namespace NTL;

int main(int argc, const char * argv[]) {


    FE_pp para;
    FE_MPK mpk;
    
    FE_SetUp(para, mpk);


    int Ax[AIK_l*AIK_l]={0, 1, -2, 0, 0,
          0, 0,  0, -2, 2,
          0, 0,  0, 0, -3,
          0, 0,  0, 0, 3,
          0, 0,  0, 0, 0};
    
    Vec<ZZ_pE> cubic_SK[para.App.num];
    Vec<ZZ_p> F[para.App.num];
    
    FE_KeyGen(para, mpk, Ax, cubic_SK, F);
    
    Vec<ZZ_p> x;
    x.SetLength(para.n);
    x[0]= 1;
    x[1]= 0;
    x[2]= 1;
    
    Vec<ZZ_pE> CT;
    FE_Enc(para, mpk, x, CT);


    ZZ_p ans;
    FE_Dec(para, mpk, CT, cubic_SK, F, ans);
    cout<<"ans = " << ans<<endl;
    
    return 0;
}

// Time example
/*
 clock_t startTime,endTime;
    startTime = clock();
 endTime = clock();
   cout << "The run time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000<< "ms" << endl;
  */
