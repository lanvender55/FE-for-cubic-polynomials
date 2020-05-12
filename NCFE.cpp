//
//  NCFE.cpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/5.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#include "NCFE.hpp"

void FE_pp_Print(FE_pp &para)
{
    cout << "n = " << para.n << endl;
    cout << "l = " << para.l << endl;
    cout << "Ring_P = " << para.Ring_P << endl;
    Cubic_pp_Print(para.Cpp);
    AIK_pp_Print(para.App);
}

void FE_MPK_Print(FE_MPK &mpk)
{
    Cubic_MPK_Print(mpk.Cmpk);
}

void FE_SetUp(FE_pp &para, FE_MPK &mpk)
{
    para.n = AIK_n;
    para.l = AIK_l;
    AIK_SetUp(para.App);
    Cubic_SetUp(para.Cpp, mpk.Cmpk, para.App.omega, para.App.num);
    para.Ring_P = para.Cpp.Ring_P;
}

void FE_Enc(FE_pp &para, FE_MPK &mpk, Vec<ZZ_p> &x, Vec<ZZ_pE> &cubic_CT)
{
    ZZ_p::init(ZZ(para.Cpp.p0));
    Vec<ZZ_p> cubic_x;
    cubic_x.SetLength(para.App.omega);
    cubic_x[0] = 1;
    for(int i=0;i<para.n;++i) cubic_x[i+1] = x[i];
    
    for(int i=para.n+1;i<para.App.omega;++i) cubic_x[i] = random_ZZ_p();
    
    ZZ_p::init(ZZ(para.Cpp.p0));
    ZZ_pE :: init(para.Cpp.Ring_P);
    Vec<ZZ_pE> Cx;
    Cx.SetLength(para.App.omega);
    for(int i=0; i<para.App.omega; ++i) Cx[i] = cubic_x[i];
    
    Cubic_Enc(para.Cpp, mpk.Cmpk, Cx, cubic_CT);
    
}

void FE_KeyGen(FE_pp &para, FE_MPK &mpk, int *Ax, Vec<ZZ_pE> *cubic_SK, Vec<ZZ_p> *F)
{
    AIK_Rebuild(para.App, Ax, F);
    
    ZZ_p::init(ZZ(para.Cpp.p0));
    ZZ_pE :: init(para.Cpp.Ring_P);
    Vec<ZZ_pE> Cg[para.App.num];
    for(int i=0; i<para.App.num; ++i)
    {
        Cg[i].SetLength(F[i].length());
        for(int j=0; j<Cg[i].length(); ++j) Cg[i][j] = F[i][j];
        Cubic_KeyGen(para.Cpp, mpk.Cmpk, Cg[i], i, cubic_SK[i]);
    }

}

void FE_Dec(FE_pp &para, FE_MPK &mpk, Vec<ZZ_pE> &cubic_CT,  Vec<ZZ_pE> *cubic_SK, Vec<ZZ_p> *F, ZZ_p &ans)
{
    ZZ_p::init(ZZ(para.Cpp.p0));
    ZZ_pE :: init(para.Cpp.Ring_P);
    Vec<ZZ_pE> Cg[para.App.num];
    for(int i=0; i<para.App.num; ++i)
    {
        Cg[i].SetLength(F[i].length());
        for(int j=0; j<Cg[i].length(); ++j) Cg[i][j] = F[i][j];
    }
    
    Vec<ZZ_pX> tans;
    tans.SetLength(para.App.num);
    for(int i=0; i<para.App.num; ++i)
    {
        Cubic_Dec(para.Cpp, mpk.Cmpk, Cg[i], i, cubic_SK[i], cubic_CT, tans[i]);
    }
    
    ZZ_p :: init(ZZ(para.Cpp.p0));
    Mat<ZZ_p> M;
    M.SetDims(para.l-1, para.l-1);
    
    for(int i=2;i<=para.l-1;++i) M(i,i-1) = -1;
    
    int tt =0;
    for(int i=1;i<=para.l-1;++i)
    {
        for(int j=i;j<=para.l-1;++j)
        {
            M(i,j) = tans[tt][0];
            tt++;
        }
    }
    
    determinant(ans,M);
    
}
