//
//  LinFE.cpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/4/28.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#include "LinFE.hpp"

void Lin_pp_Print(Lin_pp &para)
{
    cout << "----------------Lin_pp-------------" << endl;
    cout << "Lin_p0 = " << para.Lin_p0 <<endl;
    cout << "Lin_p1 = " << para.Lin_p1 <<endl;
    cout << "L      = " << para.L << endl;
    cout << "sig    = " << para.sig <<endl;
    cout << "Ring_P = " << para.Ring_P <<endl;
    cout << "-----------------------------------" << endl;
}
void Lin_PK_Print(Lin_PK &mpk)
{
    cout << "----------------Lin_PK-------------" << endl;
    cout << "Lin_u  = " << mpk.u <<endl;
    cout << "Lin_ome = " << mpk.ome <<endl;
    cout << "-----------------------------------" << endl;
}

void Turn_p1(Lin_pp &para)
{
    if(ZZ_pE::modulus() != Lin_p1)
    {
        ZZ_p :: init(ZZ(para.Lin_p1));
        ZZ_pE::init(para.Ring_P);
    }
}

void Turn_p0(Lin_pp &para)
{
    
    if(ZZ_pE::modulus() != Lin_p0)
    {
        ZZ_p :: init(ZZ(para.Lin_p0));
        ZZ_pE::init(para.Ring_P);
    }
}

void SampleD(Lin_pp &para, ZZ_pE &mu)
{
    ZZ_p :: init(ZZ(para.Lin_p1));
    ZZ_p d;
    
    bool iis =1;
    
    while(iis)
    {
        ZZ t = ZZ(3);
        t = para.Lin_p1 + para.sig*2*t;
        ZZ r1 = RandomBnd(t);
        if( r1 > para.Lin_p1)
        {
            bool is =1;
            while(is)
            {
                d = random_ZZ_p();
                if(rep(d)<= ZZ(3*para.sig) && rep(d) >= para.Lin_p1-1-(3*para.sig))
                {
                    is =0;
                }
            }
        }
        else{
            d = random_ZZ_p();
        }
        double rho;
        int tx =0;
        for(int i=0;i<para.Lin_p1;++i)
        {
            if (i == rep(d)) tx =i;
        }
        rho = exp((0.0-atan(1)*4*pow((double)tx, 2.0))/pow(para.sig, 2.0));
        
        int r2 = rand();
        if(r2 > rho) iis =0;
    }
    
    Turn_p1(para);
    mu = d;
}

void Lin_SetUp(Lin_pp &para, Lin_PK &mpk, int &L)
{
    para.Lin_p0 = Lin_p0;
    para.Lin_p1 = Lin_p1;
    para.sig = sig;
    para.L = L;
    
    ZZ_p :: init(ZZ(para.Lin_p1));
    para.Ring_P.SetLength(Poly_n);
    SetCoeff(para.Ring_P, 0, 1);
    SetCoeff(para.Ring_P, Poly_n, 1);
    ZZ_pE::init(para.Ring_P);
    
    Vec<ZZ_pE> uu, oo;
    uu.SetLength(L);
    oo.SetLength(2);
    
    for(int i=0;i<L;++i)uu[i] = random_ZZ_pE();
    oo[0] = random_ZZ_pE();
    oo[1] = random_ZZ_pE();
    
    
    mpk.u.SetLength(L);
    mpk.ome.SetLength(2);
    
    Turn_p0(para);
    for(int i=0;i<L;++i) mpk.u[i] = 0 + uu[i];
    mpk.ome[0] = 0 + oo[0];
    mpk.ome[1] = 0 + oo[1];
}

void Lin_Enc(Lin_pp &para, Lin_PK &mpk, Vec<ZZ_pE> &x, Vec<ZZ_pE> &CT)
{
    CT.SetLength(para.L + 2);
    Vec<ZZ_pE> mu,uu;
    ZZ_pE s;
    
    mu.SetLength(para.L + 2);
    uu.SetLength(para.L + 2);
    
    for(int i=0; i<para.L+2; ++i) SampleD(para, uu[i]);
    Turn_p0(para);
    for(int i=0; i<para.L+2; ++i) mu[i] = uu[i];
    s = random_ZZ_pE();
    
//    cout << "mu = "<< mu <<endl;
//    cout << "s = " << s <<endl;
    
    for(int i=0; i<para.L; ++i) CT[i] = mpk.u[i] * s + mu[i] * para.Lin_p1 + x[i];
    CT[para.L] = mpk.ome[0] * s + mu[para.L] * para.Lin_p1;
    CT[para.L + 1] = mpk.ome[1] * s + mu[para.L+1] * para.Lin_p1;
    
}

void Lin_KeyGen(Lin_pp &para, Lin_PK &mpk, Vec<ZZ_pE> &y, Vec<ZZ_pE> &SK)
{
    SK.SetLength(3);
    Turn_p0(para);
    for (int i=0; i< para.L; ++i) SK[0] += mpk.u[i]*y[i];
    
    Turn_p1(para);
    ZZ_pE kk1,kk2;
    kk1 = random_ZZ_pE();
    kk2 = (SK[0] - kk1 * mpk.ome[0]) / mpk.ome[1];
    
    Turn_p0(para);
    SK[1] = kk1;
    SK[2] = kk2;

}

void Lin_Dec(Lin_pp &para, Vec<ZZ_pE> &y, Vec<ZZ_pE> &CT, Vec<ZZ_pE> &SK, ZZ_pX &ans)
{
    Turn_p0(para);
    ZZ_pE temp;
    
    for (int i=0; i< para.L; ++i) temp += CT[i] * y[i];
    temp -= SK[1] * CT[para.L];
    temp -= SK[2] * CT[para.L+1];
    
    
    ZZ_p :: init(ZZ(para.Lin_p1));
    ZZ_pX tt;
    tt = rep(temp);
//    cout<< "tt = "<<tt<<endl;
    for(int i=0; i< deg(tt)+1; ++ i )
    {
        while(rep(tt[i]) > para.Lin_p1)  tt[i] = 0 +  tt[i];
    }
    ans = tt;
}

void Lin_test()
{
    Lin_pp para;
    Lin_PK mpk;
    int L = 3;
    Lin_SetUp(para, mpk, L);
    Lin_pp_Print(para);
    Lin_PK_Print(mpk);
    
    Turn_p1(para);
    Vec<ZZ_pE> x;
    x.SetLength(L);
    for(int i =0;i<L;++i) x[i] = random_ZZ_pE();
    cout<< "x = " << x << endl;
    
    Vec<ZZ_pE> y;
    y.SetLength(L);
    for(int i =0;i<L;++i) y[i] = random_ZZ_pE();
    cout<< "y = " << y << endl;
    
    Vec<ZZ_pE> CT;
    Lin_Enc(para, mpk, x, CT);
    cout<< "CT = " << CT << endl;
    
    
    Vec<ZZ_pE> SK;
    Lin_KeyGen(para, mpk, y, SK);
    cout<< "SK = " << SK <<endl;
    
    ZZ_pX ans;
    Lin_Dec(para, y, CT, SK, ans);
    cout << ans <<endl;
}

