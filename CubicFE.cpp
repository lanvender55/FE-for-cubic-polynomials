//
//  CubicFE.cpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/4.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#include "CubicFE.hpp"

void Cubic_pp_Print(Cubic_pp &para)
{
    cout << "----------------Cubic_pp--------------" << endl;
    cout << "Cubic_n = " << para.n <<endl;
    cout << "Cubic_Q = " << para.Q <<endl;
    cout << "Cubic_p0 = " << para.p0 <<endl;
    cout << "Cubic_p1= " << para.p1 <<endl;
    cout << "Cubic_sig0 = " << para.sig0 <<endl;
    cout << "Cubic_sig1 = " << para.sig1 <<endl;
    cout << "Cubic_Poly = " << para.Ring_P <<endl;
    Lin_pp_Print(para.Lin_pp);
    cout << "--------------------------------------" << endl;
}

void Cubic_MPK_Print(Cubic_MPK &mpk)
{
    cout << "----------------Cubic_MPK-------------" << endl;
    cout << "Cubic_u = " << mpk.u <<endl;
    Lin_PK_Print(mpk.Lmpk);
    cout << "--------------------------------------" << endl;
}

void Turn_p1(Cubic_pp &para)
{
    if(ZZ_pE::modulus() != Cubic_p1)
    {
        ZZ_p :: init(ZZ(para.p1));
        ZZ_pE::init(para.Ring_P);
    }
}

void Turn_p0(Cubic_pp &para)
{
    if(ZZ_pE::modulus() != Cubic_p0)
    {
        ZZ_p :: init(ZZ(para.p0));
        ZZ_pE::init(para.Ring_P);
    }
}

void SampleD(Cubic_pp &para, ZZ_pE &mu, double &sig)
{
    ZZ_p :: init(ZZ(para.p0));
    ZZ_p d;
    
    bool iis =1;
    
    while(iis)
    {
        ZZ t = ZZ(3);
        t = para.p0 + sig*2*t;
        ZZ r1 = RandomBnd(t);
        if( r1 > para.p0)
        {
            bool is =1;
            while(is)
            {
                d = random_ZZ_p();
                if(rep(d)<= ZZ(3*sig) && rep(d) >= para.p0-1-(3*sig))
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
        for(int i=0;i<para.p0;++i)
        {
            if (i == rep(d)) tx =i;
        }
        rho = exp((0.0-atan(1)*4*pow((double)tx, 2.0))/pow(sig, 2.0));
        
        int r2 = rand();
        if(r2 > rho) iis =0;
    }
    
    Turn_p1(para);
    mu = d;
}

void Cubic_SetUp(Cubic_pp &para, Cubic_MPK &mpk, int &n, int &Q)
{
    para.n = n;
    para.Q = Q;
    para.p0 = Cubic_p0;
    para.p1 = Cubic_p1;
    para.sig0 = sig0;
    para.sig1 = sig1;
    
    ZZ_p :: init(ZZ(para.p1));
    para.Ring_P.SetLength(Poly_n);
    SetCoeff(para.Ring_P, 0, 1);
    SetCoeff(para.Ring_P, Poly_n, 1);
    ZZ_pE::init(para.Ring_P);
    
    int L = 1 + 2*para.n + para.n*(para.n-1)/2 + para.Q;
    
    Lin_SetUp(para.Lin_pp, mpk.Lmpk, L);
    
    mpk.u.SetLength(para.n);
    for(int i=0; i<para.n; ++i) mpk.u[i] = random_ZZ_pE();
    
}

void Cubic_Enc(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &x, Vec<ZZ_pE> &CT)
{
    Vec<ZZ_pE> mu,uu;
    ZZ_pE s;
    mu.SetLength(para.n);
    uu.SetLength(para.n);
    
    Turn_p0(para);
    for(int i=0; i<para.n; ++i) SampleD(para, uu[i], para.sig0);
    Turn_p1(para);
    for(int i=0; i<para.n; ++i) mu[i] = uu[i];
    s = random_ZZ_pE();
    
    
    CT.SetLength(para.n+para.Lin_pp.L+2);
    
    for(int i=0; i<para.n; ++i) CT[i] = mpk.u[i] * s + para.p0 * mu[i] + x[i];
    
    
    Vec<ZZ_pE> Lin_CT;
    Vec<ZZ_pE> Lin_x;
    Lin_x.SetLength(para.Lin_pp.L);
    
    Lin_x[0] = power(s, 3);
    for(int i=0;i<para.n;++i) Lin_x[i+1] = CT[i]*s*s;
    int temp = para.n+1;
    for(int i=0;i<para.n;++i)
    {
        for(int j=i;j<para.n;++j)
        {
            Lin_x[temp] = CT[i]*CT[j]*s;
            temp++;
        }
    }
    
    for(int i=0;i<para.Q;++i)
    {
        ZZ_pE temp1;
        SampleD(para,temp1,para.sig1);
        Turn_p1(para);
        temp1 = para.p0*temp1;
        Lin_x[temp] = temp1;
        temp++;
    }
    
    Lin_Enc(para.Lin_pp, mpk.Lmpk, Lin_x, Lin_CT);
    for(int i = para.n; i<para.n+para.Lin_pp.L+2; ++i) CT[i] = Lin_CT[i-para.n];
}

void Cubic_Gen_Lin_y(Cubic_pp &para, Cubic_MPK &mpk, int &k, Vec<ZZ_pE> &g,  Vec<ZZ_pE> &Lin_y)
{
    Lin_y.SetLength(para.Lin_pp.L);
    set(Lin_y[para.Lin_pp.L-(para.Q-k)]);
    
    Turn_p1(para);
    int ijo = 0;
    for(int i=0;i<para.n;++i)
    {
        for(int j=i;j<para.n;++j)
        {
            for(int o=j;o<para.n;++o)
            {
                Lin_y[0] += g[ijo]*(-1)*mpk.u[i]*mpk.u[j]*mpk.u[o];
                
                Lin_y[1+i] += g[ijo] * mpk.u[j] * mpk.u[o];
                Lin_y[1+j] += g[ijo] * mpk.u[i] * mpk.u[o];
                Lin_y[1+o] += g[ijo] * mpk.u[i] * mpk.u[j];
                
                int ij = 1+para.n+ (2*para.n-i)*(i+1)/2-(para.n-j);
                int jo = 1+para.n+ (2*para.n-j)*(j+1)/2-(para.n-o);
                int io = 1+para.n+ (2*para.n-i)*(i+1)/2-(para.n-o);
                Lin_y[ij] += g[ijo] * (-1) * mpk.u[o];
                Lin_y[io] += g[ijo] * (-1) * mpk.u[j];
                Lin_y[jo] += g[ijo] * (-1) * mpk.u[i];
                
                ijo++;
            }
        }
    }
    
}

void Cubic_KeyGen(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &g, int &k, Vec<ZZ_pE> &SK)
{
    Vec<ZZ_pE> Lin_y;
    Cubic_Gen_Lin_y(para, mpk,  k, g, Lin_y);
    Lin_KeyGen(para.Lin_pp, mpk.Lmpk, Lin_y, SK);
}

void Cubic_Dec(Cubic_pp &para, Cubic_MPK &mpk, Vec<ZZ_pE> &g, int &k,  Vec<ZZ_pE> &SK,  Vec<ZZ_pE> &CT, ZZ_pX &ans)
{
    Vec<ZZ_pE> Lin_CT;
    Vec<ZZ_pE> Lin_y;
    Cubic_Gen_Lin_y(para, mpk,  k, g, Lin_y);
    
    Lin_CT.SetLength(para.Lin_pp.L+2);
    for(int i = para.n; i< para.n+para.Lin_pp.L+2; ++i)
    Lin_CT[i-para.n] = CT[i];
    
    Lin_Dec(para.Lin_pp, Lin_y, Lin_CT, SK, ans);
    
    Turn_p1(para);
    int ijo =0;
    for(int i=0;i<para.n;++i)
    {
        for(int j=i; j<para.n;++j)
        {
            for(int o=j;o<para.n;++o)
            {
                ans += rep(g[ijo]*CT[i]*CT[j]*CT[o]);
                ijo++;
            }
            
        }
    }
    
    ZZ_p :: init(ZZ(para.p0));
    ZZ_pX tt;
    tt = ans;
 //   cout<< "tt = "<<tt<<endl;
    for(int i=0; i< deg(tt)+1; ++ i )
    {
        while(rep(tt[i]) > para.p0)  tt[i] = 0 +  tt[i];
    }
    ans = tt;
    
}
/*
void Cubic_test()
{
       clock_t startTime,endTime;
    
    Cubic_pp para;
    Cubic_MPK mpk;
    int Cn = n;
    int CQ = Q;
    Cubic_SetUp(para, mpk, Cn, CQ);
//    Cubic_pp_Print(para);
 //   Cubic_MPK_Print(mpk);
    
    Turn_p0(para);
    Vec<ZZ_pE> x;
    x.SetLength(para.n);
    for(int i =0;i<para.n;++i) x[i] = random_ZZ_pE();
//    cout<< "x = " << x << endl;
    
    
    Vec<ZZ_pE> g;
    int oL = (para.n*(para.n-1)*(para.n-2))/(3*2)+ para.n*para.n;
    g.SetLength(oL);
    for(int i =0;i<oL;++i) g[i] = random_ZZ_pE();
    
    startTime = clock();
    Vec<ZZ_pE> CT;
//    Cubic_Enc(para, mpk, x, CT);
    endTime = clock();
    cout <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000<< "ms" << endl;
//    cout<< "CT = " << CT << endl;

    startTime = clock();
    int k = 2;
    Vec<ZZ_pE> SK;
 //   Cubic_KeyGen(para, mpk, g, k, SK);
    endTime = clock();
    cout  <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000<< "ms" << endl;
//    cout<< "SK = " << SK << endl;
    
    startTime = clock();
    ZZ_pX ans;
    Cubic_Dec(para, mpk, g, k, SK, CT, ans);
    endTime = clock();
    cout <<(double)(endTime - startTime) / CLOCKS_PER_SEC *1000<< "ms" << endl;
//    cout<<ans<<endl;
 
}
*/
