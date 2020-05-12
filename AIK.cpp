//
//  AIK.cpp
//  FE_implementation
//
//  Created by Zheng Zhang on 2020/5/5.
//  Copyright © 2020 Zheng Zhang. All rights reserved.
//

#include "AIK.hpp"
#include <stdio.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "pp.h"
using namespace NTL;
using namespace std;

void AIK_pp_Print(AIK_PP &para)
{
    cout<<"AIK_n = "<<para.n<<endl;
    cout<<"AIK_l = "<<para.l<<endl;
    cout<<"AIK_num = "<<para.num<<endl;
    cout<<"AIK_p = "<<para.p<<endl;
    cout<<"AIK_omega = "<<para.omega<<endl;
}

void AIK_SetUp(AIK_PP &para)
{
    para.n = AIK_n;
    para.l = AIK_l;
    para.p = ZZ(Cubic_p0);
    para.num = para.l*(para.l-1)/2;
    para.omega = 1 + para.n + para.l-2 + (para.l-1)*(para.l-2)/2;
}

void Matrix_Print(int len, int *A)
{
    for(int i=0;i<len;++i)
       {
           for(int j=0;j<len;++j)
           {
               cout<<A[i*len+j]<<" ";
           }
           cout<<endl;
       }
}

void PutinM( int &gi, int &g, int &gk)
{
    int a[3]={gi,g,gk};
    
    int temp;
    for(int i=0;i<3;++i)
    {
        temp = i;
        for(int j = i+1;j<3;++j)
        {
            if(a[temp]>a[j]) temp = j;
        }
        int tt = a[temp];
        a[temp] = a[i];
        a[i] = tt;
    }
    
    gi = a[0];
    g = a[1];
    gk = a[2];
    
}

void AIK_Rebuild(AIK_PP &para, int *Ax, Vec<ZZ_p> *F)
{
    int func_len;
    func_len = (para.omega)*(para.omega-1)*(para.omega-2)/(3*2)+para.omega*para.omega;
//    cout<<"func_len = "<<func_len<<endl;
    for(int i=0;i<para.num;++i) F[i].SetLength(func_len);
//    cout<<"func_len = "<<func_len<<endl;
    
    int Lx[(para.l-1)*(para.l-1)];
    for(int i=0;i<para.l-1;++i)
    {
        for(int j=0;j<para.l-1;++j)
        {
            Lx[i*(para.l-1)+j] = Ax[i*para.l+j+1];
        }
        if(i!=0) Lx[i*(para.l-1)+i-1] =0-para.n-1;
    }
    
    int R1[(para.l-1)*(para.l-1)];
    int R2[(para.l-1)*(para.l-1)];
    
    for(int i=0;i<(para.l-1)*(para.l-1);++i) R1[i] = R2[i] =0;
    for(int i=0;i<para.l-1;++i) R1[i*(para.l-1)+i] = R2[i*(para.l-1)+i] = 1;
   
    int temp =para.n+1;
    for(int i=0;i<para.l-1;++i)
    {
        if(i!= para.l-2) R2[(i+1)*(para.l-1)-1] = 9+i+1;
        for(int j=i+1;j<para.l-1;++j)
        {
            R1[i*(para.l-1)+j] = temp;
            temp++;
        }
    }
    for(int i =0;i<para.l-2;++i)
    {
        R2[(i+1)*(para.l-1)-1] = temp;
        temp++;
    }
/*
    cout<<"R1 = "<<endl;
    Matrix_Print(para.l-1, R1);
    
    cout<<"Lx = "<<endl;
    Matrix_Print(para.l-1, Lx);
    
    cout<<"R2 = "<<endl;
    Matrix_Print(para.l-1, R2);
*/
    temp =0;
    
    for(int i=0;i<para.l-1;++i)
    {
        for(int j=i;j<para.l-1;++j)
        {

            int matri_T[para.omega][para.omega][para.omega];
            memset(matri_T, 0, sizeof(matri_T));
            
            for(int ki = 0;ki<para.l-1;++ki)
            {
                for(int kj =0; kj<para.l-1;++kj)
                {
                    int fr1,fLx,fr2,gi=0,gk=0;
                    fr2 = kj*(para.l-1)+j;
                    fr1  = i*(para.l-1)+ki;
                    fLx = ki*(para.l-1)+kj;
                    
              //      cout<<"fr1 = "<< fr1 <<" fLx = "<< fLx << " fr2 = "<< fr2 <<endl;
                    
                    bool is =1;
                    if(R1[fr1] == 0){is=0;}
                    else if (R1[fr1] == 1){
                        gi = 0;
                    }
                    else{
                        gi = R1[fr1];
                    }
                    
                    
                    if(R2[fr2] == 0){is=0;}
                    else if(R2[fr2] == 1){
                        gk = 0;
                    }
                    else{
                        gk = R2[fr2];
                    }
                    
            //        cout<< "gi = "<< gi << " gk = "<<gk<<endl;
                    int g;
                    if(is){
                        if(Lx[fLx] == 0 ){ }
                        else if(Lx[fLx] == (0-para.n-1)){
                            g = 0;
                            PutinM( gi, g, gk);
                            matri_T[gi][g][gk] += -1;
                        //    cout<<"gi = "<<gi << " g = "<<g <<" gk = "<<gk<<" b = -1"<<endl;
                        }
                        else if(Lx[fLx] < 0 ){
                         // matri_T[gi][0][gk] += 1;
                            g = 0;
                            PutinM(gi,g,gk);
                            matri_T[gi][g][gk] += 1;
                          //  cout<<"gi = "<<gi << " g = "<<g <<" gk = "<<gk<<" b = 1"<<endl;
                         // matri_T[gi][0-Lx[fLx]][gk] += -1;
                            g = 0 - Lx[fLx];
                            PutinM( gi, g, gk);
                            matri_T[gi][g][gk] += -1;
                       //     cout<<"gi = "<<gi << " g = "<<g <<" gk = "<<gk<<" b = -1"<<endl;
                        }
                        else{
                        //    matri_T[gi][Lx[fLx]][gk] += 1;
                            g = Lx[fLx];
                            PutinM( gi, g, gk);
                            matri_T[gi][g][gk] += 1;
                         //   cout<<"gi = "<<gi << " g = "<<g <<" gk = "<<gk<<" b = 1"<<endl;
                        }
                    }
          //       cout<<"--------------------------------------------------"<<endl;
                }
            }
            
            ZZ_p::init(ZZ(Cubic_p0));
            int tt=0;
            for(int ki =0;ki<para.omega;++ki)
            {
                for(int kj=ki;kj<para.omega;++kj)
                {
                    for(int kk = kj;kk<para.omega;++kk)
                    {
                        F[temp][tt] = matri_T[ki][kj][kk];
                        tt++;
                    }
                }
            }
            temp++;
            
        }
    }
    
}
