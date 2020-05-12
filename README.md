# FE-for-cubic-polynomials
An implementation of FE scheme for cubic polynomials 
This is an implementation of functional encryption scheme, consisting of a linear function encryption, a functional encryption scheme for cubic polynomials, an AIK randomized encoding scheme and a functional encryption scheme for all circuits. 

This code is based on NTL library 11.4.3 version. 

Files
-------------------------------------------------------------------
pp.h -- the public parameters file.

LinFE.hpp LinFE.cpp -- The linear function encryption scheme with Lin_Setup, Lin_KeyGen, Lin_Enc and Lin_Dec algorithms. 

CubicFE.hpp Cubic.cpp -- The functional encryption scheme for cubic polynomials  with Cubic_Setup, Cubic_KeyGen, Cubic_Enc and Cubic_Dec algorithms. 

AIK.hpp AIK.cpp -- The AIK randomized encoding scheme with AIK_Setup, AIK_Rebuild algorithms.

NCFE.hpp NCFE.cpp -- The functional encryption scheme for cubic polynomials  with FE_Setup, FE_KeyGen, FE_Enc and FE_Dec algorithms

main.cpp -- Here is an example of the functional encryption scheme for circuit C(x_1,x_2,x_3) = x_1x_2 + x_1(1-x_2)x_3 + (1-x_2)(1-x_3).

keys.txt -- Here are some keys generated from the example in main.cpp.
