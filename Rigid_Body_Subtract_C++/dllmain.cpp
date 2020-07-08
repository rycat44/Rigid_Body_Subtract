// dllmain.cpp : Defines the entry point for the DLL application.
#include "pch.h"
#include <math.h>
#include "atlsafe.h"

CComVariant __stdcall RB_Subtract(double *refXYZ, double* refDisp, double* allXYZ, double* allDisp, int* UBound, double *scale, int *linear)
{
    // Store user parameters in DLL locally for use...
    double localScale = *scale;
    int localLinear = *linear;

    // The reference Translations & Rotatations
    double refT1, refT2, refT3, refR1, refR2, refR3; 

    // Init Matrix Arrays
    double M[3];
    double PB_prime[3];
    int localUBound = *UBound;
    double db_prime[3];
    double db_2prime[3];

    // Store the reference node Displacements & scale it ( scale usually == 1 )
    refT1 = refDisp[0] * localScale;
    refT2 = refDisp[1] * localScale;
    refT3 = refDisp[2] * localScale;
    refR1 = refDisp[3] * localScale;
    refR2 = refDisp[4] * localScale;
    refR3 = refDisp[5] * localScale;

    // RR Matrix
    double RR[3][3];

    if (localLinear)
    {
        //***************** Linearized Rotation Matrix **************************************************    
        RR[0][0] = 1;           RR[0][1] = -refR3;      RR[0][2] = refR2;
        RR[1][0] = refR3;       RR[1][1] = 1;           RR[1][2] = -refR1;
        RR[2][0] = -refR2;      RR[2][1] = refR1;       RR[2][2] = 1;

    }
    else
    {
        //***************** Normal Rotation Matrix **************************************************
        RR[0][0] = cos(refR2)*cos(refR3);                                           RR[0][1] = -cos(refR2)*sin(refR3);                                          RR[0][2] = sin(refR2);
        RR[1][0] = sin(refR1)*sin(refR2)*cos(refR3)+cos(refR1)*sin(refR3);          RR[1][1] = -sin(refR1)*sin(refR2)*sin(refR3) + cos(refR1)*cos(refR3);       RR[1][2] = -sin(refR1)*cos(refR2);
        RR[2][0] = -cos(refR1)*sin(refR2)*cos(refR3)+sin(refR1)*sin(refR3);         RR[2][1] = cos(refR1)*sin(refR2)*sin(refR3) + sin(refR1)*cos(refR3);        RR[2][2] = cos(refR1)*cos(refR2);
    }

    localUBound++;
    localUBound = (localUBound / 6)-1;

    // To store values and ultimately return to caller...
    CComSafeArray<double> sa_outVecT1(localUBound + 1);
    CComSafeArray<double> sa_outVecT2(localUBound + 1);
    CComSafeArray<double> sa_outVecT3(localUBound + 1);
    CComSafeArray<double> sa_outVecR1(localUBound + 1);
    CComSafeArray<double> sa_outVecR2(localUBound + 1);
    CComSafeArray<double> sa_outVecR3(localUBound + 1);

    for (int i = 0; i <= localUBound; i++)
    {
        //************ M = PB-PA... Position of B (PB) relative to Position of A (PA)
        M[0] = allXYZ[3 * i] - refXYZ[0];
        M[1] = allXYZ[3 * i+1] - refXYZ[1];
        M[2] = allXYZ[3 * i+2] - refXYZ[2];

        //************ PB' = [R][M] = Position of B due to rotations at A * **************
        PB_prime[0] = (RR[0][0] * M[0]) + (RR[0][1] * M[1]) + (RR[0][2] * M[2]);
        PB_prime[1] = (RR[1][0] * M[0]) + (RR[1][1] * M[1]) + (RR[1][2] * M[2]);
        PB_prime[2] = (RR[2][0] * M[0]) + (RR[2][1] * M[1]) + (RR[2][2] * M[2]);

        //************ dB' = PB' - M = Displacements at B due to rotations at A ***************
        db_prime[0] = PB_prime[0] - M[0];
        db_prime[1] = PB_prime[1] - M[1];
        db_prime[2] = PB_prime[2] - M[2];

        //************ dB''  = dB' + dA = Displacements at B due to rotations and translations at A * **************
        db_2prime[0] = db_prime[0] + refT1;
        db_2prime[1] = db_prime[1] + refT2;
        db_2prime[2] = db_prime[2] + refT3;

        //************ dB* = dB - dB'' = Displacements at B relative to A = total Nastran displacement

        sa_outVecT1.SetAt(i, allDisp[6 * i] * localScale - db_2prime[0] );
        sa_outVecT2.SetAt(i, allDisp[6 * i + 1] * localScale - db_2prime[1]);
        sa_outVecT3.SetAt(i, allDisp[6 * i + 2] * localScale - db_2prime[2]);

        //************ Rotations at B relative to A = Rotations at B - Rotations at A
        sa_outVecR1.SetAt(i, allDisp[6 * i + 3] * localScale - refR1);
        sa_outVecR2.SetAt(i, allDisp[6 * i + 4] * localScale - refR2);
        sa_outVecR3.SetAt(i, allDisp[6 * i + 5] * localScale - refR3);

    }

    CComVariant outVecT1 = sa_outVecT1.Detach();
    CComVariant outVecT2 = sa_outVecT2.Detach();
    CComVariant outVecT3 = sa_outVecT3.Detach();
    CComVariant outVecR1 = sa_outVecR1.Detach();
    CComVariant outVecR2 = sa_outVecR2.Detach();
    CComVariant outVecR3 = sa_outVecR3.Detach();

    CComSafeArray<VARIANT> sa_Disp(6);
    sa_Disp.SetAt(0, outVecT1);
    sa_Disp.SetAt(1, outVecT2);
    sa_Disp.SetAt(2, outVecT3);
    sa_Disp.SetAt(3, outVecR1);
    sa_Disp.SetAt(4, outVecR2);
    sa_Disp.SetAt(5, outVecR3);

    CComVariant vOut = sa_Disp.Detach();

    return vOut;
}



