////////////////////////////////////////////////////////////////////////////////
//Filename:main.cpp                                                           //
//Author:Yihan Wang                                                           //
//Created:4.1.2018                                                            //
//                                                                            //
//Description:                                                                //
//  test of 1D CFD for PHY688                                                 //
////////////////////////////////////////////////////////////////////////////////
#include"CFD.h"
using namespace std;
void CreateInitCond(FluVarList& initCond, FluVar& qL, FluVar& qR, int gridNum)
{
    initCond.clear();
    int halfLen = gridNum >> 2;
    
    for(int i = 0 ; i < halfLen; i++)
        initCond.push_back(qL);
    
    for(int i = 0 ; i < halfLen; i++)
        initCond.push_back(qR);
}

void Test(const char* runName, FluVarList& initCond, double tLimit)
{
    int meshNumber = initCond.size();
    fluid hydro;
    hydro.LoadInitState(initCond, 1.0/meshNumber);
    hydro.AdvanceTo(tLimit);
    hydro.Output(runName);
}

int main()
{
    FluVarList initCond;
    FluVar qL, qR;
    
    {
        qL = FluVar(1.0,   0,   1.0);
        qR = FluVar(0.125, 0.0, 0.1);
        /*CreateInitCond(initCond, qL, qR, 64);
        Test("test1_64.dat",  initCond, 0.2);
        CreateInitCond(initCond, qL, qR,128);
        Test("test1_128.dat", initCond, 0.2);*/
        CreateInitCond(initCond, qL, qR,256);
        Test("test1_256.dat", initCond, 0.2);
    }
 
    {
        qL = FluVar(1.0, -2.0, 0.4);
        qR = FluVar(1.0,  2.0, 0.4);
        /*CreateInitCond(initCond, qL, qR,  64);
        Test("test2_64.dat",  initCond, 0.15);
        CreateInitCond(initCond, qL, qR, 128);
        Test("test2_128.dat", initCond, 0.15);*/
        CreateInitCond(initCond, qL, qR, 256);
        Test("test2_256.dat", initCond, 0.15);
    }
    
    {
        qL = FluVar(1.0, 0.0, 1000.0);
        qR = FluVar(1.0, 0.0, 0.01);
        /*CreateInitCond(initCond, qL, qR, 64);
        Test("test3_64.dat",  initCond, 0.012);
        CreateInitCond(initCond, qL, qR, 128);
        Test("test3_128.dat", initCond, 0.012);*/
        CreateInitCond(initCond, qL, qR, 256);
        Test("test3_256.dat", initCond, 0.012);
    }
    
    {
        qL = FluVar(5.6698, -1.4701, 100.0);
        qR = FluVar(1.0,    -10.5,   1.0);
        /*CreateInitCond(initCond, qL, qR, 64);
        Test("test4_64.dat",  initCond, 1.0);
        CreateInitCond(initCond, qL, qR,128);
        Test("test4_128.dat", initCond, 1.0);*/
        CreateInitCond(initCond, qL, qR,256);
        Test("test4_256.dat", initCond, 1.0);
    }
    return 0;
}

