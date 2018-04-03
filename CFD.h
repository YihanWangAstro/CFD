////////////////////////////////////////////////////////////////////////////////
//Filename:CFD.h                                                              //
//Author:Yihan Wang                                                           //
//Created:4.1.2018                                                            //
//                                                                            //
//Description:                                                                //
//  rough implement of 1D CFD for PHY688                                      //
////////////////////////////////////////////////////////////////////////////////
#ifndef CFD_H
#define CFD_H
#include <vector>

struct FluVar
{
    union
    {
        struct{double rho, u, p;};
        struct{double den, rhou, rhoE;};
    };
    FluVar() : rho(0), u(0), p(0){}
    FluVar(double a, double b, double c) : rho(a), u(b), p(c){}
    FluVar& operator=(const FluVar& q)
    { rho = q.rho, u = q.u, p = q.p; return *this; }
    FluVar operator+(const FluVar& q) const
    { return FluVar(rho + q.rho, u + q.u, p + q.p); }
    FluVar operator-(const FluVar& q) const
    { return FluVar(rho - q.rho, u - q.u, p - q.p); }
    FluVar operator*(double c) const
    { return FluVar(rho * c, u * c, p * c); }
    FluVar operator/(double c) const
    { return FluVar(rho / c, u / c, p / c); }
};

typedef std::vector<FluVar> FluVarList;

class fluid
{
public:
    ////////////////////////////////Interface///////////////////////////////////
         fluid() : curTime(0), CFL(1){};
    void AdvanceTo(double tLimit);
    void LoadInitState(FluVarList& initCond, double dx);
    void Output(const char* fname);
    void Show(FluVarList& U);
    /////////////////////////////Member variables///////////////////////////////
public:
    double curTime;
    double CFL;
    //////////////////////////////Private Member////////////////////////////////
private:
    double     spaceSize;
    int        beginIndex;
    int        endIndex;
    FluVarList UGrid;//conserVars grid
    FluVarList qGrid;//FluVars grid
    FluVarList AdvectSrc;//advective source grid
    FluVarList interfaceL;
    FluVarList interfaceR;
    FluVarList tmpUGrid;//tmp conserVars grid
    FluVarList flux;
    FluVar     leftGhostCond;
    FluVar     rightGhostCond;
    ////////////////////////////Private Function////////////////////////////////
    void   CreateMesh(int gridNum, double scale);
    void   GetAdvSrc(FluVarList& src);
    void   GetInterfaces();
    FluVar GetMinmod(FluVar a, FluVar b) const;
    double GetMaxVel() const;
    double GetTimeStepsize() const;
    double Minmod(double a, double b) const;
    void   SetBoundaryCells(FluVarList& U);
    void   SolveRiemann();
    void   ToFluVars(FluVarList& U);
    void   ToConserVars(FluVarList& U);
    void   UpdateConserVars(FluVarList& dest, FluVarList& src, double dt);
};

#endif
