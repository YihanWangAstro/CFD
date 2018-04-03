////////////////////////////////////////////////////////////////////////////////
//Filename:CFD.cpp                                                            //
//Author:Yihan Wang                                                           //
//Created:4.1.2018                                                            //
//                                                                            //
//Description:                                                                //
//  rough implement of 1D CFD for PHY688. Rush code no comment :-)            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include "CFD.h"
using namespace std;
double GAMMA = 1.4;
int Ugly_riemann_(double *,double *, double *, double *,double *, double *, double *,
                  double *, double *, double *);

void fluid::AdvanceTo(double tLimit)
{
    double dt = 0;
    for(;curTime <= tLimit;)
    {
        dt = GetTimeStepsize();
        tmpUGrid = UGrid;
        GetAdvSrc(UGrid);
        UpdateConserVars(UGrid, UGrid, 0.5*dt);
        GetAdvSrc(UGrid);
        UpdateConserVars(UGrid, tmpUGrid, dt);
        curTime += dt;
    }
}

void fluid::GetAdvSrc(FluVarList& src)
{
    SetBoundaryCells(src);
    ToFluVars(src);
    GetInterfaces();
    SolveRiemann();
}

void fluid::UpdateConserVars(FluVarList& dest, FluVarList& src, double dt)
{
    for(int i = beginIndex ; i < endIndex ; ++i)
        dest[i] = src[i] + AdvectSrc[i]*dt;
}

double fluid::GetMaxVel() const
{
    double maxVel = 0;
    for(int i = beginIndex ; i < endIndex ; ++i)
        maxVel = max(maxVel, fabs(qGrid[i].u) + sqrt(GAMMA * qGrid[i].p / qGrid[i].rho) );
    return maxVel;
}

double fluid::GetTimeStepsize() const
{
    return CFL*spaceSize/GetMaxVel();
}

FluVar fluid::GetMinmod(FluVar a, FluVar b) const
{
    return FluVar(Minmod(a.rho, b.rho), Minmod(a.u, b.u), Minmod(a.p, b.p));
}

double fluid::Minmod(double a, double b) const
{
    if(a*b > 0)
        return fabs(a) < fabs(b) ? a : b;
    else
        return 0;
}

void fluid::SolveRiemann()
{
    for(int i = beginIndex ; i < endIndex + 1 ; ++i)
    {
        Ugly_riemann_(&GAMMA, &interfaceL[i].rho, &interfaceL[i].u, &interfaceL[i].p,
                              &interfaceR[i].rho, &interfaceR[i].u, &interfaceR[i].p,
                              &flux[i].den,       &flux[i].rhou,    &flux[i].rhoE);
    }
    for(int i = beginIndex ; i < endIndex ; ++i)
        AdvectSrc[i] = (flux[i] - flux[i + 1])/spaceSize;
}

void fluid::CreateMesh(int gridNum, double scale)
{
    this->spaceSize = scale;
    flux.resize(gridNum);
    UGrid.resize(gridNum);
    qGrid.resize(gridNum);
    tmpUGrid.resize(gridNum);
    AdvectSrc.resize(gridNum);
    interfaceL.resize(gridNum);
    interfaceR.resize(gridNum);
}

void fluid::LoadInitState(FluVarList& initCond, double dx)
{
    curTime    = 0;
    beginIndex = 2;
    endIndex   = beginIndex + initCond.size();
    CreateMesh(initCond.size() + 4, dx);
    
    for(int i = beginIndex ; i < endIndex; ++i)
        qGrid[i] = initCond[i - 2];
    
    ToConserVars(UGrid);
    leftGhostCond  = UGrid[beginIndex];
    rightGhostCond = UGrid[endIndex - 1];
}

void fluid::GetInterfaces()
{
    for(int i = beginIndex; i < endIndex + 1; ++i)
    {
        interfaceL[i] = qGrid[i - 1] + GetMinmod(qGrid[i] - qGrid[i - 1], qGrid[i - 1] - qGrid[i - 2]) * 0.5;
        interfaceR[i] = qGrid[i] - GetMinmod(qGrid[i + 1] - qGrid[i], qGrid[i]- qGrid[i - 1]) * 0.5;
    }
}
void fluid::Output(const char* fname)
{
    ToFluVars(UGrid);
    FILE* fstream = fopen(fname,"w");
    if(!fstream)
        printf("cannot open the output file!"), exit(0);
    for(int i = beginIndex ; i < endIndex ; ++i)
        fprintf(fstream," %lf %lf %lf %lf %lf\n", spaceSize*(i - beginIndex), qGrid[i].rho, qGrid[i].u, qGrid[i].p, qGrid[i].p/qGrid[i].rho/(GAMMA - 1));
    fclose(fstream);
}

void fluid::Show(FluVarList& U)
{
    for(int i = 0 ; i < U.size() ; ++i)
        fprintf(stdout,"%lf %lf %lf\n", U[i].rho, U[i].u, U[i].p);
}

void fluid::SetBoundaryCells(FluVarList& U)
{
    U[beginIndex - 1] = leftGhostCond;
    U[beginIndex - 2] = leftGhostCond;
    U[endIndex]       = rightGhostCond;
    U[endIndex + 1]   = rightGhostCond;
}

void fluid::ToFluVars(FluVarList& U)
{
    for(int i = 0 ; i < U.size() ; ++i)
    {
        qGrid[i].rho = U[i].den;
        qGrid[i].u   = U[i].rhou / U[i].den;
        qGrid[i].p   = (GAMMA - 1)*( U[i].rhoE - 0.5 * qGrid[i].rho * qGrid[i].u * qGrid[i].u);
    }
}

void fluid::ToConserVars(FluVarList& U)
{
    for(int i = 0 ; i < U.size() ; ++i)
    {
        U[i].den  = qGrid[i].rho;
        U[i].rhou = qGrid[i].u * qGrid[i].rho;
        U[i].rhoE = qGrid[i].p/(GAMMA - 1) + 0.5 * qGrid[i].rho * qGrid[i].u * qGrid[i].u;
    }
}
