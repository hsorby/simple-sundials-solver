
#include "odesolver.h"


void Solver::setProperties(const Properties &pProperties)
{
    // Set a property's value, but only if it is a valid property

    mProperties = pProperties;
}

void OdeSolver::initialize(double /**pVoi**/, int pRatesStatesCount,
                            double *pStates, double *pRates, double *pVariables,
                           ComputeRatesFunction pComputeRates)
{
    // Initialise the ODE solver

    mRatesStatesCount = pRatesStatesCount;

    mRates = pRates;
    mStates = pStates;
    mVariables = pVariables;

    mComputeRates = pComputeRates;
}

void OdeSolver::reinitialize(double /**pVoi**/)
{
    // Nothing to do by default...
}

