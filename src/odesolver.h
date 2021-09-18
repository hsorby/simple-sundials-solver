
#pragma once

#include <any>
#include <map>
#include <string>

#include "application_config.h"
#include "model_header.h"

using Properties = std::map<std::string, std::any>;

class Solver
{
public:
    virtual ~Solver() = default;
    void setProperties(const Properties &pProperties);

protected:
    Properties mProperties;
};


class OdeSolver : public Solver
{
public:
#if EXTERNAL_VARIABLES
    using ComputeRatesFunction = void (*)(double pVoi, double *pStates, double *pRates, double *pVariables, ExternalVariable externalVariable);
#else
    using ComputeRatesFunction = void (*)(double pVoi, double *pStates, double *pRates, double *pVariables);
#endif

    virtual void initialize(double pVoi, int pRatesStatesCount,
                            double *pStates, double *pRates, double *pVariables,
                            ComputeRatesFunction pComputeRates);
    virtual void reinitialize(double pVoi);

    virtual void solve(double &pVoi, double pVoiEnd) const = 0;

protected:
    int mRatesStatesCount = 0;

    double *mStates = nullptr;
    double *mRates = nullptr;
    double *mVariables = nullptr;

    ComputeRatesFunction mComputeRates = nullptr;
};
