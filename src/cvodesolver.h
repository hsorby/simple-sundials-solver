
#pragma once

#include "odesolver.h"

#include "nvector/nvector_serial.h"
#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nonlinearsolver.h"


static const auto MaximumStepId          = "MaximumStep";
static const auto MaximumNumberOfStepsId = "MaximumNumberOfSteps";
static const auto IntegrationMethodId    = "IntegrationMethod";
static const auto IterationTypeId        = "IterationType";
static const auto LinearSolverId         = "LinearSolver";
static const auto PreconditionerId       = "Preconditioner";
static const auto UpperHalfBandwidthId   = "UpperHalfBandwidth";
static const auto LowerHalfBandwidthId   = "LowerHalfBandwidth";
static const auto RelativeToleranceId    = "RelativeTolerance";
static const auto AbsoluteToleranceId    = "AbsoluteTolerance";
static const auto InterpolateSolutionId  = "InterpolateSolution";

//==============================================================================

static const auto AdamsMoultonMethod = "Adams-Moulton";
static const auto BdfMethod          = "BDF";

//==============================================================================

static const auto FunctionalIteration = "Functional";
static const auto NewtonIteration     = "Newton";

//==============================================================================

static const auto DenseLinearSolver    = "Dense";
static const auto BandedLinearSolver   = "Banded";
static const auto DiagonalLinearSolver = "Diagonal";
static const auto GmresLinearSolver    = "GMRES";
static const auto BiCgStabLinearSolver = "BiCGStab";
static const auto TfqmrLinearSolver    = "TFQMR";

//==============================================================================

static const auto NoPreconditioner     = "None";
static const auto BandedPreconditioner = "Banded";

//==============================================================================

// Default CVODES parameter values
// Note #1: a maximum step of 0 means that there is no maximum step as such and
//          that CVODES can use whatever step it sees fit...
// Note #2: CVODES' default maximum number of steps is 500, which ought to be
//          big enough in most cases...

static const double MaximumStepDefaultValue = 0.0;

enum {
    MaximumNumberOfStepsDefaultValue = 500
};

static const auto IntegrationMethodDefaultValue = BdfMethod;
static const auto IterationTypeDefaultValue = NewtonIteration;
static const auto LinearSolverDefaultValue = DenseLinearSolver;
static const auto PreconditionerDefaultValue = BandedPreconditioner;

enum {
    UpperHalfBandwidthDefaultValue = 0,
    LowerHalfBandwidthDefaultValue = 0
};

static const double RelativeToleranceDefaultValue = 1.0e-7;
static const double AbsoluteToleranceDefaultValue = 1.0e-7;

static const bool InterpolateSolutionDefaultValue = true;


class CvodeSolverUserData
{
public:
    explicit CvodeSolverUserData(double *pVariables, OdeSolver::ComputeRatesFunction pComputeRates);

    double *variables() const;

    OdeSolver::ComputeRatesFunction computeRates() const;

private:
    double *mVariables;

    OdeSolver::ComputeRatesFunction mComputeRates;
};

class CvodeSolver : public OdeSolver
{
public:
    ~CvodeSolver() override;

    void initialize(double pVoi, int pRatesStatesCount,
                    double *pStates, double *pRates, double *pVariables,
                    ComputeRatesFunction pComputeRates) override;
    void reinitialize(double pVoi) override;

    void solve(double &pVoi, double pVoiEnd) const override;

private:
    void *mSolver = nullptr;

    N_Vector mStatesVector = nullptr;

    SUNMatrix mMatrix = nullptr;
    SUNLinearSolver mLinearSolver = nullptr;
    SUNNonlinearSolver mNonLinearSolver = nullptr;

    CvodeSolverUserData *mUserData = nullptr;

    bool mInterpolateSolution = InterpolateSolutionDefaultValue;
};
