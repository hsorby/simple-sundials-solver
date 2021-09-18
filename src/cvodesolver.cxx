
#include "cvodesolver.h"

#include <iostream>

#include "cvodes/cvodes.h"
#include "cvodes/cvodes_bandpre.h"
#include "cvodes/cvodes_diag.h"
#include "cvodes/cvodes_direct.h"
#include "cvodes/cvodes_spils.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spbcgs.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunlinsol/sunlinsol_sptfqmr.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#include "application_config.h"
#if EXTERNAL_VARIABLES
extern "C"
{
#include "external_variables.h"
}
#endif

int rhsFunction(double pVoi, N_Vector pStates, N_Vector pRates, void *pUserData)
{
    // Compute the RHS function

    auto userData = static_cast<CvodeSolverUserData *>(pUserData);

    userData->computeRates()(
                pVoi,
                N_VGetArrayPointer_Serial(pStates),
                N_VGetArrayPointer_Serial(pRates),
                userData->variables()
                #if EXTERNAL_VARIABLES
                , computeExternalVariable
                #endif
                );

    return 0;
}


void errorHandler(int pErrorCode, const char */**pModule**/, const char */**pFunction**/,
                  char *pErrorMessage, void */**pUserData**/)
{
    // Forward errors to our CvodeSolver object

    if (pErrorCode != CV_WARNING) {
        std::cout << "Cvode Error: " << pErrorMessage << std::endl;
    }
}

CvodeSolverUserData::CvodeSolverUserData(double *pVariables,
                                         OdeSolver::ComputeRatesFunction pComputeRates) :
    mVariables(pVariables),
    mComputeRates(pComputeRates)
{
}

double *CvodeSolverUserData::variables() const
{
    // Return our constants array

    return mVariables;
}

OdeSolver::ComputeRatesFunction CvodeSolverUserData::computeRates() const
{
    // Return our compute rates function

    return mComputeRates;
}

CvodeSolver::~CvodeSolver()
{
    // Make sure that the solver has been initialised

    if (mSolver == nullptr) {
        return;
    }

    // Delete some internal objects

    N_VDestroy_Serial(mStatesVector);
    SUNLinSolFree(mLinearSolver);

    SUNNonlinSolFree(mNonLinearSolver);
    SUNMatDestroy(mMatrix);

    CVodeFree(&mSolver);

    delete mUserData;
}

void CvodeSolver::initialize(double pVoi, int pRatesStatesCount,
                             double *pStates, double *pRates, double *pVariables,
                             ComputeRatesFunction pComputeRates)
{
    // Retrieve our properties

    double maximumStep = MaximumStepDefaultValue;
    int maximumNumberOfSteps = MaximumNumberOfStepsDefaultValue;
    std::string integrationMethod = IntegrationMethodDefaultValue;
    std::string iterationType = IterationTypeDefaultValue;
    std::string linearSolver = LinearSolverDefaultValue;
    std::string preconditioner = PreconditionerDefaultValue;
    int upperHalfBandwidth = UpperHalfBandwidthDefaultValue;
    int lowerHalfBandwidth = LowerHalfBandwidthDefaultValue;
    double relativeTolerance = RelativeToleranceDefaultValue;
    double absoluteTolerance = AbsoluteToleranceDefaultValue;

    if (mProperties.count(MaximumStepId) > 0) {
        maximumStep = std::any_cast<double>(mProperties.at(MaximumStepId));
    } else {
        std::cout << "Using default value for 'maximum step' property: " << maximumStep << std::endl;
    }

    if (mProperties.count(MaximumNumberOfStepsId)) {
        maximumNumberOfSteps = std::any_cast<int>(mProperties.at(MaximumNumberOfStepsId));
    } else {
        std::cout << "Using default value for 'Maximum number of steps' property: " << maximumNumberOfSteps << std::endl;
    }

    if (mProperties.count(IntegrationMethodId)) {
        integrationMethod = std::any_cast<std::string>(mProperties.at(IntegrationMethodId));
    } else {
        std::cout << "Using default value for 'Integration method' property: " << integrationMethod << std::endl;
    }

    if (mProperties.count(IterationTypeId)) {
        iterationType = std::any_cast<std::string>(mProperties.at(IterationTypeId));

        if (iterationType == NewtonIteration) {
            // We are dealing with a Newton iteration, so retrieve and check its
            // linear solver

            if (mProperties.count(LinearSolverId)) {
                linearSolver = std::any_cast<std::string>(mProperties.at(LinearSolverId));

                bool needUpperAndLowerHalfBandwidths = false;

                if (   (linearSolver == DenseLinearSolver)
                       || (linearSolver == DiagonalLinearSolver)) {
                    // We are dealing with a dense/diagonal linear solver, so
                    // nothing more to do
                } else if (linearSolver == BandedLinearSolver) {
                    // We are dealing with a banded linear solver, so we need
                    // both an upper and a lower half bandwidth

                    needUpperAndLowerHalfBandwidths = true;
                } else {
                    // We are dealing with a GMRES/Bi-CGStab/TFQMR linear
                    // solver, so retrieve and check its preconditioner

                    if (mProperties.count(PreconditionerId)) {
                        preconditioner = std::any_cast<std::string>(mProperties.at(PreconditionerId));
                    } else {
                        std::cout << "Using default value for 'Preconditioner' property: " << preconditioner << std::endl;
                    }

                    if (preconditioner == BandedPreconditioner) {
                        // We are dealing with a banded preconditioner, so we
                        // need both an upper and a lower half bandwidth

                        needUpperAndLowerHalfBandwidths = true;
                    }
                }

                if (needUpperAndLowerHalfBandwidths) {
                    if (mProperties.count(UpperHalfBandwidthId)) {
                        upperHalfBandwidth = std::any_cast<int>(mProperties.at(UpperHalfBandwidthId));

                        if (upperHalfBandwidth >= pRatesStatesCount) {
                            std::cout << "Error: the 'Lower half-bandwidth' property value is not valid: " << upperHalfBandwidth << ", for states count: " << pRatesStatesCount << std::endl;
                            return;
                        }
                    } else {
                        std::cout << "Using default value for 'Upper half-bandwidth' property: " << upperHalfBandwidth << std::endl;
                    }

                    if (mProperties.count(LowerHalfBandwidthId)) {
                        lowerHalfBandwidth = std::any_cast<int>(mProperties.at(LowerHalfBandwidthId));

                        if (lowerHalfBandwidth >= pRatesStatesCount) {
                            std::cout << "Error: the 'Lower half-bandwidth' property value is not valid: " << lowerHalfBandwidth << ", for states count: " << pRatesStatesCount << std::endl;
                            return;
                        }
                    } else {
                        std::cout << "Using default value for 'Lower half-bandwidth' property: " << lowerHalfBandwidth << std::endl;
                    }
                }
            } else {
                std::cout << "Using default value for 'Linear solver' property: " << linearSolver << std::endl;
            }
        }
    } else {
        std::cout << "Using default value for 'Iteration type' property: " << iterationType << std::endl;
    }

    if (mProperties.count(RelativeToleranceId)) {
        relativeTolerance = std::any_cast<double>(mProperties.at(RelativeToleranceId));
    } else {
        std::cout << "Using default value for 'Relative tolerance' property: " << relativeTolerance << std::endl;
    }

    if (mProperties.count(AbsoluteToleranceId)) {
        absoluteTolerance = std::any_cast<double>(mProperties.at(AbsoluteToleranceId));
    } else {
        std::cout << "Using default value for 'Absolute tolerance' property: " << absoluteTolerance << std::endl;
    }

    if (mProperties.count(InterpolateSolutionId)) {
        mInterpolateSolution = std::any_cast<bool>(mProperties.at(InterpolateSolutionId));
    } else {
        std::cout << "Using default value for 'Interpolate solution' property: " << mInterpolateSolution << std::endl;
    }

    // Initialise our ODE solver

    OdeSolver::initialize(pVoi, pRatesStatesCount, pStates, pRates, pVariables, pComputeRates);

    // Create our states vector

    mStatesVector = N_VMake_Serial(pRatesStatesCount, pStates);

    // Create our CVODES solver

    bool newtonIteration = iterationType == NewtonIteration;

    mSolver = CVodeCreate((integrationMethod == BdfMethod) ? CV_BDF : CV_ADAMS);

    // Use our own error handler

    CVodeSetErrHandlerFn(mSolver, errorHandler, this);

    // Initialise our CVODES solver

    CVodeInit(mSolver, rhsFunction, pVoi, mStatesVector);

    // Set our user data

    mUserData = new CvodeSolverUserData(pVariables, pComputeRates);

    CVodeSetUserData(mSolver, mUserData);

    // Set our maximum step

    CVodeSetMaxStep(mSolver, maximumStep);

    // Set our maximum number of steps

    CVodeSetMaxNumSteps(mSolver, maximumNumberOfSteps);

    // Set our linear solver, if needed

    if (newtonIteration) {
        if (linearSolver == DenseLinearSolver) {
            mMatrix = SUNDenseMatrix(pRatesStatesCount, pRatesStatesCount);
            mLinearSolver = SUNLinSol_Dense(mStatesVector, mMatrix);

            CVodeSetLinearSolver(mSolver, mLinearSolver, mMatrix);
        } else if (linearSolver == BandedLinearSolver) {
            mMatrix = SUNBandMatrix(pRatesStatesCount, upperHalfBandwidth,
                                    lowerHalfBandwidth);
            mLinearSolver = SUNLinSol_Band(mStatesVector, mMatrix);

            CVodeSetLinearSolver(mSolver, mLinearSolver, mMatrix);
        } else if (linearSolver == DiagonalLinearSolver) {
            CVDiag(mSolver);
        } else {
            // We are dealing with a GMRES/Bi-CGStab/TFQMR linear solver

            if (preconditioner == BandedPreconditioner) {
                if (linearSolver == GmresLinearSolver) {
                    mLinearSolver = SUNLinSol_SPGMR(mStatesVector, PREC_LEFT, 0);
                } else if (linearSolver == BiCgStabLinearSolver) {
                    mLinearSolver = SUNLinSol_SPBCGS(mStatesVector, PREC_LEFT, 0);
                } else {
                    mLinearSolver = SUNLinSol_SPTFQMR(mStatesVector, PREC_LEFT, 0);
                }

                CVodeSetLinearSolver(mSolver, mLinearSolver, mMatrix);
                CVBandPrecInit(mSolver, pRatesStatesCount, upperHalfBandwidth,
                               lowerHalfBandwidth);
            } else {
                if (linearSolver == GmresLinearSolver) {
                    mLinearSolver = SUNLinSol_SPGMR(mStatesVector, PREC_NONE, 0);
                } else if (linearSolver == BiCgStabLinearSolver) {
                    mLinearSolver = SUNLinSol_SPBCGS(mStatesVector, PREC_NONE, 0);
                } else {
                    mLinearSolver = SUNLinSol_SPTFQMR(mStatesVector, PREC_NONE, 0);
                }

                CVodeSetLinearSolver(mSolver, mLinearSolver, mMatrix);
            }
        }
    } else {
        mNonLinearSolver = SUNNonlinSol_FixedPoint(mStatesVector, 0);

        CVodeSetNonlinearSolver(mSolver, mNonLinearSolver);
    }

    // Set our relative and absolute tolerances

    CVodeSStolerances(mSolver, relativeTolerance, absoluteTolerance);
}

void CvodeSolver::reinitialize(double pVoi)
{
    // Reinitialise our CVODES object

    CVodeReInit(mSolver, pVoi, mStatesVector);
}

void CvodeSolver::solve(double &pVoi, double pVoiEnd) const
{
    // Solve the model

    if (!mInterpolateSolution) {
        CVodeSetStopTime(mSolver, pVoiEnd);
    }

    CVode(mSolver, pVoiEnd, mStatesVector, &pVoi, CV_NORMAL);

    // Compute the rates one more time to get up to date values for the rates
    // Note: another way of doing this would be to copy the contents of the
    //       calculated rates in rhsFunction, but that's bound to be more time
    //       consuming since a call to CVode() is likely to generate at least a
    //       few calls to rhsFunction(), so that would be quite a few memory
    //       transfers while here we 'only' compute the rates one more time...

    mComputeRates(pVoiEnd, N_VGetArrayPointer_Serial(mStatesVector), mRates, mVariables
#if EXTERNAL_VARIABLES
                  , computeExternalVariable
#endif
    );
}
