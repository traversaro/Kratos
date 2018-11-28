//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"
#include "custom_constitutive/yield_criteria/bounding_surface_yield_criterion.hpp"
#include "includes/mat_variables.h"
#include "particle_mechanics_application_variables.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************
BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion()
    :MPMYieldCriterion()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion(HardeningLawPointer pHardeningLaw)
    :MPMYieldCriterion(pHardeningLaw)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

BoundingSurfaceYieldCriterion& BoundingSurfaceYieldCriterion::operator=(BoundingSurfaceYieldCriterion const& rOther)
{
    MPMYieldCriterion::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion(BoundingSurfaceYieldCriterion const& rOther)
    :MPMYieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

BoundingSurfaceYieldCriterion::~BoundingSurfaceYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************
double& BoundingSurfaceYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rPreconsolidationPressure)
{
    // Compute three invariants
    double mean_stress_p, deviatoric_q, lode_angle;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rStressVector, mean_stress_p, deviatoric_q, lode_angle);
    mean_stress_p *= -1.0;  // Using negative definition

    // Compute M_cs
    const bool fix_csl_M = this->GetHardeningLaw().GetProperties()[IS_CSL_FIX];
    double shear_M       = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
    if (!fix_csl_M)
        shear_M = this->CalculateCriticalStateLineSlope(lode_angle);

    // Get constants
    const double curvature_N = this->GetHardeningLaw().GetProperties()[BOUNDING_SURFACE_CURVATURE];
    const double ratio_R     = this->GetHardeningLaw().GetProperties()[MODEL_PARAMETER_R];

    // f = (Q/MP)^N - ln(P_c/P)/ln(R)
    rStateFunction  = std::pow(deviatoric_q/(shear_M*mean_stress_p), curvature_N);
    rStateFunction -= std::log(rPreconsolidationPressure/mean_stress_p) / std::log(ratio_R);

    return rStateFunction;
}


//*******************************CALCULATE FIRST YIELD FUNCTION DERIVATIVE *****************
//************************************************************************************
void BoundingSurfaceYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative)
{
    // Compute three invariants
    double mean_stress_p, deviatoric_q, lode_angle;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rStressVector, mean_stress_p, deviatoric_q, lode_angle);
    mean_stress_p *= -1.0;  // Using negative definition

    // Compute M_cs
    const bool fix_csl_M = this->GetHardeningLaw().GetProperties()[IS_CSL_FIX];
    double shear_M       = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
    if (!fix_csl_M)
        shear_M = this->CalculateCriticalStateLineSlope(lode_angle);

    // Get constants
    const double curvature_N = this->GetHardeningLaw().GetProperties()[BOUNDING_SURFACE_CURVATURE];
    const double ratio_R     = this->GetHardeningLaw().GetProperties()[MODEL_PARAMETER_R];
    const double alpha       = this->GetAlphaParameter();

    const double aux_multiplier = std::pow( deviatoric_q / (shear_M*mean_stress_p), curvature_N) ;

    rFirstDerivative.resize(3, false);
    // (df/dP)
    rFirstDerivative[0]  = -curvature_N / mean_stress_p * aux_multiplier + 1.0 / (mean_stress_p * std::log(ratio_R));
    // (df/dQ)
    rFirstDerivative[1]  =  curvature_N / deviatoric_q  * aux_multiplier;
    // (df/dØ)
    rFirstDerivative[2]  = -3.0/4.0 * curvature_N * aux_multiplier;
    rFirstDerivative[2] *= (1.0 - std::pow(alpha, 4)) * std::cos(3.0 * lode_angle) / (1.0 + std::pow(alpha, 4) - (1 - std::pow(alpha, 4)) * std::sin(3.0 * lode_angle) );

}


double BoundingSurfaceYieldCriterion::CalculateCriticalStateLineSlope(const double& rLodeAngle)
{
    double shear_M = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
    const double alpha = this->GetAlphaParameter();

    const double aux_multiplier = 2.0 * std::pow(alpha, 4) / (1.0 + std::pow(alpha, 4) - (1 - std::pow(alpha, 4)) * std::sin(3.0 * rLodeAngle) );
    shear_M *= std::pow(aux_multiplier, 0.25);

    KRATOS_ERROR_IF(shear_M < 0.0) << "The slope of critical state line is negative! M_cs = " << shear_M << std::endl;

    return shear_M;
}


double BoundingSurfaceYieldCriterion::GetAlphaParameter()
{
    const double shear_M = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
    const double phi_csl = (3.0 * shear_M) / (6.0 + shear_M);
    const double alpha   = (3.0 - std::sin(phi_csl)) / (3.0 + std::sin(phi_csl));

    return alpha;
}


double BoundingSurfaceYieldCriterion::GetPI()
{
    return std::atan(1.0)*4.0;
}

void BoundingSurfaceYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMYieldCriterion )
}

void BoundingSurfaceYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMYieldCriterion )
}


}  // namespace Kratos.
