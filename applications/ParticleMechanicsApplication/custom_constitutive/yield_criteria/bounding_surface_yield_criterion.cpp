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


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************
BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion()
    :YieldCriterion()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion(HardeningLawPointer pHardeningLaw)
    :YieldCriterion(pHardeningLaw)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

BoundingSurfaceYieldCriterion& BoundingSurfaceYieldCriterion::operator=(BoundingSurfaceYieldCriterion const& rOther)
{
    YieldCriterion::operator=(rOther);
    return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

BoundingSurfaceYieldCriterion::BoundingSurfaceYieldCriterion(BoundingSurfaceYieldCriterion const& rOther)
    :YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

BoundingSurfaceYieldCriterion::~BoundingSurfaceYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& BoundingSurfaceYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha, const double& rOldPreconsolidationPressure)
{
    double mean_stress_p, deviatoric_q;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rStressVector, mean_stress_p, deviatoric_q);
    deviatoric_q *= std::sqrt(3.0); //Q = sqrt(3) * J2

    const double shear_M = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

    double preconsolidation_stress = 0.0;
    preconsolidation_stress = mpHardeningLaw->CalculateHardening(preconsolidation_stress, rAlpha, rOldPreconsolidationPressure);
    
    // f = (Q/M)² + P (P - P_c)
    rStateFunction = std::pow(deviatoric_q/shear_M, 2);
    rStateFunction += (mean_stress_p * (mean_stress_p - preconsolidation_stress) );

    return rStateFunction;
}


//*******************************CALCULATE FIRST YIELD FUNCTION DERIVATIVE *****************
//************************************************************************************
void BoundingSurfaceYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative, const double& rAlpha, const double& rOldPreconsolidationPressure)
{
    double mean_stress_p, deviatoric_q;

    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants( rStressVector, mean_stress_p, deviatoric_q);
    deviatoric_q *= std::sqrt(3.0); //Q = sqrt(3) * J2

    const double shear_M = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

    double preconsolidation_stress = 0.0;
    preconsolidation_stress = mpHardeningLaw->CalculateHardening(preconsolidation_stress, rAlpha, rOldPreconsolidationPressure);

    rFirstDerivative.resize(3, false);
    rFirstDerivative[0] = 2.0 * mean_stress_p - preconsolidation_stress; // (df/dP)
    rFirstDerivative[1] = 2.0 * deviatoric_q / std::pow(shear_M, 2);         // (df/dQ)
    rFirstDerivative[2] = - mean_stress_p;                              // (df/dP_c)
}

//*******************************CALCULATE SECOND YIELD FUNCTION DERIVATIVE *****************
//************************************************************************************
void BoundingSurfaceYieldCriterion::CalculateYieldFunctionSecondDerivative(const Vector& rStressVector, Vector& rSecondDerivative)
{
    const double shear_M = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

    rSecondDerivative.resize(6, false);
    rSecondDerivative[0] = 2.0 ;                        // (df²/dP²)  
    rSecondDerivative[1] = 2.0 / std::pow(shear_M, 2) ; // (df²/dQ²)  
    rSecondDerivative[2] = 0.0 ;                        // (df²/dP_c²)
    rSecondDerivative[3] = 0.0 ;                        // (df²/dPdQ)  
    rSecondDerivative[4] = 0.0 ;                        // (df²/dQdP_c)  
    rSecondDerivative[5] =-1.0 ;                        // (df²/dPdP_c)

}

double BoundingSurfaceYieldCriterion::GetPI()
{
    return std::atan(1.0)*4.0;
}

void BoundingSurfaceYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void BoundingSurfaceYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
