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
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_bounding_surface_3D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlastic3DLaw::HenckyBoundingSurfacePlastic3DLaw()
    : HenckyElasticPlastic3DLaw()
{
    mpHardeningLaw   = MPMHardeningLaw::Pointer( new BoundingSurfaceHardeningLaw() );
    mpYieldCriterion = MPMYieldCriterion::Pointer( new BoundingSurfaceYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new BoundingSurfacePlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlastic3DLaw::HenckyBoundingSurfacePlastic3DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  MPMYieldCriterion::Pointer( new BoundingSurfaceYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyBoundingSurfacePlastic3DLaw::HenckyBoundingSurfacePlastic3DLaw(const HenckyBoundingSurfacePlastic3DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyBoundingSurfacePlastic3DLaw::Clone() const
{
    HenckyBoundingSurfacePlastic3DLaw::Pointer p_clone(new HenckyBoundingSurfacePlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlastic3DLaw::~HenckyBoundingSurfacePlastic3DLaw()
{
}

//*********************************CHECK**********************************************
//************************************************************************************

int HenckyBoundingSurfacePlastic3DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    HenckyElasticPlastic3DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

    // TODO: To be decided
    // KRATOS_ERROR_IF(PRE_CONSOLIDATION_STRESS.Key() == 0 || rProperties[PRE_CONSOLIDATION_STRESS] >= 0.00) << "PRE_CONSOLIDATION_STRESS has Key zero or invalid value (Expected negative value) " << std::endl;
    // KRATOS_ERROR_IF(OVER_CONSOLIDATION_RATIO.Key() == 0 || rProperties[OVER_CONSOLIDATION_RATIO] <= 0.00) << "OVER_CONSOLIDATION_RATIO has Key zero invalid value " << std::endl;

    // KRATOS_ERROR_IF(SWELLING_SLOPE.Key() == 0 || rProperties[SWELLING_SLOPE] <= 0.00) << "SWELLING_SLOPE has Key zero or invalid value " << std::endl;
    // KRATOS_ERROR_IF(NORMAL_COMPRESSION_SLOPE.Key() == 0 || rProperties[NORMAL_COMPRESSION_SLOPE] <= 0.00) << "NORMAL_COMPRESSION_SLOPE has Key zero or invalid value " << std::endl;
    // KRATOS_ERROR_IF(CRITICAL_STATE_LINE.Key() == 0 || rProperties[CRITICAL_STATE_LINE] <= 0.00) << "CRITICAL_STATE_LINE has Key zero or invalid value " << std::endl;
    // KRATOS_ERROR_IF(INITIAL_SHEAR_MODULUS.Key() == 0 || rProperties[INITIAL_SHEAR_MODULUS] <= 0.00) << "INITIAL_SHEAR_MODULUS has Key zero or invalid value " << std::endl;

    // KRATOS_ERROR_IF(ALPHA_SHEAR.Key() == 0 ) << "ALPHA_SHEAR has Key zero " << std::endl;

    return 0;
}

} // Namespace Kratos
