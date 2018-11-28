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
#include "custom_constitutive/hencky_bounding_surface_plane_strain_2D_law.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlasticPlaneStrain2DLaw::HenckyBoundingSurfacePlasticPlaneStrain2DLaw()
    : HenckyElasticPlasticPlaneStrain2DLaw()
{
    mpHardeningLaw   = MPMHardeningLaw::Pointer( new BoundingSurfaceHardeningLaw() );
    mpYieldCriterion = MPMYieldCriterion::Pointer( new BoundingSurfaceYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule    = MPMFlowRule::Pointer( new BoundingSurfacePlasticFlowRule(mpYieldCriterion) );
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlasticPlaneStrain2DLaw::HenckyBoundingSurfacePlasticPlaneStrain2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
{
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  MPMYieldCriterion::Pointer( new BoundingSurfaceYieldCriterion(mpHardeningLaw) );
    mpMPMFlowRule     =  pMPMFlowRule;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyBoundingSurfacePlasticPlaneStrain2DLaw::HenckyBoundingSurfacePlasticPlaneStrain2DLaw(const HenckyBoundingSurfacePlasticPlaneStrain2DLaw& rOther)
    : HenckyElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyBoundingSurfacePlasticPlaneStrain2DLaw::Clone() const
{
    HenckyBoundingSurfacePlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyBoundingSurfacePlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyBoundingSurfacePlasticPlaneStrain2DLaw::~HenckyBoundingSurfacePlasticPlaneStrain2DLaw()
{
}

//*********************************CHECK**********************************************
//************************************************************************************

int HenckyBoundingSurfacePlasticPlaneStrain2DLaw::Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    HenckyElasticPlasticPlaneStrain2DLaw::Check(rProperties, rGeometry, rCurrentProcessInfo);

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
