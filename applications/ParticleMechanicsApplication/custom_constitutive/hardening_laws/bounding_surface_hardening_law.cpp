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
#include "includes/properties.h"
#include "custom_constitutive/hardening_laws/bounding_surface_hardening_law.hpp"
#include "includes/mat_variables.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

BoundingSurfaceHardeningLaw::BoundingSurfaceHardeningLaw()
	:HardeningLaw()
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

BoundingSurfaceHardeningLaw& BoundingSurfaceHardeningLaw::operator=(BoundingSurfaceHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

BoundingSurfaceHardeningLaw::BoundingSurfaceHardeningLaw(BoundingSurfaceHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

BoundingSurfaceHardeningLaw::~BoundingSurfaceHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& BoundingSurfaceHardeningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const Variable<double>& rThisVariable)
{

    return rHardening;
}
  

void BoundingSurfaceHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )

}

void BoundingSurfaceHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )

}


}  // namespace Kratos.
