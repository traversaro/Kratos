//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez

// System includes

// External includes

// Project includes
#include "custom_constitutive/fem_dem_elastic_law.hpp"

// #include "solid_mechanics_application_variables.h"

namespace Kratos
{
    //******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDemElasticLaw::FemDemElasticLaw()
    : LinearElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDemElasticLaw::FemDemElasticLaw(const FemDemElasticLaw& rOther)
    : LinearElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer FemDemElasticLaw::Clone() const
{
    return Kratos::make_shared<FemDemElasticLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDemElasticLaw::~FemDemElasticLaw()
{
}






} // Namespace Kratos