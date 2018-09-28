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
//

#include "includes/define.h"
#include "femdem3d_hexahedron_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

FemDem3DHexahedronElement::FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: FemDem3DElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDem3DHexahedronElement::FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: FemDem3DElement(NewId, pGeometry, pProperties)
{
	//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDem3DHexahedronElement::FemDem3DHexahedronElement(FemDem3DHexahedronElement const &rOther)
	: FemDem3DElement(rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
	//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

FemDem3DHexahedronElement &FemDem3DHexahedronElement::operator=(FemDem3DHexahedronElement const &rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

	FemDem3DElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer FemDem3DHexahedronElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	//NEEDED TO CREATE AN ELEMENT
	return Element::Pointer(new FemDem3DHexahedronElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer FemDem3DHexahedronElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
	//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

	FemDem3DHexahedronElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	return Element::Pointer(new FemDem3DHexahedronElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDem3DHexahedronElement::~FemDem3DHexahedronElement()
{
}


void FemDem3DHexahedronElement::CalculateLocalSystem(
	MatrixType &rLeftHandSideMatrix,
	VectorType &rRightHandSideVector,
	ProcessInfo &rCurrentProcessInfo)
{

}






} // namespace Kratos