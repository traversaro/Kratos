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

#include "custom_elements/femdem3d_large_displacement_element.hpp"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: FemDem3DElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: FemDem3DElement(NewId, pGeometry, pProperties)
{
	//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(FemDem3DLargeDisplacementElement const &rOther)
	: FemDem3DElement(rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
	//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

FemDem3DLargeDisplacementElement &FemDem3DLargeDisplacementElement::operator=(FemDem3DLargeDisplacementElement const &rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

	FemDem3DElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer FemDem3DLargeDisplacementElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	//NEEDED TO CREATE AN ELEMENT
	return Element::Pointer(new FemDem3DLargeDisplacementElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer FemDem3DLargeDisplacementElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
	//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

	FemDem3DLargeDisplacementElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	return Element::Pointer(new FemDem3DLargeDisplacementElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::~FemDem3DLargeDisplacementElement()
{
}

void FemDem3DLargeDisplacementElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{

}

void FemDem3DLargeDisplacementElement::CalculateLocalSystem(
	MatrixType &rLeftHandSideMatrix,
	VectorType &rRightHandSideVector,
	ProcessInfo &rCurrentProcessInfo)
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0, detF, detJ;

    const SizeType mat_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);

    rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

   
    // Loop over Gauss Points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
		J = this->GetGeometry().Jacobian(J, point_number, mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, mThisIntegrationMethod);

        double IntegrationWeight = integration_points[point_number].Weight() * detJ0;
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
		Vector N = row(Ncontainer, point_number);

        Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);
		// Taking into account Volume Force into de RHS
		for (unsigned int i = 0; i < number_of_nodes; i++) {
			int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++) {
				rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
			}
		}

        
        //GeometryUtils::DeformationGradient(J, rThisKinematicVariables.InvJ0,
        //                                   rThisKinematicVariables.F);
        //CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F,
        //           rThisKinematicVariables.DN_DX);





    }
}

double FemDem3DLargeDisplacementElement::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
{
    //GeometryType& r_geom = GetGeometry();
    //GeometryUtils::JacobianOnInitialConfiguration(
    //    r_geom,
    //    r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    //double detJ0;
    //MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    //const Matrix& rDN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    //GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    //return detJ0;

	return 0;
}


void FemDem3DLargeDisplacementElement::FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo)
{
}






} // namespace Kratos