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

#if !defined(KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED

#include "custom_elements/femdem3d_element.hpp"

namespace Kratos
{
class FemDem3DHexahedronElement : public FemDem3DElement 
{

  public:
	/// Default constructors
	FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry);

	FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	FemDem3DHexahedronElement(FemDem3DHexahedronElement const &rOther);

	/// Destructor.
	virtual ~FemDem3DHexahedronElement();

	/// Assignment operator.
	FemDem3DHexahedronElement &operator=(FemDem3DHexahedronElement const &rOther);

	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	FemDem3DHexahedronElement()
	{
	}

	void CalculateLocalSystem(
		MatrixType &rLeftHandSideMatrix,
		VectorType &rRightHandSideVector,
		ProcessInfo &rCurrentProcessInfo);

	void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo);

	void GetValueOnIntegrationPoints(
		const Variable<Matrix> &rVariable,
		std::vector<Matrix> &rValues,
		const ProcessInfo &rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rOutput,
		const ProcessInfo &rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rOutput,
		const ProcessInfo &rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(
		const Variable<Matrix> &rVariable,
		std::vector<Matrix> &rOutput,
		const ProcessInfo &rCurrentProcessInfo);


  private:

	int mNumberOfEdges = 12;


}; // Class FemDem3DHexahedronElement

} // Namespace Kratos
#endif // KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED  defined