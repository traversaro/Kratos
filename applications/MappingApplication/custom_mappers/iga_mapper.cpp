//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "iga_mapper.h"
#include "custom_utilities/mapper_typedefs.h"

namespace Kratos
{

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void IGAInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                           const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();
    const SizeType num_nodes = p_geom->PointsNumber();

    const auto geom_family = p_geom->GetGeometryFamily();

    // TODO save Geometry
}

void IGALocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                                  EquationIdVectorType& rOriginIds,
                                  EquationIdVectorType& rDestinationIds,
                                  MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    KRATOS_DEBUG_ERROR_IF(mInterfaceInfos.size() != 1) << "IGALocalSystem currently only "
        << "works with on IGAInterfaceInfo!" << std::endl;
    }

    // Get the geometry (maybe change to pointer and swap them ...)
    the_geom;
    mInterfaceInfos[0].GetValue(the_geom);

    mIgaMapperCondition.GetGeometry() = the_geom;

    MatrixType integration_points_etc;
    mIgaMapperCondition.Calculate()


Calculate()
{
    mpModeler.ComputeIntegration(GetGeometry(), integration_points_etc, string: IntegrationMethod);
}


}

std::string IGASystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    buffer << " in rank " << CommRank;
    return buffer.str();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class IGAMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
