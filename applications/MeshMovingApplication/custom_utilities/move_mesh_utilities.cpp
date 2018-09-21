//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//                   Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "move_mesh_utilities.h"

namespace Kratos {

MeshVelocityCalculationUtility::MeshVelocityCalculationUtility(ModelPart& rMeshMovingModelPart,
                                   Parameters MeshVolocityCalculationParameters)
{
    // default params ...

}

void MeshVelocityCalculationUtility::CalculateMeshVelocities()
{

}

namespace MoveMeshUtilities {

//******************************************************************************
//******************************************************************************
void CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, GeometryType &rGeometry) {
  KRATOS_TRY;

  const IntegrationMethod this_integration_method =
      rGeometry.GetDefaultIntegrationMethod();
  const GeometryType::IntegrationPointsArrayType &integration_points =
      rGeometry.IntegrationPoints(this_integration_method);

  if (rInvJ0.size() != integration_points.size())
    rInvJ0.resize(integration_points.size());
  if (rDetJ0.size() != integration_points.size())
    rDetJ0.resize(integration_points.size());

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************

void CalculateMeshVelocities(ModelPart* pMeshModelPart,
                             const int TimeOrder, const double DeltaTime) {

    CalculateMeshVelocities(*pMeshModelPart, TimeOrder, DeltaTime);
}

void CalculateMeshVelocities(ModelPart &rMeshModelPart,
                             const int TimeOrder, const double DeltaTime) {
  KRATOS_TRY;

  KRATOS_ERROR_IF(DeltaTime <= 0.0) << "Invalid DELTA_TIME." << std::endl;

  const double coeff = 1 / DeltaTime;

  if (TimeOrder == 1) {
    for (ModelPart::NodeIterator i =
             rMeshModelPart.GetCommunicator().LocalMesh().NodesBegin();
         i != rMeshModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++i) {

      array_1d<double, 3> &mesh_v =
          (i)->FastGetSolutionStepValue(MESH_VELOCITY);
      array_1d<double, 3> &disp =
          (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
      array_1d<double, 3> &dispold =
          (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
      noalias(mesh_v) = disp - dispold;
      mesh_v *= coeff;
    }
  } else if (TimeOrder == 2) {
    const double c1 = 1.50 * coeff;
    const double c2 = -2.0 * coeff;
    const double c3 = 0.50 * coeff;

    for (ModelPart::NodeIterator i =
             rMeshModelPart.GetCommunicator().LocalMesh().NodesBegin();
         i != rMeshModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++i) {

      array_1d<double, 3> &mesh_v =
          (i)->FastGetSolutionStepValue(MESH_VELOCITY);
      noalias(mesh_v) = c1 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
      noalias(mesh_v) +=
          c2 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
      noalias(mesh_v) +=
          c3 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
    }
  } else {
    KRATOS_ERROR << "Wrong TimeOrder: Acceptable values are: 1 and 2"
                 << std::endl;
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void MoveMesh(const ModelPart::NodesContainerType &rNodes) {
  KRATOS_TRY;

  for (auto &rnode : rNodes) { // TODO OMP
    noalias(rnode.Coordinates()) = rnode.GetInitialPosition()
                     + rnode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void SetMeshToInitialConfiguration(
    const ModelPart::NodesContainerType &rNodes) {
  KRATOS_TRY;

  for (auto &rnode : rNodes) { // TODO OMP
    noalias(rnode.Coordinates()) = rnode.GetInitialPosition();
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
  void UpdateReferenceMesh(ModelPart &rModelPart) // TODO is this still needed?
{

  for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
          i != rModelPart.NodesEnd(); ++i){

      (i)->X0() = (i)->X();
      (i)->Y0() = (i)->Y();
      (i)->Z0() = (i)->Z();

  }
}

} // namespace Move Mesh Utilities.

} // namespace Kratos.