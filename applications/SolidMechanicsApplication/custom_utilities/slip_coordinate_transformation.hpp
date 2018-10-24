//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             October 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#ifndef KRATOS_SLIP_COORDINATE_TRANSFORMATION_H_INCLUDED
#define KRATOS_SLIP_COORDINATE_TRANSFORMATION_H_INCLUDED

// system includes

// external includes

// kratos includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/variable.h"
#include "geometries/geometry.h"

namespace Kratos {

///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// A utility to rotate the local contributions of certain nodes to the system matrix, which is required to apply slip conditions in arbitrary directions.
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
class SlipCoordinateTransformation
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of SlipCoordinateTransformation
  KRATOS_CLASS_POINTER_DEFINITION(SlipCoordinateTransformation);

  typedef Node<3>                     NodeType;
  typedef Geometry<NodeType>      GeometryType;
  typedef GeometryType::PointType    PointType;
  typedef GeometryType::SizeType      SizeType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  SlipCoordinateTransformation() {}

  /// Constructor.
  /** @param NumRowsPerNode Number of matrix or vector rows associated to each node. Velocity DOFs are assumed to be the first rows (with dimension size) in each block of rows.
   */
  SlipCoordinateTransformation(const unsigned int NumRowsPerNode) : mBlockSize(NumRowsPerNode)
  {}

  /// Destructor.
  virtual ~SlipCoordinateTransformation() {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void SetBlockSize(const unsigned int NumRowsPerNode)
  {
    mBlockSize = NumRowsPerNode;
  }

  /// Rotate the local system contributions so that they are oriented with each node's normal.
  /**
     @param rLocalMatrix Local system matrix
     @param rLocalVector Local RHS vector
     @param rGeometry A reference to the element's (or condition's) geometry
  */
  virtual void Rotate(TLocalMatrixType& rLocalMatrix,
                      TLocalVectorType& rLocalVector,
                      GeometryType& rGeometry) const
  {
    const SizeType dimension = rGeometry.WorkingSpaceDimension();

    if(mBlockSize != dimension) //Monolithic case
    {
      if(dimension == 2) RotateAux<2,3>(rLocalMatrix,rLocalVector,rGeometry);
      if(dimension == 3) RotateAux<3,4>(rLocalMatrix,rLocalVector,rGeometry);
    }
    else //fractional step case
    {
      if(dimension == 2) RotateAuxPure<2>(rLocalMatrix,rLocalVector,rGeometry);
      if(dimension == 3) RotateAuxPure<3>(rLocalMatrix,rLocalVector,rGeometry);
    }

  }

  /// RHS only version of Rotate
  virtual void Rotate(TLocalVectorType& rLocalVector,
                      GeometryType& rGeometry) const
  {
    //const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

    unsigned int Index = 0;
    const SizeType dimension = rGeometry.WorkingSpaceDimension();

    if (rLocalVector.size() > 0)
    {
      if(mBlockSize != dimension) //Monolithic case
      {
        for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
        {
          if( rGeometry[j].Is(SLIP) )
          {
            if(dimension == 3)
            {
              array_1d<double,4> aux,aux1;
              BoundedMatrix<double,4,4> rRot;
              LocalRotationOperator3D<4>(rRot,rGeometry[j]);

              for(unsigned int k=0; k<4; k++)
                aux[k] = rLocalVector[j*mBlockSize+k];

              noalias(aux1) = prod(rRot,aux);

              for(unsigned int k=0; k<4; k++)
                rLocalVector[j*mBlockSize+k] = aux1[k];
            }
            else
            {
              array_1d<double,3> aux,aux1;
              BoundedMatrix<double,3,3> rRot;
              LocalRotationOperator2D<3>(rRot,rGeometry[j]);

              for(unsigned int k=0; k<3; k++)
              {
                aux[k] = rLocalVector[j*mBlockSize+k];
              }

              noalias(aux1) = prod(rRot,aux);

              for(unsigned int k=0; k<3; k++)
                rLocalVector[j*mBlockSize+k] = aux1[k];
            }
          }
          Index += mBlockSize;
        }

      }
      else //fractional step case
      {
        for(unsigned int j = 0; j < rGeometry.PointsNumber(); ++j)
        {
          if( rGeometry[j].Is(SLIP) )
          {
            if(dimension == 3)
            {
              array_1d<double,3> aux,aux1;
              BoundedMatrix<double,3,3> rRot;
              LocalRotationOperatorPure(rRot,rGeometry[j]);

              for(unsigned int k=0; k<3; k++)
                aux[k] = rLocalVector[j*mBlockSize+k];

              noalias(aux1) = prod(rRot,aux);

              for(unsigned int k=0; k<3; k++)
                rLocalVector[j*mBlockSize+k] = aux1[k];
            }
            else
            {
              array_1d<double,2> aux,aux1;
              BoundedMatrix<double,2,2> rRot;
              LocalRotationOperatorPure(rRot,rGeometry[j]);

              for(unsigned int k=0; k<2; k++)
                aux[k] = rLocalVector[j*mBlockSize+k];

              noalias(aux1) = prod(rRot,aux);

              for(unsigned int k=0; k<2; k++)
                rLocalVector[j*mBlockSize+k] = aux1[k];
            }
          }
          Index += mBlockSize;
        }

      }

    }

  }

  /// Apply slip boundary conditions to the rotated local contributions.
  /** This function takes the local system contributions rotated so each
      node's velocities are expressed using a base oriented with its normal
      and imposes that the normal velocity is equal to the mesh velocity in
      the normal direction.
  */
  virtual void ApplySlipCondition(TLocalMatrixType& rLocalMatrix,
                                  TLocalVectorType& rLocalVector,
                                  GeometryType& rGeometry) const
  {
    const unsigned int LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

    if (LocalSize > 0)
    {
      for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
      {
        if( rGeometry[itNode].Is(SLIP) )
        {
          // We fix the first dof (normal velocity) for each rotated block
          unsigned int j = itNode * mBlockSize;
          //const double k = rLocalMatrix(j,j)+rLocalMatrix(j,j+1)+rLocalMatrix(j,j+2);

          // If the mesh is moving, we must impose v_normal = vmesh_normal
          array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
          VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
          array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
          this->Normalize(rN);

          for( unsigned int i = 0; i < j; ++i)// Skip term (i,i)
          {
            rLocalMatrix(i,j) = 0.0;
            rLocalMatrix(j,i) = 0.0;
          }
          for( unsigned int i = j+1; i < LocalSize; ++i)
          {
            rLocalMatrix(i,j) = 0.0;
            rLocalMatrix(j,i) = 0.0;
          }

          rLocalVector(j) = inner_prod(rN,VMesh);
          rLocalMatrix(j,j) = 1.0;
        }
      }
    }
  }

  /// RHS only version of ApplySlipCondition
  virtual void ApplySlipCondition(TLocalVectorType& rLocalVector,
                                  GeometryType& rGeometry) const
  {
    if (rLocalVector.size() > 0)
    {
      for(unsigned int itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
      {
        if( rGeometry[itNode].Is(SLIP) )
        {
          // We fix the first dof (normal velocity) for each rotated block
          unsigned int j = itNode * mBlockSize;

          // If the mesh is moving, we must impose v_normal = vmesh_normal
          array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
          VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
          array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
          this->Normalize(rN);

          rLocalVector[j] = inner_prod(rN,VMesh);
        }
      }
    }
  }

  /// Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
  virtual void RotateVelocities(ModelPart& rModelPart) const
  {
    const SizeType dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    TLocalVectorType Vel(dimension);
    TLocalVectorType Tmp(dimension);

    ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
    for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
    {
      ModelPart::NodeIterator itNode = it_begin+iii;
      if( itNode->Is(SLIP) )
      {
        //this->RotationOperator<TLocalMatrixType>(Rotation,);
        if(dimension == 3)
        {
          BoundedMatrix<double,3,3> rRot;
          LocalRotationOperatorPure(rRot,*itNode);

          array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
          for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
          noalias(Tmp) = prod(rRot,Vel);
          for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
        }
        else
        {
          BoundedMatrix<double,2,2> rRot;
          LocalRotationOperatorPure(rRot,*itNode);

          array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
          for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
          noalias(Tmp) = prod(rRot,Vel);
          for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
        }
      }
    }
  }

  /// Transform nodal velocities from the rotated system to the original one
  virtual void RecoverVelocities(ModelPart& rModelPart) const
  {
    const SizeType dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    TLocalVectorType Vel(dimension);
    TLocalVectorType Tmp(dimension);

    ModelPart::NodeIterator it_begin = rModelPart.NodesBegin();
#pragma omp parallel for firstprivate(Vel,Tmp)
    for(int iii=0; iii<static_cast<int>(rModelPart.Nodes().size()); iii++)
    {
      ModelPart::NodeIterator itNode = it_begin+iii;
      if( itNode->Is(SLIP) )
      {
        if(dimension == 3)
        {
          BoundedMatrix<double,3,3> rRot;
          LocalRotationOperatorPure(rRot,*itNode);

          array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
          for(unsigned int i = 0; i < 3; i++) Vel[i] = rVelocity[i];
          noalias(Tmp) = prod(trans(rRot),Vel);
          for(unsigned int i = 0; i < 3; i++) rVelocity[i] = Tmp[i];
        }
        else
        {
          BoundedMatrix<double,2,2> rRot;
          LocalRotationOperatorPure(rRot,*itNode);

          array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
          for(unsigned int i = 0; i < 2; i++) Vel[i] = rVelocity[i];
          noalias(Tmp) = prod(trans(rRot),Vel);
          for(unsigned int i = 0; i < 2; i++) rVelocity[i] = Tmp[i];
        }
      }
    }
  }

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const
  {
    std::stringstream buffer;
    buffer << "SlipCoordinateTransformation";
    return buffer.str();
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "SlipCoordinateTransformation";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {}

  ///@}
  ///@name Friends
  ///@{

  ///@}

 protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  template<unsigned int TDim, unsigned int TBlockSize, unsigned int TSkip = 0>
  void RotateAux(TLocalMatrixType& rLocalMatrix, TLocalVectorType& rLocalVector, GeometryType& rGeometry) const
  {
    const unsigned int LocalSize = rLocalVector.size();

    unsigned int Index = 0;
    int rotations_needed = 0;
    const unsigned int NumBlocks = LocalSize / TBlockSize;
    DenseVector<bool> NeedRotation( NumBlocks, false);

    std::vector< BoundedMatrix<double,TBlockSize,TBlockSize> > rRot(NumBlocks);
    for(unsigned int j = 0; j < NumBlocks; ++j)
    {
      if( rGeometry[j].Is(SLIP) )
      {
        NeedRotation[j] = true;
        rotations_needed++;

        if (TDim == 2) LocalRotationOperator2D<TBlockSize,TSkip>(rRot[j],rGeometry[j]);
        else LocalRotationOperator3D<TBlockSize,TSkip>(rRot[j],rGeometry[j]);
      }

      Index += TBlockSize;
    }

    if(rotations_needed > 0)
    {
      BoundedMatrix<double,TBlockSize,TBlockSize> mat_block, tmp;
      array_1d<double,TBlockSize> aux, aux1;

      for(unsigned int i=0; i<NumBlocks; i++)
      {
        if(NeedRotation[i] == true)
        {
          for(unsigned int j=0; j<NumBlocks; j++)
          {
            if(NeedRotation[j] == true)
            {
              ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
              noalias(tmp) = prod(mat_block,trans(rRot[j]));
              noalias(mat_block) = prod(rRot[i],tmp);
              WriteBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
            }
            else
            {
              ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
              noalias(tmp) = prod(rRot[i],mat_block);
              WriteBlockMatrix<TBlockSize>(tmp, rLocalMatrix, i*TBlockSize, j*TBlockSize);
            }
          }

          for(unsigned int k=0; k<TBlockSize; k++)
            aux[k] = rLocalVector[i*TBlockSize+k];

          noalias(aux1) = prod(rRot[i],aux);

          for(unsigned int k=0; k<TBlockSize; k++)
            rLocalVector[i*TBlockSize+k] = aux1[k];

        }
        else
        {
          for(unsigned int j=0; j<NumBlocks; j++)
          {
            if(NeedRotation[j] == true)
            {
              ReadBlockMatrix<TBlockSize>(mat_block, rLocalMatrix, i*TBlockSize, j*TBlockSize);
              noalias(tmp) = prod(mat_block,trans(rRot[j]));
              WriteBlockMatrix<TBlockSize>(tmp, rLocalMatrix, i*TBlockSize, j*TBlockSize);
            }
          }
        }

      }
    }
  }

  //to be used when there is only velocity (no additional pressure or other var block)
  template<unsigned int TDim>
  void RotateAuxPure(TLocalMatrixType& rLocalMatrix, TLocalVectorType& rLocalVector, GeometryType& rGeometry) const
  {
    const unsigned int LocalSize = rLocalVector.size();

    unsigned int Index = 0;
    int rotations_needed = 0;
    const unsigned int NumBlocks = LocalSize / mBlockSize;
    DenseVector<bool> NeedRotation( NumBlocks, false);

    std::vector< BoundedMatrix<double,TDim,TDim> > rRot(NumBlocks);
    for(unsigned int j = 0; j < NumBlocks; ++j)
    {
      if( rGeometry[j].Is(SLIP) )
      {
        NeedRotation[j] = true;
        rotations_needed++;

        LocalRotationOperatorPure(rRot[j],rGeometry[j]);
      }

      Index += mBlockSize;
    }

    if(rotations_needed > 0)
    {
      BoundedMatrix<double,TDim,TDim> mat_block, tmp;
      array_1d<double,TDim> aux, aux1;

      for(unsigned int i=0; i<NumBlocks; i++)
      {
        if(NeedRotation[i] == true)
        {
          for(unsigned int j=0; j<NumBlocks; j++)
          {
            if(NeedRotation[j] == true)
            {
              ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
              noalias(tmp) = prod(mat_block,trans(rRot[j]));
              noalias(mat_block) = prod(rRot[i],tmp);
              WriteBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
            }
            else
            {
              ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
              noalias(tmp) = prod(rRot[i],mat_block);
              WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
            }
          }

          for(unsigned int k=0; k<TDim; k++)
            aux[k] = rLocalVector[i*mBlockSize+k];

          noalias(aux1) = prod(rRot[i],aux);

          for(unsigned int k=0; k<TDim; k++)
            rLocalVector[i*mBlockSize+k] = aux1[k];

        }
        else
        {
          for(unsigned int j=0; j<NumBlocks; j++)
          {
            if(NeedRotation[j] == true)
            {
              ReadBlockMatrix<TDim>(mat_block, rLocalMatrix, i*mBlockSize, j*mBlockSize);
              noalias(tmp) = prod(mat_block,trans(rRot[j]));
              WriteBlockMatrix<TDim>(tmp, rLocalMatrix, i*mBlockSize, j*mBlockSize);
            }
          }
        }

      }
    }
  }

  template<unsigned int TBlockSize, unsigned int TSkip = 0>
  void LocalRotationOperator2D(BoundedMatrix<double,TBlockSize,TBlockSize>& rRot, const PointType& rThisPoint) const
  {
    noalias(rRot) = IdentityMatrix(TBlockSize);

    // Get the normal evaluated at the node
    const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

    double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
    aux = sqrt(aux);

    rRot(TSkip  ,TSkip  ) = rNormal[0]/aux;
    rRot(TSkip  ,TSkip+1) = rNormal[1]/aux;
    rRot(TSkip+1,TSkip  ) = -rNormal[1]/aux;
    rRot(TSkip+1,TSkip+1) = rNormal[0]/aux;
  }

  template<unsigned int TBlockSize, unsigned int TSkip = 0>
  void LocalRotationOperator3D(BoundedMatrix<double,TBlockSize,TBlockSize>& rRot, const PointType& rThisPoint) const
  {
    noalias(rRot) = IdentityMatrix(TBlockSize);

    // Get the normal evaluated at the node
    const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

    double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
    aux = sqrt(aux);
    rRot(TSkip,TSkip  ) = rNormal[0]/aux;
    rRot(TSkip,TSkip+1) = rNormal[1]/aux;
    rRot(TSkip,TSkip+2) = rNormal[2]/aux;
    // Define the new coordinate system, where the first vector is aligned with the normal

    // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
    array_1d<double,3> rT1;
    rT1(0) = 1.0;
    rT1(1) = 0.0;
    rT1(2) = 0.0;
    double dot = rRot(TSkip,TSkip);//this->Dot(rN,rT1);

    // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
    // If this is the case, repeat the procedure using (0,1,0)
    if ( fabs(dot) > 0.99 )
    {
      rT1(0) = 0.0;
      rT1(1) = 1.0;
      rT1(2) = 0.0;

      dot = rRot(TSkip,TSkip+1); //this->Dot(rN,rT1);
    }

    // calculate projection and normalize
    rT1[0] -= dot*rRot(TSkip,TSkip);
    rT1[1] -= dot*rRot(TSkip,TSkip+1);
    rT1[2] -= dot*rRot(TSkip,TSkip+2);
    this->Normalize(rT1);
    rRot(TSkip+1,TSkip  ) = rT1[0];
    rRot(TSkip+1,TSkip+1) = rT1[1];
    rRot(TSkip+1,TSkip+2) = rT1[2];

    // The third base component is choosen as N x T1, which is normalized by construction
    rRot(TSkip+2,TSkip  ) = rRot(TSkip,TSkip+1)*rT1[2] - rRot(TSkip,TSkip+2)*rT1[1];
    rRot(TSkip+2,TSkip+1) = rRot(TSkip,TSkip+2)*rT1[0] - rRot(TSkip,TSkip  )*rT1[2];
    rRot(TSkip+2,TSkip+2) = rRot(TSkip,TSkip  )*rT1[1] - rRot(TSkip,TSkip+1)*rT1[0];
  }


  void LocalRotationOperatorPure(BoundedMatrix<double,3,3>& rRot, const PointType& rThisPoint) const
  {

    // Get the normal evaluated at the node
    const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

    double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1] + rNormal[2]*rNormal[2];
    aux = sqrt(aux);
    rRot(0,0) = rNormal[0]/aux;
    rRot(0,1) = rNormal[1]/aux;
    rRot(0,2) = rNormal[2]/aux;
    // Define the new coordinate system, where the first vector is aligned with the normal

    // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
    array_1d<double,3> rT1;
    rT1(0) = 1.0;
    rT1(1) = 0.0;
    rT1(2) = 0.0;
    double dot = rRot(0,0);//this->Dot(rN,rT1);

    // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
    // If this is the case, repeat the procedure using (0,1,0)
    if ( fabs(dot) > 0.99 )
    {
      rT1(0) = 0.0;
      rT1(1) = 1.0;
      rT1(2) = 0.0;

      dot = rRot(0,1); //this->Dot(rN,rT1);
    }

    // calculate projection and normalize
    rT1[0] -= dot*rRot(0,0);
    rT1[1] -= dot*rRot(0,1);
    rT1[2] -= dot*rRot(0,2);
    this->Normalize(rT1);
    rRot(1,0) = rT1[0];
    rRot(1,1) = rT1[1];
    rRot(1,2) = rT1[2];

    // The third base component is choosen as N x T1, which is normalized by construction
    rRot(2,0) = rRot(0,1)*rT1[2] - rRot(0,2)*rT1[1];
    rRot(2,1) = rRot(0,2)*rT1[0] - rRot(0,0)*rT1[2];
    rRot(2,2) = rRot(0,0)*rT1[1] - rRot(0,1)*rT1[0];
  }

  void LocalRotationOperatorPure(BoundedMatrix<double,2,2>& rRot, PointType& rThisPoint) const
  {
    // Get the normal evaluated at the node
    const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

    double aux = rNormal[0]*rNormal[0] + rNormal[1]*rNormal[1];
    aux = sqrt(aux);

    rRot(0,0) = rNormal[0]/aux;
    rRot(0,1) = rNormal[1]/aux;
    rRot(1,0) = -rNormal[1]/aux;
    rRot(1,1) = rNormal[0]/aux;

  }

  /// Normalize a vector.
  /**
   * @param rThis the vector
   * @return Original norm of the input vector
   */
  template< class TVectorType >
  double Normalize(TVectorType& rThis) const
  {
    double Norm = 0;
    for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
      Norm += (*iComponent)*(*iComponent);
    Norm = sqrt(Norm);
    for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
      *iComponent /= Norm;
    return Norm;
  }

  ///@}
  ///@name Protected  Access
  ///@{

  unsigned int GetBlockSize() const
  {
    return mBlockSize;
  }

  ///@}
  ///@name Protected Inquiry
  ///@{

  ///@}
  ///@name Protected LifeCycle
  ///@{

  ///@}

 private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  /// Number of matrix or vector rows associated to each node.
  /** @note Velocity Dofs are assumed to be the first rows (with dimension size).
   */
  const unsigned int mBlockSize;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  template< class TVectorType >
  double Dot(const TVectorType& rV1, const TVectorType& rV2) const
  {
    double dot = 0.0;
    for( typename TVectorType::const_iterator iV1 = rV1.begin(),iV2 = rV2.begin(); iV1 != rV1.end(); ++iV1, ++iV2)
    {
      dot += (*iV1) * (*iV2);
    }
    return dot;
  }

  //auxiliary functions
  template< unsigned int TBlockSize >
  void ReadBlockMatrix(BoundedMatrix<double,TBlockSize, TBlockSize>& block, const Matrix& origin, const unsigned int i_begin, const unsigned int j_begin) const
  {
    for(unsigned int i=0; i<TBlockSize; i++)
    {
      for(unsigned int j=0; j<TBlockSize; j++)
      {
        block(i,j) = origin(i_begin+i, j_begin+j);
      }
    }
  }

  template< unsigned int TBlockSize >
  void WriteBlockMatrix( const BoundedMatrix<double,TBlockSize, TBlockSize>& block, Matrix& destination, const unsigned int i_begin, const unsigned int j_begin) const
  {
    for(unsigned int i=0; i<TBlockSize; i++)
    {
      for(unsigned int j=0; j<TBlockSize; j++)
      {
        destination(i_begin+i, j_begin+j) = block(i,j);
      }
    }
  }

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  /// Assignment operator.
  SlipCoordinateTransformation& operator=(SlipCoordinateTransformation const& rOther) {}

  /// Copy constructor.
  SlipCoordinateTransformation(SlipCoordinateTransformation const& rOther) {}

  ///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::istream& operator >>(std::istream& rIStream,
                                 SlipCoordinateTransformation<TLocalMatrixType, TLocalVectorType,
                                 TValueType>& rThis)
{
  return rIStream;
}

/// output stream function
template<class TLocalMatrixType, class TLocalVectorType, class TValueType>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SlipCoordinateTransformation<TLocalMatrixType, TLocalVectorType,
                                 TValueType>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

///@} addtogroup block

}

#endif // KRATOS_SLIP_COORDINATE_TRANSFORMATION_H_INCLUDED
