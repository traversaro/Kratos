//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein
//                   Philipp Bucher
//

#if !defined(KRATOS_BASE_MESHMOVING_STRATEGY_H_INCLUDED )
#define  KRATOS_BASE_MESHMOVING_STRATEGY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/move_mesh_utilities.h"
#include "custom_utilities/mesh_velocity_calculation_utility.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class BaseMeshMovingStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BaseMeshMovingStrategy
    KRATOS_CLASS_POINTER_DEFINITION(BaseMeshMovingStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef Kratos::unique_ptr<MeshVelocityCalculationUtility> MeshVelCalcUtilPtrType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseMeshMovingStrategy(ModelPart& rMeshMovingModelPart,
                           Parameters MeshVolocityCalculationParameters,
                           const bool ReformDofSetAtEachStep) :
        SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rMeshMovingModelPart),
        mReformDofSetAtEachStep(ReformDofSetAtEachStep)
    {
        const bool mCalculateMeshVelocities =
            MeshVolocityCalculationParameters["calculate_mesh_velocities"].GetBool()

        if (mCalculateMeshVelocities) {
            mpMeshVelocityCalculationUtility = Kratos::make_unique<MeshVelocityCalculationUtility>(
                BaseType::GetModelPart(),
                MeshVolocityCalculationParameters);
        }
    }

    /// Destructor.
    ~BaseMeshMovingStrategy() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    bool SolveSolutionStep() override
    {
        MoveMeshUtilities::SetMeshToInitialConfiguration(
            BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

        // here the Mesh-Movement-Method specific operations are performed
        ComputeMeshMovement();

        if (mCalculateMeshVelocities) {
            KRATOS_DEBUG_ERROR_IF_NOT(mpMeshVelocityCalculationUtility)
                << "mpMeshVelocityCalculationUtility is a nullptr!" << std::endl;
            mpMeshVelocityCalculationUtility->CalculateMeshVelocities();
        }

        MoveMeshUtilities::MoveMesh(
            BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

        if (mReformDofSetAtEachStep) {
            Clear();
        }

        return true;
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

    MeshVelCalcUtilPtrType mpMeshVelocityCalculationUtility = nullptr;
    bool mReformDofSetAtEachStep;
    bool mCalculateMeshVelocities;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    // here the Mesh-Movement-Method specific operations are performed
    virtual void ComputeMeshMovement() = 0;


    ///@}
    ///@name Protected  Access
    ///@{


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
    BaseMeshMovingStrategy& operator=(BaseMeshMovingStrategy const& rOther){}

    /// Copy constructor.
    BaseMeshMovingStrategy(BaseMeshMovingStrategy const& rOther){}


    ///@}

}; // Class BaseMeshMovingStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BASE_MESHMOVING_STRATEGY_H_INCLUDED  defined


