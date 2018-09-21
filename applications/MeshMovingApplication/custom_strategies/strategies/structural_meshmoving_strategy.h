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

#if !defined(KRATOS_STRUCTURAL_MESHMOVING_STRATEGY_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MESHMOVING_STRATEGY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "base_meshmoving_strategy.h"


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
class StructuralMeshMovingStrategy : public BaseMeshMovingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuralMeshMovingStrategy
    KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef BaseType::Pointer BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuralMeshMovingStrategy(ModelPart& rMeshMovingModelPart,
                                Parameters MeshVolocityCalculationParameters,
                                const bool ReformDofSetAtEachStep,
                                BaseTypePointer pStrategy) :
        BaseMeshMovingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rMeshMovingModelPart,
            MeshVolocityCalculationParameters,
            ReformDofSetAtEachStep),
            mpStrategy(pStrategy)
    { }

    /// Destructor.
    ~StructuralMeshMovingStrategy() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        mpStrategy->Initialize();
    }

    void Clear() override
    {
        mpStrategy->Clear();
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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void ComputeMeshMovement() override
    {
        mpStrategy->Solve();
    }

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

    BaseTypePointer pStrategy;

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
    StructuralMeshMovingStrategy& operator=(StructuralMeshMovingStrategy const& rOther){}

    /// Copy constructor.
    StructuralMeshMovingStrategy(StructuralMeshMovingStrategy const& rOther){}


    ///@}

}; // Class StructuralMeshMovingStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MESHMOVING_STRATEGY_H_INCLUDED  defined
