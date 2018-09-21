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

#if !defined(KRATOS_LAPLACIAN_MESHMOVING_STRATEGY_H_INCLUDED )
#define  KRATOS_LAPLACIAN_MESHMOVING_STRATEGY_H_INCLUDED


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
class LaplacianMeshMovingStrategy : public LaplacianMeshMovingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LaplacianMeshMovingStrategy
    KRATOS_CLASS_POINTER_DEFINITION(LaplacianMeshMovingStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef BaseType::Pointer BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LaplacianMeshMovingStrategy(ModelPart& rMeshMovingModelPart,
                                Parameters MeshVolocityCalculationParameters,
                                const bool ReformDofSetAtEachStep,
                                BaseTypePointer pStrategyX)
                                BaseTypePointer pStrategyY)
                                BaseTypePointer pStrategyZ = nullptr) :
        LaplacianMeshMovingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rMeshMovingModelPart,
            MeshVolocityCalculationParameters,
            ReformDofSetAtEachStep),
            mpStrategyX(pStrategyX),
            mpStrategyY(pStrategyY),
            mpStrategyZ(pStrategyZ)
    { }

    /// Destructor.
    ~LaplacianMeshMovingStrategy() override {}

    /// Deleted copy constructor.
    LaplacianMeshMovingStrategy(LaplacianMeshMovingStrategy const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    LaplacianMeshMovingStrategy& operator=(LaplacianMeshMovingStrategy const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        mpStrategyX->Initialize();
        mpStrategyY->Initialize();
        if (mpStrategyZ) mpStrategyZ->Initialize();
    }

    void Clear() override
    {
        mpStrategyX->Clear();
        mpStrategyY->Clear();
        if (mpStrategyZ) mpStrategyZ->Clear();
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
        ProcessInfo& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        const std::size_t dimension = r_process_info[DOMAIN_SIZE];

        // X DIRECTION
        r_process_info[LAPLACIAN_DIRECTION] = 1;
        mpStrategyX->Solve();
        // Y DIRECTION
        r_process_info[LAPLACIAN_DIRECTION] = 2;
        mpStrategyY->Solve();

        if (dimension == 3) { // TODO use domain size to check or check if nullptr
            // Z DIRECTION
            r_process_info[LAPLACIAN_DIRECTION] = 3;
            mpStrategyZ->Solve();
        }
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

    BaseTypePointer pStrategyX;
    BaseTypePointer pStrategyY;
    BaseTypePointer pStrategyZ;

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


    ///@}

}; // Class LaplacianMeshMovingStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_MESHMOVING_STRATEGY_H_INCLUDED  defined
