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
#include <string>
#include <iostream>


// External includes


// Project includes
#include "solving_strategies/strategies/solving_strategy.h"


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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseMeshMovingStrategy(){}

    /// Destructor.
    virtual ~BaseMeshMovingStrategy(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
        buffer << "BaseMeshMovingStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "BaseMeshMovingStrategy";}

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


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                BaseMeshMovingStrategy& rThis){}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const BaseMeshMovingStrategy& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BASE_MESHMOVING_STRATEGY_H_INCLUDED  defined


