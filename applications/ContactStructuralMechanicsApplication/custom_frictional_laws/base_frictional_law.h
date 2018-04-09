// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASE_FRICTIONAL_LAW_H_INCLUDED)
#define KRATOS_BASE_FRICTIONAL_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
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

/**
 * @class BaseFrictionalLaw
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This is the base class for any frictional law implemented in the ContactStructuralMechanicsApplication
 * @details Defines the common interface for all the frictional laws
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) BaseFrictionalLaw
    : public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BaseFrictionalLaw
    KRATOS_CLASS_POINTER_DEFINITION( BaseFrictionalLaw );

    /// The node definition
    typedef Node<3> NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BaseFrictionalLaw();

    /// Destructor.
    virtual ~BaseFrictionalLaw() {};

    /**
    * @brief Clone function (has to be implemented by any derived class)
    * @return a pointer to a new instance of this constitutive law
    */
    virtual BaseFrictionalLaw::Pointer Clone() const;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the current frictional coefficient for a given node
     * @param ThisNode The current node of analysis
     */
    virtual double GetFrictionalCoefficient(NodeType& ThisNode);

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
    //virtual std::string Info() const;

    /// Print information about this object.
    //virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //virtual void PrintData(std::ostream& rOStream) const;

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
    //BaseFrictionalLaw& operator=(BaseFrictionalLaw const& rOther);

    /// Copy constructor.
    //BaseFrictionalLaw(BaseFrictionalLaw const& rOther);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }

}; // Class BaseFrictionalLaw

///@}

///@name Type Definitions
///@{

/**
* Definition of BaseFrictionalLaw variable
*/
KRATOS_DEFINE_VARIABLE_IMPLEMENTATION( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, BaseFrictionalLaw::Pointer, FRICTIONAL_LAW )

///@}
///@name Input and output
///@{

} // namespace Kratos

#endif // KRATOS_BASE_FRICTIONAL_LAW_H_INCLUDED defined
