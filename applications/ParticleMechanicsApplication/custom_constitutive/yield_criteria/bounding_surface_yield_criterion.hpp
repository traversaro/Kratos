//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_BOUNDING_SURFACE_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_BOUNDING_SURFACE_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/yield_criteria/MPM_yield_criterion.hpp"
#include "custom_constitutive/hardening_laws/bounding_surface_hardening_law.hpp"

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
class BoundingSurfaceYieldCriterion
	: public MPMYieldCriterion
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MisesHuberYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( BoundingSurfaceYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        BoundingSurfaceYieldCriterion();

        /// Copy constructor.
        BoundingSurfaceYieldCriterion(BoundingSurfaceYieldCriterion const& rOther);

        /// Initialization constructor.
        BoundingSurfaceYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Assignment operator.
        BoundingSurfaceYieldCriterion& operator=(BoundingSurfaceYieldCriterion const& rOther);


        /// Destructor.
        virtual ~BoundingSurfaceYieldCriterion();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /*
        * @brief This function calculate either Bounding Surface or Loading Surface depends on the
        *        given input rStressVector and rPreconsolidationPressure
        * @param[in/out] rStateFunction Bounding surface or loading surface
        * @param[in] rStressVector Principal stresses
        * @param[in] rPreconsolidationPressure The value of Preconsolidation Stress
        * @return Bounding surface or loading surface
        */
        double& CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rPreconsolidationPressure) override;


        /*
        * @brief This function return the first derivative of bounding/loading surface at the given principal stress condition
        * @param[in] rStressVector Principal stresses
        * @param[in/out] rFirstDerivative First stress derivative value of yield function
        */
        void CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative) override;

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        // /// Turn back information as a string.
        // std::string Info() const override;

        // /// Print information about this object.
        // void PrintInfo(std::ostream& rOStream) const override;

        // /// Print object's data.
        // void PrintData(std::ostream& rOStream) const override;


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{
        double GetAlphaParameter();

        double CalculateCriticalStateLineSlope(const double& rLodeAngle);

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

        double GetPI();

        ///@}
        ///@name Private  Access
        ///@{


	///@}
	///@name Serialization
	///@{
	friend class Serializer;

	// A private default constructor necessary for serialization

	void save(Serializer& rSerializer) const override;

	void load(Serializer& rSerializer) override;

        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    }; // Class MisesHuberYieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    // /// input stream function
    // inline std::istream& operator >> (std::istream& rIStream,
    //                                   MisesHuberYieldCriterion& rThis);

    // /// output stream function
    // inline std::ostream& operator << (std::ostream& rOStream,
    //                                   const MisesHuberYieldCriterion& rThis)
    // {
    //     rThis.PrintInfo(rOStream);
    //     rOStream << std::endl;
    //     rThis.PrintData(rOStream);

    //     return rOStream;
    // }
    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUNDING_SURFACE_YIELD_CRITERION_H_INCLUDED  defined


