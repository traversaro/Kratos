// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
//

#if !defined(KRATOS_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED)
#define KRATOS_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;
    
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
 * @class HighCycleFatigueLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the isotropic damage theory by means of
 * linear/exponential softening.
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * The damage integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The fosftening behaviour considered (linear, exponential,etc...)
 * @tparam TYieldSurfaceType The yield surface considered
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class HighCycleFatigueLawIntegrator
{
public:
    ///@name Type Definitions
    ///@{


    /// Counted pointer of HighCycleFatigueLawIntegrator
    KRATOS_CLASS_POINTER_DEFINITION(HighCycleFatigueLawIntegrator);

    /// Initialization constructor
    HighCycleFatigueLawIntegrator()
    {
    }

    /// Copy constructor
    HighCycleFatigueLawIntegrator(HighCycleFatigueLawIntegrator const &rOther)
    {
    }

    /// Assignment operator
    HighCycleFatigueLawIntegrator &operator=(HighCycleFatigueLawIntegrator const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~HighCycleFatigueLawIntegrator()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void CalculateMaximumAndMinimumStresses(
        const double CurrentStress,
        double &MaximumStress, 
        double &MinimumStress, 
        const Vector &rPreviousStresses,
        unsigned int& rNumberOfCycles
        )
    {
        const double Stress1 = rPreviousStresses[1];
        const double Stress2 = rPreviousStresses[0];

        const double StressIncrement1 = Stress1 - Stress2;
        const double StressIncrement2 = CurrentStress - Stress1;

        if (StressIncrement1 >= 0.001 && StressIncrement2 <= 0.0) {
            MaximumStress = Stress1;
            rNumberOfCycles++;
        }
        if (StressIncrement1 <= 0.0 && StressIncrement2 >= 0.001) {
            MinimumStress = StressIncrement1;
        }
    }

    static void CalculateTensionCompressionFactor(const Vector& StressVector, double& rFactor) 
    {
        array_1d<double,3> principal_stresses;
        ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(principal_stresses, StressVector);


        double abs_component = 0.0, average_component = 0.0, sum_abs = 0.0, sum_average = 0.0;
        for (unsigned int i = 0; i < principal_stresses.size(); ++i) {
            abs_component = std::abs(principal_stresses[i]);
            average_component = 0.5 * (principal_stresses[i] + abs_component);
            sum_average += average_component;
            sum_abs += abs_component;
        }
        const double pre_indicator = sum_average / sum_abs;
        if (pre_indicator < 0.5) {
            rFactor = -1.0;
        } else {
            rFactor = 1.0;
        }
    }

    static void CalculateReversionFactor(const double MaxStress, const double MinStress, double& rReversionFactor)
    {
        rReversionFactor = MinStress / MaxStress;
    }

	static void CalculateFatigueReductionFactor(const double MaxStress,
                                                const double MinStress,
                                                double& rReversionFactor,)
	{

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

    ///@}

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif