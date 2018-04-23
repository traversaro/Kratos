// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//   

#ifndef ADJOINT_INTEGRAL_STRESS_RESPONSE_FUNCTION_H
#define ADJOINT_INTEGRAL_STRESS_RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "adjoint_structural_response_function.h"
#include "response_data.h"


// ==============================================================================

namespace Kratos
{

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

//template<class TDenseSpace>

class AdjointIntegralStressResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef AdjointStructuralResponseFunction BaseType;
    typedef array_1d<double, 3> array_3d;



    /// Pointer definition of AdjointIntegralStressResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointIntegralStressResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    AdjointIntegralStressResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
    : AdjointStructuralResponseFunction(rModelPart, rParameters)
    {
        ModelPart& r_model_part = this->GetModelPart();

        ResponseData stress_response_data;

        // Get name of model part which contains all the elements which should be traced
        mTracedElementsModelPartName = rParameters["traced_elements_model_part_name"].GetString();
        KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mTracedElementsModelPartName))
            << "No sub model part \"" << mTracedElementsModelPartName << "\"" << std::endl;

        // Tell traced elements the stress type
        TracedStressType traced_stress_type = stress_response_data.ConvertStressType(rParameters["stress_type"].GetString()); 
        KRATOS_ERROR_IF(traced_stress_type == StressTypeNotAvailible) << "Chosen stress type is not availible!" << std::endl;	

        KRATOS_ERROR_IF(r_model_part.GetSubModelPart(mTracedElementsModelPartName).Elements().size() == 0)  
            << "There is no element in  \"" << mTracedElementsModelPartName << "\"" << std::endl;

        mNumberOfTracedElements = r_model_part.GetSubModelPart(mTracedElementsModelPartName).Elements().size();

    #pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(
                r_model_part.GetSubModelPart(mTracedElementsModelPartName).Elements(),
                elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                it->SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress_type) );	

        }

        mStressValue = 0.0;
    }

    /// Destructor.
    virtual ~AdjointIntegralStressResponseFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        BaseType::Initialize();



        KRATOS_CATCH("");
    }

    // ==============================================================================
    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        // Working variables
        ProcessInfo &r_current_precess_info = rModelPart.GetProcessInfo();
        Vector element_stress;
        ModelPart& r_traced_elements_model_part =  this->GetModelPart().GetSubModelPart(mTracedElementsModelPartName);

        // Loop over traced elements and compute mean value of traced stress
        for (auto& elem_i : r_traced_elements_model_part.Elements())
        {
            elem_i.Calculate(STRESS_ON_GP, element_stress, r_current_precess_info);
       
            int stress_vec_size = element_stress.size();

            for(int i = 0; i < stress_vec_size; i++)
                mStressValue += element_stress[i];

            mStressValue /= stress_vec_size;
        } 

        mStressValue /= mNumberOfTracedElements;   
    
        return mStressValue;

        KRATOS_CATCH("");
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

    // ==============================================================================
    void CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
    {
        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        Element copy_of_element = rAdjointElem;

        if(rAdjointElem.Has(TRACED_STRESS_TYPE)) // maybe not the best way to check if there are more than one response functions in the future
        {
            Matrix stress_displ_deriv;
            copy_of_element.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displ_deriv, rProcessInfo);
 
            int num_of_dofs = stress_displ_deriv.size1();
            int num_of_deriv = stress_displ_deriv.size2();
            double stress_displ_deriv_value = 0.0;

            KRATOS_ERROR_IF(rResponseGradient.size() != stress_displ_deriv.size1())
                 << "Size of stress displacement derivative does not fit!" << std::endl;

            for (int dof_it = 0 ; dof_it < num_of_dofs; dof_it++)
            {
                for(int GP_it = 0; GP_it < num_of_deriv; GP_it++)
                    stress_displ_deriv_value += stress_displ_deriv(dof_it, GP_it);

                stress_displ_deriv_value /= (num_of_deriv * mNumberOfTracedElements);

                rResponseGradient[dof_it] = (-1) * stress_displ_deriv_value;
                stress_displ_deriv_value = 0.0;
            }
        }
    }


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

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY

        if(rAdjointElem.Has(TRACED_STRESS_TYPE)) // maybe not the best way to check if there are more than one response functions in the future
        {
            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, rVariable.Name());

            Matrix stress_DV_deriv;
            rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_DV_deriv, rProcessInfo);
 
            int num_of_DV = stress_DV_deriv.size1();
            int num_of_deriv = stress_DV_deriv.size2();
            double stress_DV_deriv_value = 0.0;

            if(rResponseGradient.size() != stress_DV_deriv.size1())
                rResponseGradient.resize(stress_DV_deriv.size1(), false);
            KRATOS_ERROR_IF(rResponseGradient.size() != rDerivativesMatrix.size1())
                 << "Size of partial stress design variable derivative does not fit!" << std::endl;

            for (int dv_it = 0 ; dv_it < num_of_DV; dv_it++)
            {
                for(int GP_it = 0; GP_it < num_of_deriv; GP_it++)
                    stress_DV_deriv_value += stress_DV_deriv(dv_it, GP_it);

                stress_DV_deriv_value /= (num_of_deriv * mNumberOfTracedElements);
            
                rResponseGradient[dv_it] =  stress_DV_deriv_value;
                stress_DV_deriv_value = 0.0;
            }

            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, "");
        }
        else
        {
            if (rResponseGradient.size() != rDerivativesMatrix.size1())
                      rResponseGradient.resize(rDerivativesMatrix.size1(), false);
            rResponseGradient.clear();
        }

        KRATOS_CATCH("")
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
                  rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if(rAdjointElem.Has(TRACED_STRESS_TYPE)) // maybe not the best way to check if there are more than one response functions in the future
        {
            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, rVariable.Name());

            Matrix stress_DV_deriv;
            rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_DV_deriv, rProcessInfo);

            int num_of_DV = stress_DV_deriv.size1();
            int num_of_deriv = stress_DV_deriv.size2();
            double stress_DV_deriv_value = 0.0;

            if(rResponseGradient.size() != stress_DV_deriv.size1())
                rResponseGradient.resize(stress_DV_deriv.size1(), false);
            KRATOS_ERROR_IF(rResponseGradient.size() != rDerivativesMatrix.size1())
                << "Size of partial stress design variable derivative does not fit!" << std::endl;

            for (int dv_it = 0 ; dv_it < num_of_DV; dv_it++)
            {
                for(int GP_it = 0; GP_it < num_of_deriv; GP_it++)
                    stress_DV_deriv_value += stress_DV_deriv(dv_it, GP_it);

                stress_DV_deriv_value /= (num_of_deriv * mNumberOfTracedElements);

                rResponseGradient[dv_it] = stress_DV_deriv_value;
                stress_DV_deriv_value = 0.0;
            }

            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, "");
        }
        else
        {
            if (rResponseGradient.size() != rDerivativesMatrix.size1())
                      rResponseGradient.resize(rDerivativesMatrix.size1(), false);
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if(rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    // ==============================================================================

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

    double mStressValue;
    std::string mTracedElementsModelPartName;
    int mNumberOfTracedElements;

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
    //      AdjointIntegralStressResponseFunction& operator=(AdjointIntegralStressResponseFunction const& rOther);

    /// Copy constructor.
    //      AdjointIntegralStressResponseFunction(AdjointIntegralStressResponseFunction const& rOther);

    ///@}

}; // Class AdjointIntegralStressResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_INTEGRAL_STRESS_RESPONSE_FUNCTION_H
