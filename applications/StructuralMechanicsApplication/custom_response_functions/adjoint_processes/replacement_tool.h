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


#if !defined(KRATOS_REPLACEMENT_TOOL_H_INCLUDED)
#define  KRATOS_REPLACEMENT_TOOL_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{


class ReplacementTool
{
public:

    ///@name Type Definitions
    ///@{
    /// Pointer definition of ReplacementTool
    KRATOS_CLASS_POINTER_DEFINITION(ReplacementTool);
    ///@}
    ///@name Life Cycle
    ///@{
     // Default constructor.
    ReplacementTool()
    {
    }
    /// Destructor.
    virtual ~ReplacementTool()
    {
    }
    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    bool GetElementNameAndCheckModelPart(ModelPart& rModelPart, std::string& rElementName)
    {
        bool has_element = false;
        std::string reference_element_name;
        std::string element_name; 

        if(rModelPart.Elements().size() > 0)
        {
            CompareElementsAndConditionsUtility::GetRegisteredName(*rModelPart.ElementsBegin(), reference_element_name);

        #pragma omp parallel for                              
            for(int i=0; i< (int)rModelPart.Elements().size(); i++)
            {
                ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin() + i;

                
                CompareElementsAndConditionsUtility::GetRegisteredName(*it, element_name);

                KRATOS_ERROR_IF(reference_element_name != element_name)
                    << "There are more than one element types in sub model part! Replacement not possible!" << std::endl;
            }
            rElementName = reference_element_name;
            has_element = true;
        }

        return has_element;
    }   

    //-----------------------------------------------------------------------------------------------------------------
    bool GetConditionNameAndCheckModelPart(ModelPart& rModelPart, std::string& rConditionName)
    {
        bool has_condition = false;
        std::string reference_condition_name;
        std::string condition_name; 

        if(rModelPart.Conditions().size() > 0)
        {
            CompareElementsAndConditionsUtility::GetRegisteredName(*rModelPart.ConditionsBegin(), reference_condition_name);

        #pragma omp parallel for                              
            for(int i=0; i< (int)rModelPart.Conditions().size(); i++)
            {
                ModelPart::ConditionsContainerType::iterator it = rModelPart.ConditionsBegin() + i;

                
                CompareElementsAndConditionsUtility::GetRegisteredName(*it, condition_name);

                KRATOS_ERROR_IF(reference_condition_name != condition_name)
                    << "There are more than one condition types in sub model part! Replacement not possible!" << std::endl;
            }
            rConditionName = reference_condition_name;
            has_condition =  true;
        }

        return has_condition;
    }  

private:
};// class ReplacementTool

    
}  // namespace Kratos.

#endif // KRATOS_REPLACEMENT_TOOL_H_INCLUDED  defined


