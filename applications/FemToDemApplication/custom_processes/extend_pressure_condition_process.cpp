//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#include "custom_processes/extend_pressure_condition_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include <string>

namespace Kratos {

template <SizeType TDim>
ExtendPressureConditionProcess<TDim>::ExtendPressureConditionProcess(
    ModelPart &r_model_part)
    : mr_model_part(r_model_part)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::Execute() 
{
	int maximum_condition_id, counter_of_affected_nodes = 0;
    this->GetMaximumConditionIdOnSubmodelPart(maximum_condition_id);

    for (ModelPart::ElementsContainerType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it) {
        bool condition_is_active = true;
        if ((*it)->IsDefined(ACTIVE)) {
            condition_is_active = (*it)->Is(ACTIVE);
        }
        // It's going to be removed
        if (condition_is_active == false && (*it)->GetValue(SMOOTHING) == false) {
            unsigned int local_id, counter = 0, pressure_id;
            // Loop over nodes in order to check if there's pressure on nodes
            for (IndexType i = 0; i < (*it)->GetGeometry().PointsNumber(); ++i) {
                if ((*it)->GetGeometry().GetPoint(i).GetValue(PRESSURE_ID) != 0) {
                    pressure_id = (*it)->GetGeometry().GetPoint(i).GetValue(PRESSURE_ID);
                    counter++;
                } else {
                    local_id = i;
                }
            }
            if (counter == 2) {
                KRATOS_WATCH((*it)->Id())
                this->CreateAndAddPressureConditions2(it, local_id, pressure_id, maximum_condition_id);
                counter_of_affected_nodes++;
                // We use this flag to enter once on each element
                (*it)->SetValue(SMOOTHING, true);
            } else if (counter == 3) {
                KRATOS_WATCH((*it)->Id())
                this->CreateAndAddPressureConditions3(it, pressure_id, maximum_condition_id);
                counter_of_affected_nodes++;
                // We use this flag to enter once on each element
                (*it)->SetValue(SMOOTHING, true);
            }
        }
    }
    auto& process_info = mr_model_part.GetProcessInfo();
    process_info[ITER] = counter_of_affected_nodes;

    for (IndexType i = 0; i < mNodeIdContainer.size(); ++i) {
        mr_model_part.GetNode(mNodeIdContainer[i]).SetValue(PRESSURE_ID, mNodePressureIdContainer[i]);
    }
    mNodeIdContainer.clear();
    mNodePressureIdContainer.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions2(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const unsigned int LocalId,
    const int PressureId,
	int& MaximumConditionId
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mr_model_part.GetSubModelPart(sub_model_name);

    std::vector<IndexType> condition_nodes_id(2);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);

    auto& r_geom = (*itElem)->GetGeometry();
    r_sub_model_part.AddNode(mr_model_part.pGetNode(r_geom[LocalId].Id()));
	
    // Set the flag
    // r_geom[LocalId].SetValue(PRESSURE_ID, PressureId);
    mNodeIdContainer.push_back(r_geom[LocalId].Id());
    mNodePressureIdContainer.push_back(PressureId);

    const IndexType id_1 = LocalId == 0 ? 0 : LocalId == 1 ? 1 : 2;
    const IndexType id_2 = LocalId == 0 ? 1 : LocalId == 1 ? 2 : 0;
    const IndexType id_3 = LocalId == 0 ? 2 : LocalId == 1 ? 0 : 1;

    condition_nodes_id[0] = r_geom[id_2].Id();
    condition_nodes_id[1] = r_geom[id_1].Id();
	MaximumConditionId++;
    auto& line_cond1 = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    MaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);

    condition_nodes_id[0] = r_geom[id_1].Id();
    condition_nodes_id[1] = r_geom[id_3].Id();
    MaximumConditionId++;
    auto& line_cond2 = r_sub_model_part.CreateNewCondition(
					                    "LineLoadCondition2D2N",
					                    MaximumConditionId,
					                    condition_nodes_id,
					                    p_properties, 0);
    
    // adding the conditions to the computing model part
    mr_model_part.GetSubModelPart("computing_domain").AddCondition(line_cond1);
    mr_model_part.GetSubModelPart("computing_domain").AddCondition(line_cond2);

    // We remove the condition regarding the erased edge...
    for (ModelPart::ConditionsContainerType::ptr_iterator it = mr_model_part.Conditions().ptr_begin(); it != mr_model_part.Conditions().ptr_end(); ++it) {
        // Nodes of the condition
        const IndexType Id1 = (*it)->GetGeometry()[0].Id();
        const IndexType Id2 = (*it)->GetGeometry()[1].Id();
        if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_3].Id()) || (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_3].Id())) {
			mr_model_part.RemoveCondition((*it)->Id());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions3(
    ModelPart::ElementsContainerType::ptr_iterator itElem,
    const int PressureId,
	int& MaximumConditionId
    )
{
    std::string sub_model_name;
	sub_model_name = "Normal_Load-auto-" + std::to_string(PressureId);
    auto& r_sub_model_part = mr_model_part.GetSubModelPart(sub_model_name);

    std::vector<IndexType> condition_nodes_id(2);
    ModelPart::PropertiesType::Pointer p_properties = r_sub_model_part.pGetProperties(1);
    auto& r_geom = (*itElem)->GetGeometry();

    IndexType local_id;
    for (IndexType i = 0; i < (*itElem)->GetGeometry().size(); ++i) {
        if ((*itElem)->GetGeometry()[i].GetValue(NUMBER_OF_ACTIVE_ELEMENTS) == 1) {
            local_id = i;
        }
    }
    const IndexType id_1 = local_id == 0 ? 0 : local_id == 1 ? 1 : 2;
    const IndexType id_2 = local_id == 0 ? 1 : local_id == 1 ? 2 : 0;
    const IndexType id_3 = local_id == 0 ? 2 : local_id == 1 ? 0 : 1;

    condition_nodes_id[0] = r_geom[id_2].Id();
    condition_nodes_id[1] = r_geom[id_3].Id();
	MaximumConditionId++;
    auto& line_cond = r_sub_model_part.CreateNewCondition(
                                       "LineLoadCondition2D2N",
                                       MaximumConditionId,
                                       condition_nodes_id,
                                       p_properties, 0);

    // adding the conditions to the computing model part
    mr_model_part.GetSubModelPart("computing_domain").AddCondition(line_cond);

    // We remove the condition regarding the erased edges...
    int counter = 0;
    for (ModelPart::ConditionsContainerType::ptr_iterator it = mr_model_part.Conditions().ptr_begin(); it != mr_model_part.Conditions().ptr_end(); ++it) {
        const IndexType Id1 = (*it)->GetGeometry()[0].Id();
        const IndexType Id2 = (*it)->GetGeometry()[1].Id();
        KRATOS_WATCH((*it)->Id())
        KRATOS_WATCH(Id1)
        KRATOS_WATCH(Id2)
        
        if ((Id1 == r_geom[id_2].Id() && Id2 == r_geom[id_1].Id()) || (Id2 == r_geom[id_2].Id() && Id1 == r_geom[id_1].Id())) {
			mr_model_part.RemoveCondition((*it)->Id());
            counter++;
        } else if ((Id1 == r_geom[id_1].Id() && Id2 == r_geom[id_3].Id()) || (Id2 == r_geom[id_1].Id() && Id1 == r_geom[id_3].Id())) {
			mr_model_part.RemoveCondition((*it)->Id());
            counter++;
        }
    }
    KRATOS_ERROR_IF(counter < 2) << "NO PUEDE SER..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::GetMaximumConditionIdOnSubmodelPart(
      int& MaximumConditionId
)
{
    MaximumConditionId = 0;
    for (ModelPart::ConditionIterator itCond = mr_model_part.ConditionsBegin();
        itCond != mr_model_part.ConditionsEnd();
        itCond++) {

        if (((*itCond)).Id() > MaximumConditionId) MaximumConditionId = ((*itCond)).Id();
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<3>::Execute() 
{

}

/***********************************************************************************/
/***********************************************************************************/

template class ExtendPressureConditionProcess<2>;
template class ExtendPressureConditionProcess<3>;

}  // namespace Kratos
