//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
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

    for (ModelPart::ElementsContainerType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it) {

		bool condition_is_active = true;
		if ((*it)->IsDefined(ACTIVE)) {
			condition_is_active = (*it)->Is(ACTIVE);
		}

		// It's going to be removed
		if (condition_is_active == false) {
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
				this->CreateAndAddPressureConditions(it, local_id, pressure_id);
			}
		}
	}
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ExtendPressureConditionProcess<2>::CreateAndAddPressureConditions(
	ModelPart::ElementsContainerType::ptr_iterator it,
	unsigned int LocalId,
	int PressureId
) 
{
	int maximum_cond_id;
	std::string sub_model_name;
	this->GetMaximumConditionIdOnSubmodelPart(PressureId, maximum_cond_id, sub_model_name);

	int node_1_id, node_2_id;
	std::vector<IndexType> condition1_nodes_id, condition2_nodes_id;
	ModelPart::PropertiesType::Pointer p_properties = mr_model_part.GetSubModelPart(sub_model_name).pGetProperties(1);

	if (LocalId == 0) {
		// Set the flag
		(*it)->GetGeometry().GetPoint(0).SetValue(PRESSURE_ID, PressureId);
		const int id = (*it)->GetGeometry().GetPoint(0).Id();
		mr_model_part.GetSubModelPart(sub_model_name).AddNode(mr_model_part.pGetNode(id));
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(0).Id());
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(1).Id());
		ModelPart::ConditionType::Pointer condition1 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 1,
																	  condition1_nodes_id,
																	  p_properties);

		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(0).Id());
		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(2).Id());
		ModelPart::ConditionType::Pointer condition2 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 2,
																	  condition2_nodes_id,																	  
																	  p_properties);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition1);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition2);
		mr_model_part.AddCondition(condition1);
		mr_model_part.AddCondition(condition2);
	} else if (LocalId == 1) {
		// Set the flag
		(*it)->GetGeometry().GetPoint(1).SetValue(PRESSURE_ID, PressureId);
		const int id = (*it)->GetGeometry().GetPoint(1).Id();
		mr_model_part.GetSubModelPart(sub_model_name).AddNode(mr_model_part.pGetNode(id));
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(0).Id());
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(1).Id());
		ModelPart::ConditionType::Pointer condition1 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 1,
																	  condition1_nodes_id,
																	  p_properties);

		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(1).Id());
		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(2).Id());
		ModelPart::ConditionType::Pointer condition2 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 2,
																	  condition2_nodes_id,																	  
																	  p_properties);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition1);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition2);
		mr_model_part.AddCondition(condition1);
		mr_model_part.AddCondition(condition2);
	} else if (LocalId == 2) {
		// Set the flag
		(*it)->GetGeometry().GetPoint(2).SetValue(PRESSURE_ID, PressureId);
		const int id = (*it)->GetGeometry().GetPoint(2).Id();
		mr_model_part.GetSubModelPart(sub_model_name).AddNode(mr_model_part.pGetNode(id));
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(0).Id());
		condition1_nodes_id.push_back((*it)->GetGeometry().GetPoint(2).Id());
		ModelPart::ConditionType::Pointer condition1 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 1,
																	  condition1_nodes_id,
																	  p_properties);

		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(1).Id());
		condition2_nodes_id.push_back((*it)->GetGeometry().GetPoint(2).Id());
		ModelPart::ConditionType::Pointer condition2 = mr_model_part.CreateNewCondition(
																	 "LineLoadCondition2D2N",
																	  maximum_cond_id + 2,
																	  condition2_nodes_id,																	  
																	  p_properties);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition1);
		mr_model_part.GetSubModelPart(sub_model_name).AddCondition(condition2);
		mr_model_part.AddCondition(condition1);
		mr_model_part.AddCondition(condition2);
	}
}
/***********************************************************************************/
/***********************************************************************************/
template <>
void ExtendPressureConditionProcess<2>::GetMaximumConditionIdOnSubmodelPart(
	const int PressureId,
	int& MaximumConditionId,
	std::string& SubModelName
)
{
	MaximumConditionId = 0;
	SubModelName = "Normal_Load-auto-" + std::to_string(PressureId);
	// for (ModelPart::ConditionIterator itCond = mr_model_part.GetSubModelPart(SubModelName).ConditionsBegin();
	// 	itCond != mr_model_part.GetSubModelPart(SubModelName).ConditionsEnd();
	// 		itCond++) {

	// 	if (((*itCond)).Id() > MaximumConditionId) MaximumConditionId = ((*itCond)).Id();
	// }
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