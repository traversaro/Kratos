from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AdjointReplacementProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class AdjointReplacementProcess(KM.Process):
    """Replace elements and conditions."""
    def __init__(self, main_model_part, Parameters):
        self.main_model_part = main_model_part
        
        self.structural_model_part_names = Parameters["problem_domain_sub_model_part_list"]
        self.processes_model_part_names = Parameters["processes_sub_model_part_list"]

    def Execute(self, from_primal_to_adjoint = True):
        
        structural_parts = []
        for i in range(self.structural_model_part_names.size()):
            structural_parts.append(self.main_model_part.GetSubModelPart(self.structural_model_part_names[i].GetString()))
        
        process_parts = []
        for i in range(self.processes_model_part_names.size()):
            process_parts.append(self.main_model_part.GetSubModelPart(self.processes_model_part_names[i].GetString()))
        
        has_element = False
        self.element_name = None
        has_condition = False
        self.condition_name = None

        if(from_primal_to_adjoint):
            for part in structural_parts:
                has_element = StructuralMechanicsApplication.ReplacementTool.GetElementNameAndCheckModelPart(part, self.element_name)
                if(has_element):
                    if (self.element_name is "CrLinearBeamElement3D2N"):
                        new_elem_name = "CrLinearBeamAdjointElement3D2N"
                    elif (self.element_name is "ShellThinElement3D3N"):
                        new_elem_name = "ShellThinAdjointElement3D3N"  
                    else:       
                        raise Exception("There is no corresponding adjoint element to", element_name, "!")

                    self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString("Condition") # dummy condition
                    ## Call the replace elements and conditions process
                    KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            for part in process_parts:
                has_condition = StructuralMechanicsApplication.ReplacementTool.GetConditionNameAndCheckModelPart(part, self.condition_name)
                if(has_condition):
                    if (self.condition_name is "PointLoadCondition2D1N"):
                        new_condition_name = "PointLoadAdjointCondition2D1N"
                    elif (self.element_name is "PointLoadCondition3D1N"):
                        new_condition_name = "PointLoadAdjointCondition3D1N"  
                    elif (self.element_name is "SurfaceLoadCondition3D3N"):
                        new_condition_name = "SurfaceLoadAdjointCondition3D3N"  
                    elif (self.element_name is "SurfaceLoadCondition3D4N"):
                        new_condition_name = "SurfaceLoadAdjointCondition3D4N"  
                    else:       
                        raise Exception("There is no corresponding adjoint condition to", condition_name, "!")

                    self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString("Element") # dummy element
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_condition_name)
                    ## Call the replace elements and conditions process
                    KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        else:
            for part in structural_parts:
                has_element = StructuralMechanicsApplication.ReplacementTool.GetElementNameAndCheckModelPart(part, self.element_name)
                if(has_element):
                    if (self.element_name is "CrLinearBeamAdjointElement3D2N"):
                        new_elem_name = "CrLinearBeamElement3D2N"
                    elif (self.element_name is "ShellThinAdjointElement3D3N"):
                        new_elem_name = "ShellThinElement3D3N"  
                    else:       
                        raise Exception("There is no corresponding primal element to", element_name, "!")

                    self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString("Condition") # dummy condition
                    ## Call the replace elements and conditions process
                    KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            for part in process_parts:
                has_condition = StructuralMechanicsApplication.ReplacementTool.GetConditionNameAndCheckModelPart(part, self.condition_name)
                if(has_condition):
                    if (self.condition_name is "PointLoadAdjointCondition2D1N"):
                        new_condition_name = "PointLoadCondition2D1N"
                    elif (self.element_name is "PointLoadAdjointCondition3D1N"):
                        new_condition_name = "PointLoadCondition3D1N"  
                    elif (self.element_name is "SurfaceLoadAdjointCondition3D3N"):
                        new_condition_name = "SurfaceLoadACondition3D3N"  
                    elif (self.element_name is "SurfaceLoadAdjointCondition3D4N"):
                        new_condition_name = "SurfaceLoadCondition3D4N"    
                    else:       
                        raise Exception("There is no corresponding primal condition to", condition_name, "!")

                    self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString("Element") # dummy element
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_condition_name)
                    ## Call the replace elements and conditions process
                    KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            #TODO: zu klären, darf core replacement process verändert werden (Übergabewerte an Create(..) )
            #TODO: klappt Verwendung von 'problem_domain_sub_model_part_list' und 'processes_sub_model_part_list'?  diese müssten den adjoint
            #      parameters hinzugefügt werden. 


##--------------------------------------------------------------------------------------------------------------
#def _replace_elements_and_conditions(self):
#
#    submodelpartslist = self.__generate_submodelparts_list_from_input(self.settings["model_parts_to_replace"])
#  
#    for submodelpart in submodelpartslist:
#
#        StructuralMechanicsApplication.ReplacementTool.GetNamesAndCheckModelPart(submodelpart, self.element_name, self.condition_name)
#
#
#
#   
#      
#
##--------------------------------------------------------------------------------------------------------------
#def __generate_submodelparts_list_from_input(self,param):
#    '''Parse a list of variables from input.'''
#    # At least verify that the input is a string
#    if not param.IsArray():
#        raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))
#    # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
#    return [self.model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]
#
#
##--------------------------------------------------------------------------------------------------------------
#def _get_element_num_nodes(self):
#    if self.main_model_part.NumberOfElements() != 0:
#        if sys.version_info[0] >= 3: # python3 syntax
#            element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
#        else: # python2 syntax
#            element_num_nodes = len(self.main_model_part.Elements.__iter__().next().GetNodes())
#    else:
#        element_num_nodes = 0
#    element_num_nodes = self.main_model_part.GetCommunicator().MaxAll(element_num_nodes)
#    return element_num_nodes
#def _get_condition_num_nodes(self):
#    if self.main_model_part.NumberOfConditions() != 0:
#        if sys.version_info[0] >= 3: # python3 syntax
#            condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
#        else: # python2 syntax
#            condition_num_nodes = len(self.main_model_part.Conditions.__iter__().next().GetNodes())
#    else:
#        condition_num_nodes = 0
#    condition_num_nodes = self.main_model_part.GetCommunicator().MaxAll(condition_num_nodes)
#    return condition_num_nodes