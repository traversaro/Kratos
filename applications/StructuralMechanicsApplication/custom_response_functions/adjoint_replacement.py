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
        
        has_element = True
        self.element_name = None
        has_condition = True
        self.condition_name = None

        self.settings = KM.Parameters("{}")

        if(from_primal_to_adjoint):
            for part in structural_parts:
                self.element_name  = SMA.ReplacementTool().GetElementNameAndCheckModelPart(part)
        
                if(str(self.element_name) != ""):
                    if (str(self.element_name) == "CrLinearBeamElement3D2N"):
                        new_elem_name = "CrLinearBeamAdjointElement3D2N"
                    elif (str(self.element_name) == "ShellThinElement3D3N"):
                        new_elem_name = "ShellThinAdjointElement3D3N"  
                    else:       
                        raise Exception("There is no corresponding adjoint element to", str(self.element_name), "!")

                    self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString("PointLoadCondition3D1N") # dummy condition
                    ## Call the replace elements and conditions process
                    #KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
                    self.settings = KM.Parameters("{}")

            for part in process_parts:
                self.condition_name = SMA.ReplacementTool().GetConditionNameAndCheckModelPart(part)
            
                if(str(self.condition_name) != ""):
                    if (str(self.condition_name) == "PointLoadCondition2D1N"):
                        new_condition_name = "PointLoadAdjointCondition2D1N"
                    elif (str(self.condition_name) == "PointLoadCondition3D1N"):
                        new_condition_name = "PointLoadAdjointCondition3D1N"  
                    elif (str(self.condition_name) == "SurfaceLoadCondition3D3N"):
                        new_condition_name = "SurfaceLoadAdjointCondition3D3N"  
                    elif (str(self.condition_name) == "SurfaceLoadCondition3D4N"):
                        new_condition_name = "SurfaceLoadAdjointCondition3D4N"  
                    else:       
                        raise Exception("There is no corresponding adjoint condition to", str(self.condition_name), "!")    

                    self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString("CrLinearBeamElement3D2N") # dummy element
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_condition_name)
                    ## Call the replace elements and conditions process
                    #KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
                    self.settings = KM.Parameters("{}")  

        else:
            for part in structural_parts:
                has_element = SMA.ReplacementTool().GetElementNameAndCheckModelPart(part,self.test_string)
                if(str(self.element_name) == ""):
                    if (str(self.element_name) == "CrLinearBeamAdjointElement3D2N"):
                        new_elem_name = "CrLinearBeamElement3D2N"
                    elif (str(self.element_name) == "ShellThinAdjointElement3D3N"):
                        new_elem_name = "ShellThinElement3D3N"  
                    else:       
                        raise Exception("There is no corresponding primal element to", str(self.element_name), "!")

                    self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString("PointLoadCondition3D1N") # dummy condition
                    ## Call the replace elements and conditions process
                    #KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
                    self.settings = KM.Parameters("{}")

            for part in process_parts:
                has_condition = SMA.ReplacementTool().GetConditionNameAndCheckModelPart(part)
                if(str(self.condition_name) == ""):
                    if (str(self.condition_name) ==  "PointLoadAdjointCondition2D1N"):
                        new_condition_name = "PointLoadCondition2D1N"
                    elif (str(self.condition_name) == "PointLoadAdjointCondition3D1N"):
                        new_condition_name = "PointLoadCondition3D1N"  
                    elif (str(self.condition_name) ==  "SurfaceLoadAdjointCondition3D3N"):
                        new_condition_name = "SurfaceLoadACondition3D3N"  
                    elif (str(self.condition_name) == "SurfaceLoadAdjointCondition3D4N"):
                        new_condition_name = "SurfaceLoadCondition3D4N"    
                    else:       
                        raise Exception("There is no corresponding primal condition to", str(self.condition_name), "!")

                    self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
                    self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString("CrLinearBeamElement3D2N") # dummy element
                    self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_condition_name)
                    ## Call the replace elements and conditions process
                    #KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
                    self.settings = KM.Parameters("{}")