from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignFlagsProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignFlagsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help": "This process assigns a flag to entities",
             "model_part_name": "MODEL_PART_NAME",
             "entity_type_list": [],
             "flags_list": []
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        flags_list = []
        for j in range(0, self.settings["flags_list"].size() ):
            flags_list.append(KratosMultiphysics.KratosGlobals.GetVariable(self.settings["flags_list"][j].GetString()))

        for i in range(0, self.settings["entity_type_list"].size() ):
            entity_type = self.settings["entity_type_list"][i].GetString()
            KratosSolid.AssignFlagsToEntitiesProcess(self.model_part,entity_type,flags_list).Execute()

