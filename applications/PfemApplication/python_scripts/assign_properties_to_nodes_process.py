from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignPropertiesToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignPropertiesToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "model_part_name": "main_domain",
             "fluid_mixture"  : false,
             "solid_mixture"  : false
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

    @classmethod
    def GetVariables(self):
        nodal_variables = []
        return nodal_variables

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        self.settings.RemoveValue("model_part_name")
        self.AssignPropertiesProcess = KratosPfem.AssignPropertiesToNodes(self.model_part, self.settings)

        for Element in self.model_part.Elements:

            density = Element.Properties.GetValue(KratosMultiphysics.DENSITY)
            bulk_modulus = Element.Properties.GetValue(KratosMultiphysics.BULK_MODULUS)
            viscosity = Element.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            gravity = Element.Properties.GetValue(KratosMultiphysics.GRAVITY)

        for Nodes in self.model_part.Nodes:

            Nodes.SetSolutionStepValue(KratosMultiphysics.DENSITY,density)
            Nodes.SetSolutionStepValue(KratosMultiphysics.BULK_MODULUS,bulk_modulus)
            Nodes.SetSolutionStepValue(KratosMultiphysics.VISCOSITY,viscosity)
            Nodes.SetSolutionStepValue(KratosMultiphysics.GRAVITY,gravity)


        self.AssignPropertiesProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.AssignPropertiesProcess.ExecuteInitializeSolutionStep()
