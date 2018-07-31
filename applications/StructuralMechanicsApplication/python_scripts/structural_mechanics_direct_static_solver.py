from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import structural_mechanics_solver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsDirectStaticSolver(model, custom_settings)

class StructuralMechanicsDirectStaticSolver(structural_mechanics_solver.MechanicalSolver):

    def __init__(self, model, custom_settings):

        direct_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "direct_structural"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, direct_settings)
        self.scheme_settings = direct_settings["scheme_settings"]

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        # Construct the base solver.
        super(StructuralMechanicsDirectStaticSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[DirectMechanicalSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(StructuralMechanicsDirectStaticSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)
        self.print_on_rank_zero("::[DirectMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D elements for sensitivity analysis available")
        super(StructuralMechanicsDirectStaticSolver, self).PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(self.main_model_part).Execute()
        self.print_on_rank_zero("::[DirectMechanicalSolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z, self.main_model_part)
        self.print_on_rank_zero("::[DirectMechanicalSolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        if self.response_function_settings["response_type"].GetString() == "state_response":
            self.response_function = StructuralMechanicsApplication.DirectStructuralStateResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        super(StructuralMechanicsDirectStaticSolver, self).Initialize()

        self.print_on_rank_zero("::[DirectMechanicalSolver]:: ", "Finished initialization.")

    def Solve(self):
        super(StructuralMechanicsDirectStaticSolver, self).Solve()

    def SolveSolutionStep(self):
        super(StructuralMechanicsDirectStaticSolver, self).SolveSolutionStep()

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            if self.settings["compute_reactions"].GetBool():
                raise Exception("\"compute_reactions\" is not possible for sensitivity analysis models parts")
            if self.settings["move_mesh_flag"].GetBool():
                raise Exception("\"move_mesh_flag\" is not allowed for adjoint models parts")
            mechanical_solution_strategy = self._create_linear_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available for direct sensitivity analysis!\n"
            err_msg += "Available options are: \"linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_solution_scheme(self):
        return StructuralMechanicsApplication.DirectStructuralStaticScheme(self.scheme_settings, self.response_function)
