from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
import structural_mechanics_static_solver


def CreateSolver(model, custom_settings):
    return StaticMechanicalPrimalSolver(model, custom_settings)


class StaticMechanicalPrimalSolver(structural_mechanics_static_solver.StaticMechanicalSolver):

    def __init__(self, model, custom_settings):
        super(StaticMechanicalPrimalSolver, self).__init__(model, custom_settings)

    def AddVariables(self):
        super(StaticMechanicalPrimalSolver, self).AddVariables()
        # Add variables of the upcoming adjoint sensitivity analysis
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)

        self.print_on_rank_zero("::[StaticMechanicalPrimalSolver]:: ", "Variables ADDED")


