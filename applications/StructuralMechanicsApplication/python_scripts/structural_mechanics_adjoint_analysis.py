from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsAdjointAnalysis(StructuralMechanicsAnalysis):
    """
    This class is the special-script of the StructuralMechanicsApplication put in a class

    It is used to solve the adjoint problem based on the model of the primal problem.
    """
    def __init__(self, model, project_parameters):
        #solver_settings = project_parameters["solver_settings"]

        super(StructuralMechanicsAdjointAnalysis, self).__init__(model, project_parameters)

    def Initialize(self):
        #self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")
