from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/cpuigbo/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");
import time as timer
# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.PfemApplication               as KratosPfem
#import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid
import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import MainPfem
import MainSolid

class Solution(MainPfem.PfemSolution):

    def __init__(self, Model, file_parameters = "ProjectParameters.json", file_name = None):

        self.pp = self.ProblemParameters()
        super(Solution,self).__init__(Model, file_parameters, file_name)

        self.main_model_part = self.model.GetMainModelPart()

        self.fluid_model_part = self.main_model_part
        self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY,self.ProjectParameters["problem_data"]["gravity_vector"].GetVector())

        self.pp.nodal_results = []
        output_settings = self.ProjectParameters["output_configuration"]["result_file_configuration"]
        for i in range(output_settings["nodal_results"].size()):
            self.pp.nodal_results.append(output_settings["nodal_results"][i].GetString())

        self.pp.gauss_points_results = []
        for i in range(output_settings["gauss_point_results"].size()):
            self.pp.nodal_results.append(output_settings["gauss_point_results"][i].GetString())

        self.pp.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()
        self.pp.Start_time = self.ProjectParameters["time_settings"]["start_time"].GetDouble()
        self.pp.VolumeOutput = True

        if output_settings["gidpost_flags"]["GiDPostMode"].GetString() == "GiD_PostBinary":
            self.pp.GiDPostMode  = "Binary"
        else:
            self.pp.GiDPostMode  = "Ascii"

        if output_settings["gidpost_flags"]["MultiFileFlag"].GetString() == "MultipleFiles":
            self.pp.GiDMultiFileFlag  = "Multiples"
        else:
            self.pp.GiDMultiFileFlag  = "Single"

        if output_settings["gidpost_flags"]["WriteDeformedMeshFlag"].GetString() == "WriteDeformed":
            self.pp.GiDWriteMeshFlag  = True
        else:
            self.pp.GiDWriteMeshFlag  = False

        if output_settings["gidpost_flags"]["WriteConditionsFlag"].GetString() == "WriteConditions":
            self.pp.GiDWriteConditionsFlag  = True
        else:
            self.pp.GiDWriteConditionsFlag  = False


    def _get_solver(self):
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.AddFluidVariablesBySwimmingDEMAlgorithm()
        return (solver_module.CreateSolver(self.ProjectParameters["solver_settings"]["Parameters"], self.model.GetModel()))


    def CalculateNodalArea(self):

        KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part,self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]).Execute()

    def AddFluidVariablesBySwimmingDEMAlgorithm(self):
        self.vars_man.AddNodalVariables(self.main_model_part, self.pp.fluid_vars)

    def GetDeltaTimeFromParameters(self):
        return self.ProjectParameters["time_settings"]["time_step"].GetDouble()

    class ProblemParameters:
        def __init__(self):
            pass


if __name__ == "__main__":
    Solution().Run()
