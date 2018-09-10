from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import os
#os.environ['OMP_NUM_THREADS'] = "1"
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis
import structural_mechanics_adjoint_analysis

with open("Shell_Q3_Thin_nonlinear_static_test_parameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()
my_analysis_model = Model()

primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(my_analysis_model, ProjectParametersPrimal)
primal_analysis.Run()

print("Finished primal problem!")

with open("Shell_Q3_Thin_nonlinear_static_adjoint_test_parameters.json",'r') as parameter_file:
    ProjectParametersAdjoint = Parameters( parameter_file.read())

#model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
#model_adjoint = Model()

adjoint_analysis = structural_mechanics_adjoint_analysis.StructuralMechanicsAdjointAnalysis(my_analysis_model, ProjectParametersAdjoint)
adjoint_analysis.Run()


#print("Finished primal problem!")
#
#with open("AdjointParameters.json",'r') as parameter_file:
#    ProjectParametersAdjoint = Parameters( parameter_file.read())
#
#model_part_name_adjoint = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
#model_truss_adjoint = Model()
#
#adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_truss_adjoint, ProjectParametersAdjoint)
#adjoint_analysis.Run()


