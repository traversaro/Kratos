from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import os
#os.environ['OMP_NUM_THREADS'] = "1"
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis

# read parameter files
with open("Shell_Q3_Thin_nonlinear_static_test_parameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())
with open("Shell_Q3_Thin_nonlinear_static_adjoint_test_parameters.json",'r') as parameter_file:
    ProjectParametersAdjoint = Parameters( parameter_file.read())

my_analysis_model = Model()

# Add model part in order to unsure that primal and adjoint solver using the same model part
solver_settings = ProjectParametersPrimal["solver_settings"]
model_part_name = solver_settings ["model_part_name"].GetString()
if not my_analysis_model.HasModelPart(model_part_name):
    model_part = ModelPart(model_part_name)
    domain_size = solver_settings["domain_size"].GetInt()
    if domain_size < 0:
        raise Exception('Please specify a "domain_size" >= 0!')
    model_part.ProcessInfo.SetValue(DOMAIN_SIZE, domain_size)
    my_analysis_model.AddModelPart(model_part)

# construct primal and adjoint analysis
primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(my_analysis_model, ProjectParametersPrimal)
adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(my_analysis_model, ProjectParametersAdjoint)


# run analysis
primal_analysis.Run()


print("")
print("------------------------------------")
print("Finished primal problem!")
print("------------------------------------")
print("")

adjoint_analysis.Run()

print("")
print("------------------------------------")
print("Finished adjoint problem!")
print("------------------------------------")
print("")


    #if not model.HasModelPart(model_part_name):
    #    model_part = ModelPart(model_part_name)
    #    domain_size = solver_settings["domain_size"].GetInt()
    #    if domain_size < 0:
    #        raise Exception('Please specify a "domain_size" >= 0!')
    #    model_part.ProcessInfo.SetValue(DOMAIN_SIZE, domain_size)
    #    model.AddModelPart(model_part)
    #else:
    #    model_part = model.GetModelPart(model_part_name)