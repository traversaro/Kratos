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

# construct primal and adjoint analysis
my_analysis_model = Model()
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


