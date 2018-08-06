from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis

with open("primal_parameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_primal = Model()

primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
primal_analysis.Run()

print("Finsished primal problem! Start to solve adjoint problem .. ******************************************************************************")

with open("adjoint_parameters.json",'r') as parameter_file:
    ProjectParametersAdjoint = Parameters( parameter_file.read())

model_adjoint = Model()

adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
adjoint_analysis.Run()

model_part_name_adjoint = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
adj_mp = model_adjoint.GetModelPart(model_part_name_adjoint)

influnce_fct_1 = adj_mp.GetNode(16).GetSolutionStepValue(ADJOINT_DISPLACEMENT_Y)
primal_dis = adj_mp.GetNode(16).GetSolutionStepValue(DISPLACEMENT_Y)
print("")
print("primal state variable = ", primal_dis)
print("adjoint variable = ", influnce_fct_1)
print("")

print("shape sensitivities = ")
for node in adj_mp.Nodes:
    shape_sens = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
    print("Node id #",node.Id,": ",shape_sens)
print("")    
#print("stiffness sensitivities = ")
#for element in adj_mp.Elements:
#    E_sens = element.GetValue(YOUNG_MODULUS_SENSITIVITY)
#    print("Element id #",element.Id,": ", E_sens) 
