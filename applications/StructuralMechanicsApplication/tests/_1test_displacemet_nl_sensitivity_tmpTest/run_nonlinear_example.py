from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_mechanics_analysis
import structural_mechanics_analysis_nonlinear_sensitivity
import numpy as np
import matplotlib.pyplot as plt

def compute_EF_disp_curvature(model_part, load_factor):
    for node in model_part.Nodes:
        disp_1_x = node.GetSolutionStepValue(DISPLACEMENT_X, 3)
        disp_2_x = node.GetSolutionStepValue(DISPLACEMENT_X, 2)
        disp_3_x = node.GetSolutionStepValue(DISPLACEMENT_X, 1)
        disp_1_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 3)
        disp_2_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 2)
        disp_3_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 1)
        disp_1_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 3)
        disp_2_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 2)
        disp_3_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 1)

        disp_x = np.array([disp_1_x, disp_2_x, disp_3_x])
        px = np.polyfit(load_factor, disp_x, 2)
        disp_y = np.array([disp_1_y, disp_2_y, disp_3_y])
        py = np.polyfit(load_factor, disp_y, 2)
        disp_z = np.array([disp_1_z, disp_2_z, disp_3_z])
        pz = np.polyfit(load_factor, disp_z, 2)

        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_X, 2 * px[0])
        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_Y, 2 * py[0])
        node.SetValue(DISPLACEMENT_NL_SENSITIVITY_Z, 2 * pz[0])
        #print("curvature x = ", 2 * px[0])
        #print("curvature y = ", 2 * py[0])
        #print("curvature z = ", 2 * pz[0])
        #print("")

def compute_first_and_second_order_nl_sensitivity_factors(model_part, load_factor):
    lambda_0 = load_factor[0]
    lambda_1 = load_factor[1]
    lambda_2 = load_factor[2]
    f_1 = lambda_1 / lambda_0
    #f_2 = lambda_2 / lambda_0
    delta_10 = lambda_1 - lambda_0
    delta_20 = lambda_2 - lambda_0
    for node in model_part.Nodes:
        disp_1_x = node.GetSolutionStepValue(DISPLACEMENT_X, 2)
        disp_2_x = node.GetSolutionStepValue(DISPLACEMENT_X, 1)
        disp_3_x = node.GetSolutionStepValue(DISPLACEMENT_X, 0)
        disp_1_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 2)
        disp_2_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 1)
        disp_3_y = node.GetSolutionStepValue(DISPLACEMENT_Y, 0)
        disp_1_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 2)
        disp_2_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 1)
        disp_3_z = node.GetSolutionStepValue(DISPLACEMENT_Z, 0)
        print("Displacement = ", disp_1_x)
        print("Displacement = ", disp_2_x)
        print("Displacement = ", disp_3_x)

        if abs(disp_1_x) > 1e-8:
            sensitivity_first_order_1_x = disp_2_x / ( disp_1_x * f_1 )
            #sensitivity_first_order_2_x = disp_3_x / ( disp_1_x * f_2 )
            #sensitivity_second_order_x = sensitivity_first_order_2_x / sensitivity_first_order_1_x
            slope_10_x = (disp_2_x-disp_1_x)/delta_10
            slope_20_x = (disp_3_x-disp_1_x)/delta_20
            sensitivity_second_order_x = slope_20_x / slope_10_x
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_X, sensitivity_first_order_1_x)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_X, sensitivity_second_order_x)
        else:
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_X, 0.0)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_X, 0.0)

        if abs(disp_1_y) > 1e-8:
            sensitivity_first_order_1_y = disp_2_y / ( disp_1_y * f_1 )
            #sensitivity_first_order_2_y = disp_3_y / ( disp_1_y * f_2 )
            #sensitivity_second_order_y = sensitivity_first_order_2_y / sensitivity_first_order_1_y
            slope_10_y = (disp_2_y-disp_1_y)/delta_10
            slope_20_y = (disp_3_y-disp_1_y)/delta_20
            sensitivity_second_order_y = slope_20_y / slope_10_y
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_Y, sensitivity_first_order_1_y)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_Y, sensitivity_second_order_y)
        else:
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_Y, 0.0)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_Y, 0.0)

        if abs(disp_1_z) > 1e-8:
            sensitivity_first_order_1_z = disp_2_z / ( disp_1_z * f_1 )
            #sensitivity_first_order_2_z = disp_3_z / ( disp_1_z * f_2 )
            #sensitivity_second_order_z = sensitivity_first_order_2_z / sensitivity_first_order_1_z
            slope_10_z = (disp_2_z-disp_1_z)/delta_10
            slope_20_z = (disp_3_z-disp_1_z)/delta_20
            sensitivity_second_order_z = slope_20_z / slope_10_z
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_Z, sensitivity_first_order_1_z)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_Z, sensitivity_second_order_z)
        else:
            node.SetValue(NL_SENSITIVITY_FIRST_ORDER_Z, 0.0)
            node.SetValue(NL_SENSITIVITY_SECOND_ORDER_Z, 0.0)

# begin first analysis ****************************************************************************************************

with open("PrimalParameters.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()
model_truss_primal = Model()

primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_truss_primal, ProjectParametersPrimal)

primal_analysis.Initialize()

primal_analysis.RunSolutionLoop()

primal_analysis.Finalize()

# end first analysis ****************************************************************************************************

print("")

# begin second analysis ****************************************************************************************************
model_part = model_truss_primal.GetModelPart(model_part_name_primal)
#time = model_part.ProcessInfo[TIME]
#print("time after primal analysis = ", time , " ****************************************************")

with open("PrimalParameters_curvature.json",'r') as parameter_file:
    ProjectParametersPrimal = Parameters( parameter_file.read())

model_part_name_primal = ProjectParametersPrimal["problem_data"]["model_part_name"].GetString()

#model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True

curvature_analysis = structural_mechanics_analysis_nonlinear_sensitivity.StructuralMechanicsAnalysisNLSensitivity(model_truss_primal, ProjectParametersPrimal)

for node in model_part.Nodes:
    disp_3_x = node.GetSolutionStepValue(DISPLACEMENT_X, 0)
    print("Displacement = ", disp_3_x, " ***************************************************************************************")

curvature_analysis.Initialize()

for node in model_part.Nodes:
    disp_3_x = node.GetSolutionStepValue(DISPLACEMENT_X, 2)
    print("Displacement = ", disp_3_x, " ***************************************************************************************")

model_part.SetBufferSize(5)

curvature_analysis.RunSolutionLoop()

# compute sensitivity measures
load_factors = np.array([1.0, 1.175, 1.35])
compute_EF_disp_curvature(model_part, load_factors)
compute_first_and_second_order_nl_sensitivity_factors(model_part, load_factors)

curvature_analysis.Finalize()
model_part = model_truss_primal.GetModelPart(model_part_name_primal)
#time_2 = model_part.ProcessInfo[TIME]
#print("time after 3-step analysis = ", time_2 , " ****************************************************")

# end second analysis ****************************************************************************************************

for node in model_part.Nodes:
    print("curvature x = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_X))
    print("curvature y = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Y))
    print("curvature z = ", node.GetValue(DISPLACEMENT_NL_SENSITIVITY_Z))
    print("")
    print("first order sen x = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_X))
    print("first order sen y = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_Y))
    print("first order sen z = ", node.GetValue(NL_SENSITIVITY_FIRST_ORDER_Z))
    print("")
    print("second order sen x = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_X))
    print("second order sen y = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_Y))
    print("second order sen z = ", node.GetValue(NL_SENSITIVITY_SECOND_ORDER_Z))
    print("")
    print("********************************************")
    print("")



#print("")
#print("polynomial coeff = ", z)
#print("curvature = ", 2 * z[0])
#print("")
#print(disp_1)
#print(disp_2)
#print(disp_3)
#print("Finished polynomial interpolation!")

#p = np.poly1d(z)
#print(np.poly1d(p))
#xp = np.linspace(0, 1.0, 100)
#_ = plt.plot(x, y, '.', xp, p(xp), '-')
#plt.ylim(0,0.2)
#plt.show()



