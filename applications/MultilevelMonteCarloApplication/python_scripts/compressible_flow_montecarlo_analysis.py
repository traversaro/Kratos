from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
# from analysis_stage import AnalysisStage

# Importing derived classes
from potential_flow_analysis import PotentialFlowAnalysis

# Import pycompss
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *

# Import Monte Carlo library
import mc as mc

class MonteCarloAnalysis(PotentialFlowAnalysis):
    """Main script for Monte Carlo simulations"""
    
    def __init__(self,model,parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(model,parameters)
        # self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
    
    def _GetSimulationName(self):
        return "Monte Carlo Analysis"

    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        Mach = self.sample[0]
        a_infinity = 340 # [m/s] velocity of sound at infinity
        alpha =  self.sample[1]
        v_norm = Mach * a_infinity
        velocity = [v_norm*np.cos(alpha),v_norm*np.sin(alpha),0]
        boundary_processes = self.project_parameters["processes"]["boundary_conditions_process_list"]       
        for i in range(0,boundary_processes.size()):
            python_module = boundary_processes[i]["python_module"].GetString()
            if python_module == "apply_far_field_process":
                self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["velocity_infinity"].SetVector(velocity)


    
##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################


'''
function generating the random sample
here the sample has a normal distribution
'''
def GenerateSample():
    sample = []
    mean_Mach = 0.7
    std_deviation_Mach = 0.01
    number_samples = 1
    sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
    mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
    std_deviation_angle_attack = 0.01
    sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
    print("MACH NUMBER = ",sample[0],"ANGLE ATTACK = ",sample[1])
    if sample[0] >= 1.0 or sample[0] <= 0.0 :
        raise Exception ("Mach randomly computed > 1 or < 0")
    return sample



'''
function evaluating the QoI of the problem
'''
def EvaluateQuantityOfInterest(simulation):
    """here we evaluate the QoI of the problem: the lift coefficient"""
    Q = simulation._GetSolver().main_model_part.GetValue(KratosMultiphysics.FRICTION_COEFFICIENT)
    return Q


'''
function executing the problem
input:
        model_part_file_name : path of the model part file (still to implement how to import in efficient way in a loop where I have different model part files and different ProjectParameters files, thus for now read model part name from the ProjectParameters.json file)
        parameter_file_name  : path of the Project Parameters file
        sample               : stochastic random variable
output:
        QoI                  : Quantity of Interest
'''
@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN, returns=1)
def execution_task(model_part_file_name, parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()
    local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    sample = GenerateSample()
    simulation = MonteCarloAnalysis(model,local_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    return QoI


'''
function executing the problem for sample = 1.0
input:
        model_part_file_name  : path of the model part file
        parameter_file_name   : path of the Project Parameters file
output:
        ExactExpectedValueQoI : Quantity of Interest for sample = 1.0
'''
@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN,returns=1)
def exact_execution_task(model_part_file_name, parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()      
    local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    sample = [1.0,1.0]
    simulation = MonteCarloAnalysis(model,local_parameters, sample)
    simulation.Run() 
    ExactExpectedValueQoI = EvaluateQuantityOfInterest(simulation)
    return ExactExpectedValueQoI



'''
function computing the relative error between the Multilevel Monte Carlo expected value and the exact expected value
input :
        AveragedMeanQoI       : Multilevel Monte Carlo expected value
        ExactExpectedValueQoI : exact expected value
output :
        relative_error        : relative error
'''
@task(returns=1)
def compare_mean(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg = 'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python montecarlo_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python3 montecarlo_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "/home/kratos105b/DataDisk/MultilevelMonteCarloApplication/CompressiblePotentialFlow/ProjectParameters2.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them

    number_samples = 10
    Qlist = []

    for instance in range (0,number_samples):
        Qlist.append(execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name))

    '''Compute mean, second moment and sample variance'''
    MC_mean = 0.0
    MC_second_moment = 0.0
    for i in range (0,number_samples):
        nsam = i+1
        MC_mean, MC_second_moment, MC_variance = mc.update_onepass_M_VAR(Qlist[i], MC_mean, MC_second_moment, nsam)

    MC_mean = compss_wait_on(MC_mean)
    print("\nlist of lift coefficients computed = ",Qlist)
    print("\nMC mean = ",MC_mean)
    
