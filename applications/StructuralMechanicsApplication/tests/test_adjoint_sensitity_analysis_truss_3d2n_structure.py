from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
import structural_mechanics_adjoint_analysis

try:
    from KratosMultiphysics.HDF5Application import *
    has_hdf5_application = True
except ImportError:
    has_hdf5_application = False

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(this_file_dir, "adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n")

@KratosUnittest.skipUnless(has_hdf5_application,"Missing required application: HDF5Application")
class TestAdjointSensitivityAnalysisBeamStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        pass

    def test_local_stress_response(self):
        # solve primal problem
        with controlledExecutionScope(_get_test_working_dir()):
            with open("nonlinear_truss_test_parameters.json",'r') as parameter_file:
                ProjectParametersPrimal = Parameters( parameter_file.read())

            model_primal = Model()
            primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
            primal_analysis.Run()

        #Create the adjoint solver
        with controlledExecutionScope(_get_test_working_dir()):
            with open("nonlinear_truss_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
                ProjectParametersAdjoint = Parameters( parameter_file.read())

            #model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
            #model_adjoint = Model()

            adjoint_analysis = structural_mechanics_adjoint_analysis.StructuralMechanicsAdjointAnalysis(model_primal, ProjectParametersAdjoint)
            #adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
            adjoint_analysis.Run()

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        with controlledExecutionScope(_get_test_working_dir()):
            kratos_utilities.DeleteFileIfExisting("nonlinear_truss_structure.time")
            for file_name in os.listdir():
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()
