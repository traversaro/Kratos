from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication as KratosEmpire
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

if "EMPIRE_API_LIBSO_ON_MACHINE" in os.environ:
    have_empire = True
else:
    have_empire = False

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

@UnitTest.skipUnless(have_empire,"EMPIRE is not available")
class TestEmpireWrapper(KratosUnittest.TestCase):
    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("Structure_EigenResults_0.post.msh")
        kratos_utils.DeleteFileIfExisting("Structure_EigenResults_0.post.res") # usually this is deleted by the check process but not if it fails


    def test_ArrayBasic(self):
        pass

    def test_FieldBasic(self):
        pass

    def test_FieldCoSim(self):
        pass

    def test_FieldCosimIterative(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
