from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import CouplingFemDem3D

def Wait():
	input("Press Something")

# Main script of the coupled FEM-DEM Application 3D for hexahedrons
class FEMDEM3DHexahedrons_Solution(CouplingFemDem3D.FEMDEM3D_Solution):

#============================================================================================================================
	def GenerateDEM(self): # 3D version for hexahedrons
		pass
#============================================================================================================================
	def CheckForPossibleIndentations(self): # Verifies if an element has indentations between its DEM
		pass
#============================================================================================================================
