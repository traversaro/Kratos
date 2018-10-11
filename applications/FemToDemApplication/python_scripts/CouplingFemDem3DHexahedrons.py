from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import CouplingFemDem3D
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.FemToDemApplication   as KratosFemDem

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

	def CheckInactiveNodes(self):
		pass
#============================================================================================================================

	def InitializeSolutionStep(self):

        # modified for the remeshing
		self.FEM_Solution.delta_time = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
		self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
		self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
		self.FEM_Solution.step = self.FEM_Solution.step + 1
		self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

		if self.DoRemeshing:
			is_remeshing = self.CheckIfHasRemeshed()
			
			if is_remeshing:
				# Extrapolate the VonMises normalized stress to nodes (remeshing)
				KratosFemDem.StressToNodesProcess(self.FEM_Solution.main_model_part, 2).Execute()

			# Perform remeshing
			self.RemeshingProcessMMG.ExecuteInitializeSolutionStep()

			self.nodal_neighbour_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.FEM_Solution.main_model_part, 4, 5)
			self.nodal_neighbour_finder.Execute()

			if is_remeshing:

				# Initialize the "flag" IS_DEM in all the nodes
				KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, self.FEM_Solution.main_model_part.Nodes)
				# Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
				KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, self.FEM_Solution.main_model_part.Nodes)
				# Initialize the "flag" RADIUS in all the nodes
				KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.RADIUS, False, self.FEM_Solution.main_model_part.Nodes)

				# Remove DEMS from previous mesh
				self.SpheresModelPart.Elements.clear()
				self.SpheresModelPart.Nodes.clear()

				self.InitializeMMGvariables()
				self.FEM_Solution.model_processes = self.FEM_Solution.AddProcesses()
				self.FEM_Solution.model_processes.ExecuteInitialize()
				self.FEM_Solution.model_processes.ExecuteBeforeSolutionLoop()
				self.FEM_Solution.model_processes.ExecuteInitializeSolutionStep()

		self.FEM_Solution.InitializeSolutionStep()

		# Create initial skin of DEM's
		self.create_initial_dem_skin = False  # Hard Coded TODO
		if self.create_initial_dem_skin and self.FEM_Solution.step == 1:
			self.CreateInitialSkinDEM()

		# Create the DEM after the remeshing
		if self.DoRemeshing and is_remeshing:
			self.GenerateDemAfterRemeshing()

#============================================================================================================================
