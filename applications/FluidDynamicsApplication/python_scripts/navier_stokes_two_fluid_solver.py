from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidMonolithicSolver(model, custom_settings)

class NavierStokesTwoFluidMonolithicSolver(FluidSolver):

    def _ValidateSettings(self, settings):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluid_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "convection_settings":{
                "linear_solver_settings": {
                    "solver_type": "AMGCL",
                    "max_iteration": 200,
                    "gmres_krylov_space_dimension": 50,
                    "tolerance": 1e-6,
                    "provide_coordinates": false,
                    "smoother_type": "spai0",
                    "krylov_type": "bicgstabl",
                    "coarsening_type": "aggregation",
                    "scaling": true,
                    "verbosity": 0
                    },
                "max_cfl" : 10.0,
                "max_substeps":25,
                "cross_wind_stabilization_factor" : 0.0
            },
            "redistance_settings":{
                "linear_solver_settings": {
                    "solver_type": "AMGCL",
                    "max_iteration": 100,
                    "tolerance": 1e-6,
                    "provide_coordinates": false,
                    "smoother_type": "spai0",
                    "krylov_type": "lgmres",
                    "coarsening_type": "aggregation",
                    "scaling": true,
                    "verbosity": 1
                    },
                "max_iterations": 4
            },
            "maximum_iterations": 7,
            "dynamic_tau": 1.0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "AMGCL"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0
            },
            "move_mesh_flag": false,
            "is_slip": false,
            "slip_length": 1e+8,
            "penalty_coefficient": 10.0
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        super(NavierStokesTwoFluidMonolithicSolver,self).__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.min_buffer_size = 3

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidMonolithicSolver", "Construction of NavierStokesTwoFluidMonolithicSolver finished.")


    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidMonolithicSolver", "Fluid solver variables added correctly.")


    def ImportModelPart(self):
        super(NavierStokesTwoFluidMonolithicSolver, self).ImportModelPart()

    def PrepareModelPart(self):
        super(NavierStokesTwoFluidMonolithicSolver, self).PrepareModelPart()
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets DENSITY, DYNAMIC_VISCOSITY
            self._set_physical_properties()
            ## Sets the constitutive law
            self._set_constitutive_law()
            ## Setting the nodal distance
            self._set_distance_function()

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
                                                                            self.settings["time_order"].GetInt())

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        
        ## Construct Level Set Process
        import linear_solver_factory
        convection_settings = self.settings["convection_settings"]
        convection_linear_solver = linear_solver_factory.ConstructSolver(convection_settings["linear_solver_settings"])
        self.levelset_convector = KratosMultiphysics.LevelSetConvectionProcess3D(KratosMultiphysics.DISTANCE,
                                                                                 self.main_model_part,
                                                                                 convection_linear_solver,
                                                                                 convection_settings["max_cfl"].GetDouble(),
                                                                                 convection_settings["cross_wind_stabilization_factor"].GetDouble(),
                                                                                 convection_settings["max_substeps"].GetInt())


        ## Construct Redistance Process
        redistance_settings = self.settings["redistance_settings"]
        redistance_linear_solver = linear_solver_factory.ConstructSolver(redistance_settings["linear_solver_settings"])
        self.redistance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(self.main_model_part,
                                                                                             redistance_linear_solver,
                                                                                             redistance_settings["max_iterations"].GetInt())

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidMonolithicSolver", "Solver initialization finished.")


    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.bdf_process).Execute()
            self.levelset_convector.Execute()
            self.redistance_process.Execute()
            self._set_physical_properties()
            (self.solver).InitializeSolutionStep()


    def _set_physical_properties(self):
        # Transfer density and (dynamic) viscostity to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            air_rho = el.Properties.GetValue(KratosMultiphysics.DENSITY_AIR)
            if air_rho <= 0.0:
                raise Exception("DENSITY_AIR set to {0} in Properties {1}, positive number expected.".format(air_rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            break
        
        for node in self.main_model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho) 
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, air_rho) 
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity*air_rho/rho)
        # # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        # KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        # KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)


    def _set_constitutive_law(self):
        ## Construct the constitutive law needed for the embedded element
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.NewtonianTwoFluid3DLaw()
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            raise Exception("NewtonianTwoFluid2DLaw not implemented")

    def _set_distance_function(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            import read_distance_from_file
            DistanceUtility = read_distance_from_file.DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver","Distance function taken from the .mdpa input file.")
