 
## Adjoint Sensitivity Analysis

### General remarks:

This feature of the Structural Mechanics Applicattion provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy, stresses) with respect to different kinds of design variables (e.g. nodal coordinates, material or cross-sectional properties) with the adjoint approach. Therefore for each response function an adjoint problem has to be solved, which is the solution of a simple linear static problem. The sensitivies are than computed in a post-processing step. The implemented sensitivity analysis uses a so called semi-analytic approach which means that the derivatives at element level are computed by finite differences.

*Please note:* 
- This feature currently only works for linear problems.
- This feature makes use of the HDF5Application

### Features:  

- Availible response utilities (response functions):
     * Base class of structural response functions
     * Strain energy
     * Displacement or rotation of a node 
     * Stress resultant of a single element
  
- Schemes:
	* Scheme to solve the adjoint problem
	* Eigen solver scheme

- Processes:
    * replacement process (replaces all elements and conditions of a model with its adjoint equivalents and vice versa)

- A set of adjoint *Neumann* conditions:
     * Point loads (derived from PointLoadCondition)
     * Surface load (derived from SurfaceLoadCondition3D)
   
- Structural adjoint elements:
    * Uni-dimensional elements :
       	* Linear 3D beam element (derived from CrBeamElementLinear3D2N)
    * Two-dimensional elements :
        * Thin triangular shell (derived from ShellThinElement3D3N)

*Please note:* The adjoint elements and conditions are derived from elements/conditions of the Structural Mechanics Application and can not be seen independently from them. Rather they have to be traced as additions of their parents in order to use the parents in the context of adjoint sensitivity analysis. So basic tasks like the computation of the stiffness matrix are not overwritten. The main task of the adjoint elements/conditions is the derive different quantities with respect to the design variable or state (e.g. the right hand side or post-processing results like stresses).

### Usage: 
In order to perform a sensitivity analysis for one response function with the adjoint approach, the solutuions of two linear static problems are necessary: The primal and the adjoint problem. 

*Please note:* For the solution of the two problems differnt kind of variables are used in order to store them. For the primal prblem the usual variables ```DISPLACEMENT``` and ```ROTATION``` and for the adjoint problem ```ADJOINT_DISPLACEMENT``` and ```ADJOINT_ROTATION```

#### Definition of the Primal Problem
The primal problem can be defind by the regular input files which are needed for an usual linear static analysis. As only difference the output process of the HDF5Application has to be added to the ```list_other_processes``` in the project parameters:

```python
    "list_other_processes" :[{
        "kratos_module" : "KratosMultiphysics.HDF5Application",
        "python_module" : "single_mesh_primal_output_process",
        "help"          : "",
        "process_name"  : "",
        "Parameters" : {
            "model_part_name" : "Structure",
            "file_settings" : {
                "file_access_mode" : "truncate"
            },
            "model_part_output_settings" : {
                "prefix" : "/ModelData"
            },
            "nodal_results_settings" : {
                "list_of_variables": ["DISPLACEMENT", "ROTATION"]
            }
        }
    } 
```

#### Definition of the Adjoint Problem
In order to define the adjoint problem an additional *.json-file for the adjoint project parameters is necessary. This input file is in principle very similar to the respective file of the primal analysis. In comparsion to a regular file for a linear static analysis three points have to be modified:
- ```solver_settings``` by using the ```adjoint_structural_solver``` as ```solver_type``` and by the definion of the ```response_function_settings```
- The input process of the HDF5Application has to be added to the ```list_other_processes``` in order to read the primal solution
- When defining *Dirichlet* conditions in the ```constraints_process_list``` the ```variable_name``` has to be modified to ```ADJOINT_DISPLACEMENT``` respective ```ADJOINT_ROTATION```

For example the ```solver_settings``` can be look like this (Hints for the ```response_function_settings``` are given below):

```python
    "solver_settings"                  : {
        "solver_type"                  : "adjoint_structural_solver",
        "scheme_settings" : {
            "scheme_type"              : "structural"
            },
        "response_function_settings" : {
                "response_type"     : "adjoint_nodal_displacement",
                "use_kratos"        : true,
                "gradient_mode"     : "semi_analytic",
                "sensitivity_model_part_name" : "Parts_Beam",
                "nodal_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
                "element_sensitivity_variables"  : ["I22"],
                "condition_sensitivity_variables"  : [],
                "step_size"         : 1e-6,
                "traced_node"       : 6,
                "traced_dof"        : "DISPLACEMENT_Z"

            },
        "echo_level"                   : 0,
        "problem_domain_sub_model_part_list" : ["Parts_Beam"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_support","ROTATION_support"],
        "computing_model_part_name" : "computing_domain",
        "rotation_dofs"                      : true,
        "linear_solver_settings"       : {
            "solver_type"         : "Super_LU"
        },
        "model_import_settings"        : {
            "input_type"     : "mdpa",
            "input_filename" : "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/Beam_structure"
        },
        "material_import_settings" :{
            "materials_filename": "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/materials_beam.json"
        }
    }
```

and the ```list_other_processes``` like

```python
 "loads_process_list"       : [],
    "list_other_processes" :[{
        "kratos_module" : "KratosMultiphysics.HDF5Application",
        "python_module" : "single_mesh_primal_input_process",
        "help"          : "",
        "process_name"  : "",
        "Parameters" : {
	        "model_part_name" : "Structure",
            "file_settings" : {
                "file_access_mode" : "read_only"
            },
            "nodal_results_settings" : {
                "list_of_variables": ["DISPLACEMENT", "ROTATION"]
            }
        }
     }
```     



