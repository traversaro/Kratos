 
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

#### Primal Problem
The primal problem can be defind by the regular input files which are needed for an usual linear static analysis. As only difference the output process of the HDF5Application has to be added in the ```list_other_processes``` in the project parameters:

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

#### Adjoint Problem
In order to define the adjoint problem an additional *.json-file for the adjoint project parameters is necessary. This input file is in principle very similar to the respective file of the primal analysis. In comparsion to a regular file for a linear static analysis three points have to be modified:
- ```solver_settings``` by using the ```adjoint_structural_solver``` as ```solver_type``` and by the definion of the ```response_function_settings```
- 

