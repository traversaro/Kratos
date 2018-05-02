 
## Adjoint Sensitivity Analysis

### General remarks:

This feature of the Structural Mechanics Applicattion provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy, stresses) with respect to different kinds of design variables (e.g. nodal coordinates, material or cross-sectional properties) with the adjoint approach. Therefore for each response function an adjoint problem has to be solved, which is the solution of a simple linear static problem. The sensitivies are than computed in a post-processing step. The implemented sensitivity analysis uses a so called semi-analytic approach which means that the derivatives at element level are computed by finite differences.
  
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
    * replacement process (replaces all elements and conditions of a model with its adjoint equivalents and and vice versa)

- A set of adjoint *Neumann* conditions:
     * Point loads (derived from PointLoadCondition)
     * Surface load (derived from SurfaceLoadCondition3D)
   
- Structural adjoint elements:
    * Uni-dimensional elements :
       	* Linear 3D beam element (derived from CrBeamElementLinear3D2N)
    * Two-dimensional elements :
        * Thin triangular shell (derived from ShellThinElement3D3N)

*Please note:* The adjoint elements and conditions are derived from elements/conditions of the Structural Mechanics Application and can not be seen independently from them. Rather they have to be traced as additions of their parents in order to use the parents in the context of adjoint sensitivity analysis. So basic tasks like the computation of the stiffness matrix are not overwritten. The main task of the adjoint elements/conditions is the derive different quantities with respect to the design varibale or state (e.g. the right hand side or post-processing results like stresses).

### Usage:  	    

