 
## Adjoint Sensitivity Analysis

### Theory:


  
  
### Features:  

- Availible response utilities (response functions):
     * Strain energy
     * Displacement or rotation of a node 
     * Stress resultant of a single element
  
- A set of adjoint *Neumann* conditions:
     * Point loads (loads applied directly on the nodes)
     * Surface load (a distributed load applied over a face)
   
- Structural adjoint elements:
    * Uni-dimensional elements :
       	* Linear beam element (3D)
    * Two-dimensional elements :
        * Thin shell (triangular)
       		
- Schemes:
	* Scheme to solve the adjoint problem
	* Eigen solver scheme

- Processes:
    * replacement process (replaces all elements and conditions of a model with its adjoint equivalents and and vice versa)
