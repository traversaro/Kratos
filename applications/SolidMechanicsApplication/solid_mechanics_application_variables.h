//------------------------------------------------------------------
//           ___      _ _    _                                     .
//   KRATOS / __| ___| (_)__| |                                    .
//          \__ \/ _ \ | / _` |                                    .
//          |___/\___/_|_\__,_| MECHANICS                          .
//			                                           .
//   License:(BSD)	  SolidMechanicsApplication/license.txt    .
//   Main authors:        Josep Maria Carbonell                    .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED)
#define  KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/mat_variables.h"
#include "custom_solvers/time_integration_methods/time_integration_methods_container.hpp"
#include "custom_utilities/shell_cross_section.hpp"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;
  typedef TimeIntegrationMethodsContainer                                TimeIntegrationContainerType;      
  typedef TimeIntegrationContainerType::Pointer                   TimeIntegrationContainerPointerType;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Define Variables

  // Generalized eigenvalue problem
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int, BUILD_LEVEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, EIGENVALUE_VECTOR)
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix , EIGENVECTOR_MATRIX )

  //for integration methods
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, TimeIntegrationContainerPointerType, TIME_INTEGRATION_METHODS )  
  
  //for explicit schemes
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, MIDDLE_VELOCITY )

  //solution
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int, WRITE_ID )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int, TIME_INTEGRATION_ORDER )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, RAYLEIGH_ALPHA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, RAYLEIGH_BETA )

  //geometrical
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, GEOMETRIC_STIFFNESS )

  //beam cross section    
  //KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, BeamCrossSection::Pointer, SBEAM_CROSS_SECTION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, CROSS_SECTION_AREA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, CROSS_SECTION_RADIUS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int,    CROSS_SECTION_SIDES )
    
  //shell cross section
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int, SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )
      
  //shell generalized variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_STRAIN )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_STRAIN_GLOBAL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_CURVATURE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_CURVATURE_GLOBAL )      
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_FORCE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_FORCE_GLOBAL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_MOMENT )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, SHELL_MOMENT_GLOBAL )
    
  //nodal load variables (legacy)
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, POINT_LOAD )

  //force loads
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FORCE_LOAD )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, FORCE_LOAD_VECTOR )

  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FOLLOWER_FORCE_LOAD )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, FOLLOWER_FORCE_LOAD_VECTOR )

  //moment loads
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, MOMENT_LOAD )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, MOMENT_LOAD_VECTOR )

  //elastic loads
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, ELASTIC_LOAD )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, ELASTIC_LOAD_VECTOR )

  //force pressure
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, POSITIVE_FACE_PRESSURE_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, NEGATIVE_FACE_PRESSURE_VECTOR )

  //moment pressures
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, PLANE_MOMENT_LOAD )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, PLANE_MOMENT_LOAD_VECTOR )
  
  //elastic pressures
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, BALLAST_COEFFICIENT )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, BALLAST_COEFFICIENT_VECTOR )
    
  //element  
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, VON_MISES_STRESS )

  //nodal dofs
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, PRESSURE_REACTION )    

  //explicit beam
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, EXTERNAL_MOMENT )
    
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, POSITION_MOMENTUM )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, ROTATION_MOMENTUM )  

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, INERTIA_DYADIC ) 

  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, RESIDUAL_LYAPUNOV )  
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, TANGENT_MATRIX )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, TANGENT_LYAPUNOV )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, ALPHA_TRAPEZOIDAL_RULE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, bool, POSITION_UPDATE_LABEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, bool, ROTATION_UPDATE_LABEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, bool, MOMENTUM_UPDATE_LABEL )

  //reading beam section properties
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SECTION_HEIGHT )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SECTION_WIDTH  )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, INERTIA_X )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, INERTIA_Y )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SECTION_SIZE )
    
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, YOUNGxAREA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, YOUNGxINERTIA_X )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, YOUNGxINERTIA_Y )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SHEARxREDUCED_AREA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, SHEARxPOLAR_INERTIA )
    
  ///@}

}

#endif	/* KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
