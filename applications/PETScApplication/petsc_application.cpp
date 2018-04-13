//
//   Project Name:        KratosPETScApplication      $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:           April 2018 $
//   Revision:            $Revision:              0.0 $
//
//

// System includes

// External includes

// Project includes
#include "petsc_application.h"
#include "petsc_application_variables.h"


namespace Kratos {

KratosPETScApplication::KratosPETScApplication():
    KratosApplication("PETScApplication")
    {}

void KratosPETScApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosPETScApplication... " << std::endl;


}
}  // namespace Kratos.
