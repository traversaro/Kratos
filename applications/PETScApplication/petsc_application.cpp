//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
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
