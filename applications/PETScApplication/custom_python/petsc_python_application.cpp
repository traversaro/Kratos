//
//   Project Name:        KratosPETScApplication      $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:           April 2018 $
//   Revision:            $Revision:              0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "petsc_application.h"

namespace Kratos
{

    namespace Python
    {

	using namespace pybind11;


	PYBIND11_MODULE(KratosPETScApplication,m)
	{

	    class_<KratosPETScApplication,
		   KratosPETScApplication::Pointer,
		   KratosApplication>(m,"KratosPETScApplication")
		.def(init<>())
		;

	    AddCustomStrategiesToPython(m);
	    AddCustomUtilitiesToPython(m);

	    //registering variables in pythonn ( if must to be seen from python )

	}


    }  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
