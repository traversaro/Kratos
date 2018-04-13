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
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_strategies_to_python.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace pybind11;

void AddCustomStrategiesToPython(pybind11::module& m)
{
  typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

  typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
  typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
  typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

  
}

}  // namespace Python.

} // Namespace Kratos
