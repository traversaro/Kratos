//
//   Project Name:
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $

// External includes

// Project includes

#include "includes/model_part.h"
#include "fem_to_dem_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/stress_to_nodes_process.hpp"
#include "custom_processes/damage_to_nodes_process.hpp"
#include "custom_processes/dem_after_remesh_identificator_process.hpp"
#include "custom_processes/initial_dem_skin_process.hpp"
//#include "custom_processes/skin_nodes_detection_process.hpp"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module &m)
{
	using namespace pybind11;

	typedef Process ProcessBaseType;

	// Stress extrapolation to Nodes
	class_<StressToNodesProcess, StressToNodesProcess::Pointer, Process>(m, "StressToNodesProcess")
		.def(init<ModelPart &, unsigned int>())
		.def("Execute", &StressToNodesProcess::Execute);

	// Damage extrapolation to Nodes
	class_<DamageToNodesProcess, DamageToNodesProcess::Pointer, Process>(m, "DamageToNodesProcess")
		.def(init<ModelPart &, unsigned int>())
		.def("Execute", &DamageToNodesProcess::Execute);

	class_<DemAfterRemeshIdentificatorProcess, DemAfterRemeshIdentificatorProcess::Pointer, Process>(m, "DemAfterRemeshIdentificatorProcess")
		.def(init<ModelPart &>())
		.def("Execute", &DemAfterRemeshIdentificatorProcess::Execute);

	class_<InitialDemSkinProcess, InitialDemSkinProcess::Pointer, Process>(m, "InitialDemSkinProcess")
		.def(init<ModelPart &>())
		.def("Execute", &InitialDemSkinProcess::Execute);

	//class_<SkinNodesDetectionProcess2D, SkinNodesDetectionProcess2D::Pointer, Process>(m, "SkinNodesDetectionProcess2D")
	//	.def(init<ModelPart &>())
	//	.def("Execute", &SkinNodesDetectionProcess2D::Execute);
}
} // namespace Python.
} // Namespace Kratos