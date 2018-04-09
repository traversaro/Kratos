// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#include "custom_frictional_laws/base_frictional_law.h"

namespace Kratos
{
BaseFrictionalLaw::BaseFrictionalLaw()
{

}

/***********************************************************************************/
/***********************************************************************************/

BaseFrictionalLaw::Pointer BaseFrictionalLaw::Clone() const
{
    KRATOS_ERROR << "Called the virtual function for Clone" << std::endl;;
}

/***********************************************************************************/
/***********************************************************************************/

double BaseFrictionalLaw::GetFrictionalCoefficient(NodeType& ThisNode)
{
    return 0.0;
}

} // namespace Kratos

