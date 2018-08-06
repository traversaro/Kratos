// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_total_lagrangian.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceTotalLagrangian::AdjointFiniteDifferenceTotalLagrangian(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement)
{
}

AdjointFiniteDifferenceTotalLagrangian::~AdjointFiniteDifferenceTotalLagrangian()
{
}

void AdjointFiniteDifferenceTotalLagrangian::EquationIdVector(EquationIdVectorType& rResult,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTotalLagrangian::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTotalLagrangian::GetValuesVector(Vector& rValues, int Step)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
    }
}

void AdjointFiniteDifferenceTotalLagrangian::Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "Adjoint total langrangian is currently not aible to serve with local stresses!" << std::endl;

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceTotalLagrangian::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    GeometryType& r_geom = GetGeometry();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < r_geom.size(); ++i)
    {
        auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node)


        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    }

    // Check constitutive law
    //if ( mConstitutiveLawVector.size() > 0 ) {
    //    return mConstitutiveLawVector[0]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );

    return 0;

    KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTotalLagrangian::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        double L = 0.0;

        const SizeType number_of_nodes = GetGeometry().size();
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            for (IndexType j = i+1; j < number_of_nodes; ++j)
            {
                dx = this->GetGeometry()[i].X0() - this->GetGeometry()[j].X0();
                dy = this->GetGeometry()[i].Y0() - this->GetGeometry()[j].Y0();
                dz = this->GetGeometry()[i].Z0() - this->GetGeometry()[j].Z0();
                L += std::sqrt(dx*dx + dy*dy + dz*dz);
            }
        }
        L /= number_of_nodes*(number_of_nodes-1)/2.0; 
        std::cout << L << std::endl;
        
        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

/*void AdjointFiniteDifferenceTotalLagrangian::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    mpPrimalElement->CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rOutput, rCurrentProcessInfo);
}*/


void AdjointFiniteDifferenceTotalLagrangian::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceTotalLagrangian::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


