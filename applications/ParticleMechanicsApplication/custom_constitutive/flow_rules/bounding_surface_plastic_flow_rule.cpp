//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes
#include <iostream>
#include <cmath>

// External includes
#include "includes/ublas_interface.h"
#include "includes/mat_variables.h"

// Project includes
#include "custom_constitutive/flow_rules/bounding_surface_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "particle_mechanics_application.h"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
BoundingSurfacePlasticFlowRule::BoundingSurfacePlasticFlowRule()
    :MPMFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BoundingSurfacePlasticFlowRule::BoundingSurfacePlasticFlowRule(YieldCriterionPointer pYieldCriterion)
    :MPMFlowRule(pYieldCriterion)
{

}

//********* ASSIGMENT OPERATOR
BoundingSurfacePlasticFlowRule& BoundingSurfacePlasticFlowRule::operator=(BoundingSurfacePlasticFlowRule const& rOther)
{
    MPMFlowRule::operator=(rOther);
    return *this;

}



//********** COPY CONSTRUCTOR *********
BoundingSurfacePlasticFlowRule::BoundingSurfacePlasticFlowRule(BoundingSurfacePlasticFlowRule const& rOther)
    :MPMFlowRule(rOther)
{
}

//*******   CLONE ********
MPMFlowRule::Pointer BoundingSurfacePlasticFlowRule::Clone() const
{
    MPMFlowRule::Pointer p_clone(new BoundingSurfacePlasticFlowRule(*this));
    return p_clone;
}



// ********** DESTRUCTOR **************
BoundingSurfacePlasticFlowRule::~BoundingSurfacePlasticFlowRule()
{
}

void BoundingSurfacePlasticFlowRule::InitializeMaterial(YieldCriterionPointer& pYieldCriterionPointer, HardeningLawPointer& pHardeningPointer, const Properties& rProp)
{
    MPMFlowRule::InitializeMaterial(pYieldCriterionPointer, pHardeningPointer, rProp);

    mElasticPrincipalStrain = ZeroVector(3);
    mPlasticPrincipalStrain = ZeroVector(3);

    mPrincipalStressUpdated = ZeroVector(3);
    mLargeStrainBool = true;
    mRegion = 0;

    // Used to calculate Omega
    mInitialVolumetricStrain = 0.0;

    // Saved state function and its derivative
    mStateFunction = 0.0;
    mStateFunctionFirstDerivative  = ZeroVector(3);
    mStateFunctionSecondDerivative = ZeroVector(6);

    this->InitializeMaterialParameters();
}

// Initiate Material Parameters which are allowed to change
void BoundingSurfacePlasticFlowRule::InitializeMaterialParameters(){

}


bool BoundingSurfacePlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, 
    Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool plasticity_active = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
    
    Vector principal_stress = ZeroVector(3);
    Vector main_strain      = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);

    for(unsigned int i=0; i<3; i++)
    {
        // the rStressMatrix is the precomputed principal stress or trial principal stress
        principal_stress[i] = rStressMatrix(i,i);
    }

    // Sorting Principal Stress and Strain - "0" is the largest one and "2" is the lowest one
    // MPMStressPrincipalInvariantsUtility::SortPrincipalStress(principal_stress, main_strain, rReturnMappingVariables.MainDirections);

    // Assigning to local variables
    mElasticPrincipalStrain = main_strain;

    // Check for the yield Condition -- calling the yield criterion
    // rReturnMappingVariables.TrialStateFunction = 0.0;
    // rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, principal_stress, 0.0, mMaterialParameters.PreconsolidationPressure);
    
    // If yield is reached, do return mapping
    // if (rReturnMappingVariables.TrialStateFunction <= 0.0)
    // {
    //     mRegion = 0;
    //     mPrincipalStressUpdated = principal_stress;
    //     plasticity_active = false;
    //     rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

    //     this->UpdateStateVariables(mPrincipalStressUpdated);
        
    // }
    // else
    // {
    //     unsigned int region = 0;
    //     Vector principal_stress_updated = ZeroVector(3);

    //     // Perform return mapping to the yield surface: Will update mElasticPrincipalStrain, Region, and principal_stress_updated
    //     bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, principal_stress, mElasticPrincipalStrain, region, principal_stress_updated);
    //     KRATOS_ERROR_IF(!converged) << "Warning:: Constitutive Law does not converge! "<<std::endl;

    //     mRegion = region;
    //     mPrincipalStressUpdated = principal_stress_updated;

    //     plasticity_active = true;
    //     rReturnMappingVariables.Options.Set(PLASTIC_REGION,true);
    // }

    // rStressMatrix is the matrix of the updated stress in cartesian configuration -- this function perform back transformation
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, mPrincipalStressUpdated, rStressMatrix);
 
    // Delta plastic strain
    mPlasticPrincipalStrain = main_strain - mElasticPrincipalStrain ;

    // We're saving the updated info in terms of principal strain and stress in these matrix
    // these information will be used for the evaluation of the second contribution in the
    // consistent tangent matrix
    for (unsigned int i=0; i<3; i++)
    {
        rReturnMappingVariables.StrainMatrix(i,i) = mElasticPrincipalStrain[i];
        rReturnMappingVariables.TrialIsoStressMatrix(i,i) = mPrincipalStressUpdated[i];
    }

    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

    return plasticity_active;

}

bool BoundingSurfacePlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, 
    Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated)
{

    // Calculate stress return in principal stress space
    // The flow rule is written for non-associated plasticity and explicit assumption using image point
    // Refer to paper by (Russel&Khalili, 2003) for the theoretical description
    bool converged = false;

    return converged;
}

void BoundingSurfacePlasticFlowRule::CalculatePlasticPotentialDerivatives(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector, Vector& rFirstDerivative)
{
    double mean_stress_p, deviatoric_q, lode_angle;
    MPMStressPrincipalInvariantsUtility::CalculateStressInvariants(rPrincipalStressVector, mean_stress_p, deviatoric_q, lode_angle);
    mean_stress_p *= -1.0;  // p is defined negative

    // Get material parameters
    const double parameter_A = GetProperties()[MODEL_PARAMETER_A];
    const double shear_M     = this->CalculateCriticalStateLineSlope(lode_angle);
    const double direction_T = this->GetDirectionParameter(rPrincipalStressVector, rImagePointPrincipalStressVector); 
    const double alpha       = this->GetAlphaParameter();

    rFirstDerivative = ZeroVector(3);
    rFirstDerivative[0]  = parameter_A * (shear_M - direction_T * (deviatoric_q/mean_stress_p));
    rFirstDerivative[1]  = direction_T;
    rFirstDerivative[2]  = - direction_T * 3.0/4.0 * deviatoric_q;
    rFirstDerivative[2] *= (1.0 - std::pow(alpha, 4)) * std::cos(3.0 * lode_angle) / (1.0 + std::pow(alpha, 4) - (1 - std::pow(alpha, 4)) * std::sin(3.0 * lode_angle) );

}

void BoundingSurfacePlasticFlowRule::ComputeElasticMatrix_3X3(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, Matrix& rElasticMatrix)
{
    // Assemble rElasticMatrix matrix
    rElasticMatrix = ZeroMatrix(3);
}

void BoundingSurfacePlasticFlowRule::ComputePlasticMatrix_3X3(const Vector& rPrincipalStressVector, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Matrix& rElasticMatrix, Matrix& rPlasticMatrix)
{
    // Arrange rPlasticMatrix matrix
    rPlasticMatrix = ZeroMatrix(3);
}

// Compute Trial elastic principal stress matrix from Trial elastic principal strain matrix
void BoundingSurfacePlasticFlowRule::CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{
    Vector main_strain      = ZeroVector(3);
    Vector principal_stress = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
    {
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    this->CalculatePrincipalStressVector(main_strain, principal_stress);

    // Evalute the Kirchhoff principal stress
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = principal_stress[i];
    }

}

// Function to compute Principal Stress Vector from Principal Strain Vector
void BoundingSurfacePlasticFlowRule::CalculatePrincipalStressVector(Vector& rPrincipalStrain, Vector& rPrincipalStress)
{

}

// Function which returns principal strains from volumetric and deviatoric strain components
void BoundingSurfacePlasticFlowRule::CalculatePrincipalStrainFromStrainInvariants(Vector& rPrincipalStrain, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Vector& rDirectionVector)
{
    rPrincipalStrain = ZeroVector(3);
    
    for (unsigned int i = 0; i<3; ++i)
    {
        rPrincipalStrain[i] += 1.0/3.0 * rVolumetricStrain;
    }
    rPrincipalStrain += std::sqrt(3.0/2.0) * rDeviatoricStrain * rDirectionVector;
}


// Function which returns volumetric and deviatoric strain components from principal strain
void BoundingSurfacePlasticFlowRule::CalculateStrainInvariantsFromPrincipalStrain(const Vector& rPrincipalStrain, double& rVolumetricStrain, double& rDeviatoricStrain, Vector& rDeviatoricStrainVector)
{
    rDeviatoricStrainVector = rPrincipalStrain;
    rVolumetricStrain = sum(rPrincipalStrain);
    for (unsigned int i = 0; i<3; ++i)
    {
        rDeviatoricStrainVector[i] -= 1.0/3.0 * rVolumetricStrain;
    }
    rDeviatoricStrain = std::sqrt(2.0/3.0) * norm_2(rDeviatoricStrainVector);
}

// Function to return matrix from principal space to normal space
void BoundingSurfacePlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
{
    rStressMatrix = ZeroMatrix(3,3); 
    Vector aux_N  = ZeroVector(3);
    Matrix aux_M  = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
            aux_N[j] = rEigenVectors(j,i);
        aux_M = MathUtils<double>::TensorProduct3(aux_N, aux_N);
        rStressMatrix += rPrincipalStress[i]*aux_M;
    }
}

// Function that compute the consistent tangent stiffness matrix (in normal space) considering both elastic and elasto-plastic case
void BoundingSurfacePlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix)
{
    
    // Compute Consistent Tangent Stiffness matrix in principal space
    Matrix D_elasto_plastic = ZeroMatrix(6,6);

    // Return constitutive matrix from principal space to normal space
    Matrix A = ZeroMatrix(6,6);
    Matrix A_trans = ZeroMatrix(6,6); 
    this->CalculateTransformationMatrix(rReturnMappingVariables.MainDirections, A);
    A_trans = trans(A);

    Matrix aux_mat = ZeroMatrix(6,6);
    aux_mat = prod(A_trans, D_elasto_plastic);
    rConsistMatrix = prod(aux_mat, A);

}

void BoundingSurfacePlasticFlowRule::CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA)
{
    Matrix A1 = ZeroMatrix(3,3);
    Matrix A2 = ZeroMatrix(3,3);
    Matrix A3 = ZeroMatrix(3,3);
    Matrix A4 = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3 ; i++)
    {
        for(unsigned int j = 0; j<3 ; j++)
        {
            A1(i,j) = rMainDirection(i,j) * rMainDirection(i,j);
            rA(i,j) = A1(i,j);
        }
    }
    Vector Hj1 = ZeroVector(3);
    Hj1[0] = 0;
    Hj1[1] = 2;
    Hj1[2] = 1;

    Vector Hj2 = ZeroVector(3);
    Hj2[0] = 1;
    Hj2[1] = 0;
    Hj2[2] = 2;

    for(unsigned int k = 0; k<3; k++)
    {
        for(unsigned int l = 0; l<3; l++)
        {
            A2(k,l) = rMainDirection(k,Hj1[l]) * rMainDirection(k,Hj2[l]);
            A3(k,l) = rMainDirection(Hj1[k],l) * rMainDirection(Hj2[k],l);
            A4(k,l) = rMainDirection(Hj1[k],Hj1[l]) * rMainDirection(Hj2[k],Hj2[l]) + rMainDirection(Hj2[k],Hj1[l]) * rMainDirection(Hj1[k],Hj2[l]);
        }
    }

    for(unsigned int i = 0 ; i<3 ; i++)
    {
        for(unsigned int j = 3; j<6 ; j++)
        {
            unsigned int index_j = j - 3;
            unsigned int index_i = i + 3;

            rA(i,j) = A2(i, index_j);
            rA(index_i, index_j) = A3(i, index_j);
            rA(index_i, j) = A4(i, index_j);
        }
    }

    rA = trans(rA);
}


// Function that gives the elastic left cauchy green tensor B
Matrix BoundingSurfacePlasticFlowRule::GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables)
{
    Vector landa_2 = ZeroVector(3);

    for (unsigned int i = 0; i<3; ++i)
        landa_2[i] = std::exp(2.0*mElasticPrincipalStrain[i]);

    Matrix output = ZeroMatrix(3,3);
    this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, landa_2, output);

    return output;
}

// Function that updates internal variable at every time step once the nonlinear iteration converges
bool BoundingSurfacePlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    // Compute Delta Plastic Strain
    double norm_plastic_principal_strain = norm_2(mPlasticPrincipalStrain);

    // Compute Strain Components and its invariants
    double volumetric_strain, deviatoric_strain;
    Vector deviatoric_strain_vector;
    this->CalculateStrainInvariantsFromPrincipalStrain(mPlasticPrincipalStrain, volumetric_strain, deviatoric_strain, deviatoric_strain_vector);

    // Update Equivalent Plastic Strain
    mInternalVariables.DeltaPlasticStrain = norm_plastic_principal_strain;
    mInternalVariables.EquivalentPlasticStrain += norm_plastic_principal_strain;

    // Update Accumulated Plastic Volumetric Strain
    mInternalVariables.DeltaPlasticVolumetricStrain = volumetric_strain;
    mInternalVariables.AccumulatedPlasticVolumetricStrain += volumetric_strain;

    // Update Accumulated Plastic Deviatoric Strain
    mInternalVariables.DeltaPlasticDeviatoricStrain = deviatoric_strain;
    mInternalVariables.AccumulatedPlasticDeviatoricStrain += deviatoric_strain;

    return true;
}

double BoundingSurfacePlasticFlowRule::CalculateCriticalStateLineSlope(const double& rLodeAngle)
{
    double shear_M = GetProperties()[CRITICAL_STATE_LINE];
    const double alpha = this->GetAlphaParameter();

    const double aux_multiplier = 2.0 * std::pow(alpha, 4) / (1.0 + std::pow(alpha, 4) - (1 - std::pow(alpha, 4)) * std::sin(3.0 * rLodeAngle) );
    shear_M *= std::pow(aux_multiplier, 0.25);

    KRATOS_ERROR_IF(shear_M < 0.0) << "The slope of critical state line is negative! M_cs = " << shear_M << std::endl;

    return shear_M;
}

double BoundingSurfacePlasticFlowRule::GetDirectionParameter(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector)
{
    double direction_T;

    // Compute angle of stress reference axis and the stress point
    const double gamma_angle    = this->GetDirectionAngle(rPrincipalStressVector);
    const double gamma_angle_IP = this->GetDirectionAngle(rImagePointPrincipalStressVector);

    if (std::abs(gamma_angle - gamma_angle_IP) < 0.5 * this->GetPI())
        direction_T = 1.0;
    else if (std::abs(gamma_angle - gamma_angle_IP) >= 0.5 * this->GetPI())
        direction_T = -1.0;
    else
        KRATOS_ERROR << "There is a problem with angles gamma_angle = " << gamma_angle <<  ", gamma_angle_IP = " << gamma_angle_IP << std::endl;

    return direction_T;
}

double BoundingSurfacePlasticFlowRule::GetDirectionAngle(const Vector& rPrincipalStressVector)
{
    double angle = std::acos(1/std::sqrt(3.0));
    const double aux   = 2.0 * rPrincipalStressVector[2] - rPrincipalStressVector[1] - rPrincipalStressVector[0];
    const double aux_2 = rPrincipalStressVector[1] - rPrincipalStressVector[0];
    if(std::abs(aux) > 0.0)
    {
        if(rPrincipalStressVector[2] >= 0.0)
            angle = std::atan(std::sqrt(3.0)*aux_2/aux) + 2.0 * GetPI();
        else
            angle = std::atan(std::sqrt(3.0)*aux_2/aux) + GetPI();
    }

    return angle;
}

double BoundingSurfacePlasticFlowRule::GetAlphaParameter()
{
    const double shear_M = GetProperties()[CRITICAL_STATE_LINE];
    const double phi_csl = (3.0 * shear_M) / (6.0 + shear_M);
    const double alpha   = (3.0 - std::sin(phi_csl)) / (3.0 + std::sin(phi_csl));
    
    return alpha;
}

double BoundingSurfacePlasticFlowRule::GetPI()
{
    return std::atan(1.0)*4.0;
}

unsigned int BoundingSurfacePlasticFlowRule::GetPlasticRegion()
{
    return mRegion;
}

void BoundingSurfacePlasticFlowRule::save( Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMFlowRule )
}

void BoundingSurfacePlasticFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMFlowRule )

}

} //end namespace kratos
