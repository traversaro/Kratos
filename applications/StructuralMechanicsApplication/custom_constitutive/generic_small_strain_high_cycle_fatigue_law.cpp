// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
//  Collaborator:    
//

// System includes

// External includes

// Project includes
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/generic_small_strain_high_cycle_fatigue_law.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/high_cycle_fatigue_law_integrator.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    Vector& integrated_stress_vector = rValues.GetStressVector();
    array_1d<double, VoigtSize> auxiliar_integrated_stress_vector = integrated_stress_vector;
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    // We compute the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        // Converged values
        double& threshold = this->GetThreshold();
        double& damage = this->GetDamage();

        // S0 = C:E
        array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        double min_stress = 0.0, max_stress = 0.0, sign_factor;
        HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactor(predictive_stress_vector, sign_factor);
        uniaxial_stress *= sign_factor;
        unsigned int number_of_cycles = this->GetNumberOfCycles();
        HighCycleFatigueLawIntegrator<6>::CalculateMaximumAndMinimumStresses(uniaxial_stress, max_stress, min_stress, 
                                                                            this->GetPreviousStresses(), number_of_cycles);
        //this->SetNumberOfCycles(number_of_cycles);

        unsigned int aux = 0;
        if (std::abs(max_stress) > 0.0 || max_stress != this->GetMaxStress()) {
            this->SetMaxStress(max_stress);
        }
        if (std::abs(min_stress) > 0.0 || min_stress != this->GetMinStress()) {
            this->SetMinStress(min_stress);
        }

        double fatigue_reduction_factor = this->GetFatigueReductionFactor();
        if (number_of_cycles > this->GetNumberOfCycles() && max_stress > 0.0 && std::abs(min_stress) > 0.0) {
            double reversion_factor = this->GetReversionFactor();
            double B0 = this->GetFatigueReductionParameter();
            HighCycleFatigueLawIntegrator<6>::CalculateFatigueReductionFactor(this->GetMaxStress(),
                                                                              this->GetMinStress(),
                                                                              reversion_factor,
                                                                              rValues.GetMaterialProperties(),
                                                                              this->GetNumberOfCycles(),
                                                                              fatigue_reduction_factor,
                                                                              B0);
            this->SetReversionFactor(reversion_factor);
            this->SetFatigueReductionParameter(B0);
            this->SetFatigueReductionFactor(fatigue_reduction_factor);
        }
        
        uniaxial_stress /= fatigue_reduction_factor;  // Fatigue contribution
        const double F = uniaxial_stress - threshold;

        if (F <= 0.0) { // Elastic case
            this->SetNonConvDamage(damage);
            this->SetNonConvThreshold(threshold);
            noalias(auxiliar_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;
			noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;

                TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(auxiliar_integrated_stress_vector, r_strain_vector, uniaxial_stress, rValues);

                this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
            }
        } else { // Damage case
            const double characteristic_length = rValues.GetElementGeometry().Length();
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, damage, threshold, rValues, characteristic_length);

            // Updated Values
            noalias(auxiliar_integrated_stress_vector) = predictive_stress_vector;
            this->SetNonConvDamage(damage);
            this->SetNonConvThreshold(uniaxial_stress);

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
                noalias(r_tangent_tensor) = rValues.GetConstitutiveMatrix();

                TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(auxiliar_integrated_stress_vector, r_strain_vector, uniaxial_stress, rValues);

                this->SetValue(UNIAXIAL_STRESS, uniaxial_stress, rValues.GetProcessInfo());
            }
            noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->SetDamage(this->GetNonConvDamage());
    this->SetThreshold(this->GetNonConvThreshold());

    Vector previous_stresses = ZeroVector(2);
    previous_stresses[1] = this->GetValue(UNIAXIAL_STRESS,previous_stresses[1]);
    previous_stresses[0] = this->GetPreviousStresses()[1];
    this->SetPreviousStresses(previous_stresses);
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<TrescaPlasticPotential<6>>>>;

} // namespace Kratos