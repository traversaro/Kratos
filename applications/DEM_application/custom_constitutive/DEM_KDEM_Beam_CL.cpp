// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Beam_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Beam::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Beam(*this));
        return p_clone;
    }

    void DEM_KDEM_Beam::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Beam to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM_Beam::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                         SphericContinuumParticle* neighbor,
                                                         double equiv_young,
                                                         double distance,
                                                         double calculation_area,
                                                         double LocalCoordSystem[3][3],
                                                         double ElasticLocalRotationalMoment[3],
                                                         double ViscoLocalRotationalMoment[3],
                                                         double equiv_poisson,
                                                         double indentation) {

        KRATOS_TRY
        double LocalDeltaRotatedAngle[3]    = {0.0};
        double LocalDeltaAngularVelocity[3] = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3> GlobalDeltaAngularVelocity;
        noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);

        const double element_mass  = element->GetMass();
        const double neighbor_mass = neighbor->GetMass();
        const double equiv_mass    = element_mass * neighbor_mass / (element_mass + neighbor_mass);
        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double Inertia_Ix = std::max(element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X], neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X]);
        const double Inertia_Iy = std::max(element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y], neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y]);
        const double Inertia_J  = Inertia_Ix + Inertia_Iy;

        const double my_gamma    = element->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = neighbor->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        double aux = (element->GetRadius() + neighbor->GetRadius()) / distance; // This is necessary because if spheres are not tangent the DeltaAngularVelocity has to be interpolated

        //Viscous parameter taken from Olmedo et al., 'Discrete element model of the dynamic response of fresh wood stems to impact'
        array_1d<double, 3> visc_param;
        visc_param[0] = 2.0 * equiv_gamma * std::sqrt(equiv_mass * equiv_young * Inertia_Ix / aux); // OLMEDO
        visc_param[1] = 2.0 * equiv_gamma * std::sqrt(equiv_mass * equiv_young * Inertia_Iy / aux); // OLMEDO
        visc_param[2] = 2.0 * equiv_gamma * std::sqrt(equiv_mass * equiv_shear * Inertia_J  / aux); // OLMEDO

        array_1d<double, 3> LocalEffDeltaRotatedAngle;
        LocalEffDeltaRotatedAngle[0] = LocalDeltaRotatedAngle[0] * aux;
        LocalEffDeltaRotatedAngle[1] = LocalDeltaRotatedAngle[1] * aux;
        LocalEffDeltaRotatedAngle[2] = LocalDeltaRotatedAngle[2] * aux;

        array_1d<double, 3> LocalEffDeltaAngularVelocity;
        LocalEffDeltaAngularVelocity[0] = LocalDeltaAngularVelocity[0] * aux;
        LocalEffDeltaAngularVelocity[1] = LocalDeltaAngularVelocity[1] * aux;
        LocalEffDeltaAngularVelocity[2] = LocalDeltaAngularVelocity[2] * aux;

        ElasticLocalRotationalMoment[0] = -equiv_young * Inertia_Ix * LocalEffDeltaRotatedAngle[0] / distance;
        ElasticLocalRotationalMoment[1] = -equiv_young * Inertia_Iy * LocalEffDeltaRotatedAngle[1] / distance;
        ElasticLocalRotationalMoment[2] = -equiv_shear * Inertia_J  * LocalEffDeltaRotatedAngle[2] / distance;

        ViscoLocalRotationalMoment[0] = -visc_param[0] * LocalEffDeltaAngularVelocity[0];
        ViscoLocalRotationalMoment[1] = -visc_param[1] * LocalEffDeltaAngularVelocity[1];
        ViscoLocalRotationalMoment[2] = -visc_param[2] * LocalEffDeltaAngularVelocity[2];

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
