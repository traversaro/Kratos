//
// Authors: Joaquín Irazábal jirazabal@cimne.upc.edu
//

#include "beam_particle.h"
#include <cmath>

namespace Kratos {

    void BeamParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbors of my neighbors.
    {
        double alpha = 1.0;
        //double external_sphere_area = 4 * Globals::Pi * GetRadius()*GetRadius();
        double effectiveVolumeRadius = EffectiveVolumeRadius();  //calculateEffectiveVolumeRadius
        double external_sphere_area = 4 * Globals::Pi * effectiveVolumeRadius * effectiveVolumeRadius;
        double total_equiv_area = 0.0;
        double total_mContIniNeighArea = 0.0;
        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;
        Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);
        bool print_debug_files = false;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();
            double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, cont_ini_neigh_area); //This call fills the vector of areas only if the Constitutive Law wants.
            total_equiv_area += area;

        }
        if (print_debug_files) {
            std::ofstream outputfile("external_sphere_area-total_equiv_area.txt", std::ios_base::out | std::ios_base::app);
            outputfile << external_sphere_area << "  " << total_equiv_area << "\n";
            outputfile.close();
        }
        if (cont_ini_neighbours_size >= 6) {
            if (!IsSkin()) {
                AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha);
                if (print_debug_files) {
                    std::ofstream outputfile("alpha.txt", std::ios_base::out | std::ios_base::app);
                    outputfile << alpha << "\n";
                    outputfile.close();
                }
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                    total_mContIniNeighArea += cont_ini_neigh_area[i];
                } //for every neighbor

                if (print_debug_files) {
                    std::ofstream outputfile2("total_mContIniNeighArea-total_equiv_area.txt", std::ios_base::out | std::ios_base::app);
                    outputfile2 << total_mContIniNeighArea << "  " << total_equiv_area << "\n";
                    outputfile2.close();
                }
            }

            else {//skin sphere
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11.0);

                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                } //for every neighbor
            }
        }//if more than 3 neighbors
    }

    void BeamParticle::ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                 ProcessInfo& r_process_info,
                                                                 array_1d<double, 3>& rElasticForce,
                                                                 array_1d<double, 3>& rContactForce,
                                                                 double& RollingResistance)
    {
        KRATOS_TRY

        const int time_steps = r_process_info[TIME_STEPS];
        const int& search_control = r_process_info[SEARCH_CONTROL];
        DenseVector<int>& search_control_vector = r_process_info[SEARCH_CONTROL_VECTOR];

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        Vector& cont_ini_neigh_area            = this->GetValue(NEIGHBOURS_CONTACT_AREAS);
        int NeighbourSize = mNeighbourElements.size();
        GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            data_buffer.mpOtherParticle = neighbour_iterator;

            unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();

            noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();

            const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

            data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
            double radius_sum = GetRadius() + other_radius;

            double initial_delta = GetInitialDelta(i);

            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - data_buffer.mDistance;
            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el = 0.0;
            double kt_el = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            double LocalCoordSystem[3][3]    = {{0.0}, {0.0}, {0.0}};
            double OldLocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbor properties
            double other_young = data_buffer.mpOtherParticle->GetYoung();
            double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
                mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator);
            }

            EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, vel, delta_displ);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, data_buffer.mDt, ang_vel, neighbour_iterator);
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, OldLocalCoordSystem, LocalCoordSystem, neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double LocalElasticExtraContactForce[3] = {0.0};
            double GlobalElasticContactForce[3] = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};

            FilterNonSignificantDisplacements(DeltDisp, RelVel, indentation);

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);

            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticContactForces[i]);
            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};
            double cohesive_force =  0.0;
            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

            if (i < mContinuumInitialNeighborsSize) {

                mContinuumConstitutiveLawArray[i]->CheckFailure(i, this, neighbour_iterator);

                mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info,
                                                                   OldLocalElasticContactForce,
                                                                   LocalElasticContactForce,
                                                                   LocalElasticExtraContactForce,
                                                                   LocalCoordSystem,
                                                                   LocalDeltDisp,
                                                                   kn_el,
                                                                   kt_el,
                                                                   contact_sigma,
                                                                   contact_tau,
                                                                   failure_criterion_state,
                                                                   equiv_young,
                                                                   equiv_shear,
                                                                   indentation,
                                                                   calculation_area,
                                                                   acumulated_damage,
                                                                   this,
                                                                   neighbour_iterator,
                                                                   i,
                                                                   r_process_info[TIME_STEPS],
                                                                   sliding,
                                                                   search_control,
                                                                   search_control_vector,
                                                                   equiv_visco_damp_coeff_normal,
                                                                   equiv_visco_damp_coeff_tangential,
                                                                   LocalRelVel,
                                                                   ViscoDampingLocalContactForce);

            } else if (indentation > 0.0) {
                const double previous_indentation = indentation + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce,
                        LocalDeltDisp, LocalRelVel, indentation, previous_indentation,
                        ViscoDampingLocalContactForce, cohesive_force, this, data_buffer.mpOtherParticle, sliding, LocalCoordSystem);
            } else { //Not bonded and no idata_buffer.mpOtherParticlendentation
                LocalElasticContactForce[0] = 0.0;      LocalElasticContactForce[1] = 0.0;      LocalElasticContactForce[2] = 0.0;
                ViscoDampingLocalContactForce[0] = 0.0; ViscoDampingLocalContactForce[1] = 0.0; ViscoDampingLocalContactForce[2] = 0.0;
                cohesive_force= 0.0;
            }

            // Transforming to global forces and adding up
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it.
                mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, LocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator, r_process_info, i, indentation);
            }

            array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces);

            AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                                  GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce,ViscoDampingLocalContactForce, 0.0, other_ball_to_ball_forces, rElasticForce, rContactForce, i, r_process_info); //TODO: replace the 0.0 with an actual cohesive force for discontinuum neighbours

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                ComputeMoments(LocalContactForce[2], TotalGlobalElasticContactForce, RollingResistance, LocalCoordSystem[2], data_buffer.mpOtherParticle, indentation);
                if (i < mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {
                    mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, data_buffer.mDistance, calculation_area,
                                                                                        LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment, equiv_poisson, indentation);
                }

                AddUpMomentsAndProject(LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
                double total_local_elastic_contact_force[3] = {0.0};
                total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
                total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
                total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
                CalculateOnContactElements(i, total_local_elastic_contact_force, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
            }

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) /*&& (i < mContinuumInitialNeighborsSize)*/) {
                AddNeighbourContributionToStressTensor(TotalGlobalElasticContactForce, LocalCoordSystem[2], data_buffer.mDistance, radius_sum, this);
            }

            AddContributionToRepresentativeVolume(data_buffer.mDistance, radius_sum, calculation_area);

            ComputeForceWithNeighbourFinalOperations();

            /*if (i < mContinuumInitialNeighborsSize) {
                DEM_COPY_SECOND_TO_FIRST_3(mArrayOfDeltaDisplacements[i], DeltDisp);
            }*/
        } // for each neighbor

        ComputeBrokenBondsRatio();

        KRATOS_CATCH("")
    } //  ComputeBallToBallContactForce

} // namespace Kratos
