//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Josep Maria Carbonell
//


#if !defined(KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED )
#define  KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "custom_solvers/solution_schemes/dynamic_scheme.hpp"
#include "utilities/coordinate_transformation_utilities.h"
//#include "utilities/dof_updater.h"

namespace Kratos {

    /**@name Kratos Globals */
    /*@{ */


    /*@} */
    /**@name Type Definitions */
    /*@{ */

    /*@} */


    /**@name  Enum's */
    /*@{ */


    /*@} */
    /**@name  Functions */
    /*@{ */



    /*@} */
    /**@name Kratos Classes */
    /*@{ */

    /// Bossak time scheme for the incompressible flow problem.
    /** This class provides a second order time scheme of the generalized-alpha Newmark
        family of methods. It also includes code required to implement slip conditions
        on the incompressible flow problem and provides the possibility of using a RANS
        model by passing a turbulence model as an argument to the constructor.
        This time scheme is intended to be used in combination with elements of type
        ASGS2D, ASGS3D, VMS or derived classes.
        To use the slip condition, assign IS_STRUCTURE != 0.0 to the non-historic database
        of the relevant boundary nodes (that is, use SetValue(IS_STRUCTURE,...)). To use
        a wall law in combination with the slip condition, use MonolithicWallCondition to
        mesh the boundary
        @see ASGS2D, ASGS3D, VMS, MonolithicWallConditon
     */
    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class AleSolutionScheme : public DynamicScheme<TSparseSpace, TDenseSpace> {
    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION(AleSolutionScheme);

        typedef DynamicScheme<TSparseSpace, TDenseSpace> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef Element::GeometryType  GeometryType;


        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /** Constructor without a turbulence model
         */
        AleSolutionScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods, Flags& rOptions, unsigned int DomainSize)
            : BaseType(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods, Flags& rOptions),
            mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
        {

        }

        /** Destructor.
         */
        ~AleSolutionScheme() override {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        /**
        * Performing the update of the solution
        * @param rModelPart: The model of the problem to solve
        * @param rDofSet: Set of all primary variables
        * @param rDx: incremental update of primary variables
        */

        void Update(ModelPart& rModelPart,
            DofsArrayType& rDofSet,
            SystemVectorType& rDx) override
        {
            KRATOS_TRY;

            mRotationTool.RotateVelocities(rModelPart);

            this->UpdateDofs(rModelPart,rDofSet,rDx);
            //mpDofUpdater->UpdateDofs(rDofSet,rDx); //TODO: delete

            mRotationTool.RecoverVelocities(rModelPart);

            this->UpdateVariables(rModelPart);

            //TODO: update mesh_acceleration, mesh_velocity, mesh_displacement

            this->MoveMesh(rModelPart);

            KRATOS_CATCH( "" );
        }

        //***************************************************************************

        // void AdditionalUpdateOperations(ModelPart& rModelPart,
        //                                 DofsArrayType& rDofSet,
        //                                 TSystemMatrixType& A,
        //                                 TSystemVectorType& Dv,
        //                                 TSystemVectorType& b)
        // {
        //     KRATOS_TRY

        //     int NumThreads = OpenMPUtils::GetNumThreads();
        //     OpenMPUtils::PartitionVector NodePartition;
        //     OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

        //     //updating time derivatives (nodally for efficiency)
        //     #pragma omp parallel
        //     {
        //         array_1d<double, 3 > DeltaVel;

        //         int k = OpenMPUtils::ThisThread();

        //         ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
        //         ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

        //         for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {
        //             noalias(DeltaVel) = (itNode)->FastGetSolutionStepValue(VELOCITY) - (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

        //             array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
        //             array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);

        //             UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

        //             if (mMeshVelocity == 2)//Lagrangian
        //             {
        //                 if((itNode)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) < 1e-15)
        //                 {
        //                     array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);
        //                     array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
        //                     array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

        //                     noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
        //                     UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
        //                 }
        //                 else
        //                 {
        //                     noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = ZeroVector(3);
        //                     noalias(itNode->FastGetSolutionStepValue(DISPLACEMENT)) = ZeroVector(3);
        //                 }
        //             }
        //         }
        //     }

        //     KRATOS_CATCH("")

        // }

        /**
        * Performing the prediction of the solution
        * @param rModelPart: The model of the problem to solve
        * @param rDofSet set of all primary variables
        * @param rDx: Incremental update of primary variables
        */

        void Predict(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemVectorType& rDx) override
        {
            KRATOS_TRY;

            this->PredictVariables(rModelPart);

            //TODO: predict mesh_acceleration, mesh_velocity, mesh_displacement

            this->MoveMesh(rModelPart);

            KRATOS_CATCH( "" );
        }

        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold
//         void Predict(ModelPart& rModelPart,
//                              DofsArrayType& rDofSet,
//                              TSystemMatrixType& A,
//                              TSystemVectorType& Dv,
//                              TSystemVectorType& b) override
//         {
//             // if (rModelPart.GetCommunicator().MyPID() == 0)
//             //     std::cout << "prediction" << std::endl;

//             int NumThreads = OpenMPUtils::GetNumThreads();
//             OpenMPUtils::PartitionVector NodePartition;
//             OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

//             #pragma omp parallel
//             {
//                 //array_1d<double, 3 > DeltaDisp;

//                 int k = OpenMPUtils::ThisThread();

//                 ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
//                 ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

//                 for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {
//                     array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);
//                     double& OldPressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 1);

//                     //predicting velocity
//                     //ATTENTION::: the prediction is performed only on free nodes
//                     array_1d<double, 3 > & CurrentVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY);
//                     double& CurrentPressure = (itNode)->FastGetSolutionStepValue(PRESSURE);

//                     if ((itNode->pGetDof(VELOCITY_X))->IsFree())
//                         (CurrentVelocity[0]) = OldVelocity[0];
//                     if (itNode->pGetDof(VELOCITY_Y)->IsFree())
//                         (CurrentVelocity[1]) = OldVelocity[1];
//                     if (itNode->HasDofFor(VELOCITY_Z))
//                         if (itNode->pGetDof(VELOCITY_Z)->IsFree())
//                             (CurrentVelocity[2]) = OldVelocity[2];

//                     if (itNode->pGetDof(PRESSURE)->IsFree())
//                         CurrentPressure = OldPressure;

//                     // updating time derivatives ::: please note that displacements and
//                     // their time derivatives can not be consistently fixed separately
//                     array_1d<double, 3 > DeltaVel;
//                     noalias(DeltaVel) = CurrentVelocity - OldVelocity;
//                     array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
//                     array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);

//                     UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

//                     if (mMeshVelocity == 2) //Lagrangian
//                     {
//                         array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
//                         array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);

//                   if((itNode)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) < 1e-15)
// 			{
// 			    noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
// 			    UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
// 			}
// 			else
// 			{
// 			  itNode->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
// 			  itNode->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
// 			  itNode->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
// 			  itNode->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
// 			}
//                     }
//                 }
//             }

// //              if (rModelPart.GetCommunicator().MyPID() == 0)
// //                  std::cout << "end of prediction" << std::endl;

//         }

        /**
        * @brief This function is designed to move the mesh
        * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
        */
        void MoveMesh(ModelPart& rModelPart) override
        {
            KRATOS_TRY

            if( this->mOptions.Is(LocalFlagType::MOVE_MESH) ){

                if (rModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false)
                {
                    KRATOS_ERROR << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the Model Part. Add DISPLACEMENT to the list of variables" << std::endl;
                }

                bool DisplacementIntegration = false;
                for(typename IntegrationMethodsVectorType::iterator it=mTimeVectorIntegrationMethods.begin();
                    it!=mTimeVectorIntegrationMethods.end(); ++it)
                {
                    if( "DISPLACEMENT" == (*it)->GetVariableName() ){
                        DisplacementIntegration = true;
                        break;
                    }
                }

                if(DisplacementIntegration == true){

                    // Update mesh positions : node coordinates
                    const int nnodes = rModelPart.NumberOfNodes();
                    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();

                #pragma omp parallel for
                    for(int i = 0; i<nnodes; i++)
                    {
                        ModelPart::NodesContainerType::iterator it_node = it_begin + i;

                        // TODO: review
                        noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
                        if( it_node->IsNot(SLIP) ){
                            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
                        } else {
                            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
                        }

                    }
                }
            }

            KRATOS_CATCH("")
        }

        //***************************************************************************

        /** this function is designed to be called in the builder and solver
        to introduce
        the selected time integration scheme. It "asks" the matrix needed to
        the element and
        performs the operations needed to introduce the seected time
        integration scheme.

          this function calculates at the same time the contribution to the
        LHS and to the RHS
          of the system
         */
        void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                          LocalSystemMatrixType& rLHS_Contribution,
                                          LocalSystemVectorType& rRHS_Contribution,
                                          Element::EquationIdVectorType& EquationId,
                                          ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            int thread = OpenMPUtils::ThisThread();

            //Initializing the non linear iteration for the current element
            // (pCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo); //TODO: delete

            (pCurrentElement) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

            (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

            // TODO: this will be generally false
            if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
                (pCurrentElement) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
                (pCurrentElement) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

                AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);
                AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
            }
            else{
                (pCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);
                (pCurrentElement)->CalculateLocalVelocityContribution(mMatrix.D[thread], rRHS_Contribution, rCurrentProcessInfo);

                AddDynamicsToLHS (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
                AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
            }

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,pCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,pCurrentElement->GetGeometry());

            KRATOS_CATCH("")
        }

        void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                        LocalSystemVectorType& rRHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,
                                        ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            int thread = OpenMPUtils::ThisThread();

            //Initializing the non linear iteration for the current element
            // (pCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo); //TODO: delete

            // Basic operations for the element considered
            (pCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

            (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

            if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

                (pCurrentElement) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);
                (pCurrentElement) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

                AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

            }
            else{

                (pCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);
                (pCurrentElement)->CalculateLocalVelocityContribution(mMatrix.D[thread], rRHS_Contribution, rCurrentProcessInfo);

                AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
            }

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(rRHS_Contribution,pCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(rRHS_Contribution,pCurrentElement->GetGeometry());

            KRATOS_CATCH("")
        }

        /** functions totally analogous to the precedent but applied to
        the "condition" objects
         */
        void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                            LocalSystemMatrixType& rLHS_Contribution,
                                                            LocalSystemVectorType& rRHS_Contribution,
                                                            Element::EquationIdVectorType& EquationId,
                                                            ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            int thread = OpenMPUtils::ThisThread();

            // (pCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo); //TODO: delete

            // Basic operations for the condition considered
            (pCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

            (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

            if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
                (pCurrentCondition) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);
                (pCurrentCondition) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

                AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);
                AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
            }
            else{
                (pCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);
                (pCurrentCondition)->CalculateLocalVelocityContribution(mMatrix.D[thread], rRHS_Contribution, rCurrentProcessInfo);

                AddDynamicsToLHS  (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
                AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
            }

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,pCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,pCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                                          LocalSystemVectorType& rRHS_Contribution,
                                                          Element::EquationIdVectorType& EquationId,
                                                          ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            int thread = OpenMPUtils::ThisThread();

            // (pCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo); //TODO: delete

            // Basic operations for the condition considered
            (pCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

            (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

            if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){
                (pCurrentCondition) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);
                (pCurrentCondition) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

                AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);
            }
            else{
                (pCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);
                (pCurrentCondition)->CalculateLocalVelocityContribution(mMatrix.D[thread], rRHS_Contribution,rCurrentProcessInfo);

                AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
            }

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(rRHS_Contribution,pCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(rRHS_Contribution,pCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }


        //TODO: is this necessary in PFEM?
        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
            Element::EquationIdVectorType EquationId;
            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)

            #pragma omp parallel for
            for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
            {
                auto itNode = rModelPart.NodesBegin() + k;
                (itNode->FastGetSolutionStepValue(REACTION)).clear();

                // calculating relaxed acceleration
                const array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
                const array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                const array_1d<double, 3> relaxed_acceleration = (1 - mAlphaBossak) * CurrentAcceleration
                                                                    + mAlphaBossak * OldAcceleration;
                (itNode)->SetValue(RELAXED_ACCELERATION, relaxed_acceleration);
            }

            //for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)

            #pragma omp parallel for firstprivate(EquationId,RHS_Contribution,LHS_Contribution)
            for(int k = 0; k<static_cast<int>(rModelPart.Elements().size()); k++)
            {
                auto itElem = rModelPart.Elements().ptr_begin()+k;
                int thread_id = OpenMPUtils::ThisThread();

                (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
                //KRATOS_WATCH(LHS_Contribution);
                //basic operations for the element considered
                (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

                //std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
                (*itElem)->CalculateMassMatrix(mMass[thread_id], CurrentProcessInfo);
                (*itElem)->CalculateLocalVelocityContribution(mDamp[thread_id], RHS_Contribution, CurrentProcessInfo);

                (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);

                //adding the dynamic contributions (statics is already included)
                AddDynamicsToLHS(LHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);
                AddDynamicsToRHS((*itElem), RHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);

                GeometryType& rGeom = (*itElem)->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int Dimension = rGeom.WorkingSpaceDimension();

                unsigned int index = 0;
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    auto& reaction = rGeom[i].FastGetSolutionStepValue(REACTION);

                    double& target_value0 = reaction[0];
                    const double& origin_value0 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value0 -= origin_value0;

                    double& target_value1 = reaction[1];
                    const double& origin_value1 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value1 -= origin_value1;

                    if (Dimension == 3)
                    {
                      double& target_value2 = reaction[2];
                      const double& origin_value2 = RHS_Contribution[index++];
                      #pragma omp atomic
                      target_value2 -= origin_value2;
                    }
            //        rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
             //          rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
            //        if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
                    index++; // skip pressure dof
                }
            }

            rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

            // Base scheme calls FinalizeSolutionStep method of elements and conditions
            Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
        }

        //************************************************************************************************
        //************************************************************************************************

        /*@} */
        /**@name Operations */
        /*@{ */


        /*@} */
        /**@name Access */
        /*@{ */


        /*@} */
        /**@name Inquiry */
        /*@{ */


        /*@} */
        /**@name Friends */
        /*@{ */


        /*@} */

    protected:
        /**@name Protected static Member Variables */
        /*@{ */


        /*@} */
        /**@name Protected member Variables */
        /*@{ */

        /*@} */
        /**@name Protected Operators*/
        /*@{ */

        /*@} */
        /**@name Protected Operations*/
        /*@{ */


        /*@} */
        /**@name Protected  Access */
        /*@{ */


        /*@} */
        /**@name Protected Inquiry */
        /*@{ */


        /*@} */
        /**@name Protected LifeCycle */
        /*@{ */



        /*@} */

    private:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */

        CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;

        // typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater(); //TODO: delete

        /*@} */
        /**@name Private Operators*/
        /*@{ */

        /*@} */
        /**@name Private Operations*/
        /*@{ */


        /*@} */
        /**@name Private  Access */
        /*@{ */


        /*@} */
        /**@name Private Inquiry */
        /*@{ */


        /*@} */
        /**@name Un accessible methods */
        /*@{ */


        /*@} */

    }; /* Class Scheme */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_ALE_SOLUTION_SCHEME_H_INCLUDED  defined */
