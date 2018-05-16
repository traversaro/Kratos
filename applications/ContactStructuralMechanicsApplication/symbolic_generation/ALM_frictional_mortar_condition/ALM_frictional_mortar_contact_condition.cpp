// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom ) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation> >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::~AugmentedLagrangianMethodFrictionalMortarContactCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Initialize( )
{
    KRATOS_TRY;

    BaseType::Initialize();

    // We initailize the mortar operators
    mCurrentMortarOperators.Initialize();
    mPreviousMortarOperators.Initialize();
    mPreviousMortarOperatorsInitialized = false;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    // We "save" the mortar operator in case not initialized
    if (mPreviousMortarOperatorsInitialized == false) {
        ComputeStandardMortarOperators(mPreviousMortarOperators, rCurrentProcessInfo);
        mPreviousMortarOperatorsInitialized = true;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    // We compute the current standard mortar operators
    ComputeStandardMortarOperators(mCurrentMortarOperators, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    // We "save" the mortar operator for the next step
    ComputeStandardMortarOperators(mPreviousMortarOperators, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // We compute the current standard mortar operators
    ComputeStandardMortarOperators(mCurrentMortarOperators, rCurrentProcessInfo);
    if (mPreviousMortarOperatorsInitialized == false) {
        ComputeStandardMortarOperators(mPreviousMortarOperators, rCurrentProcessInfo);
        mPreviousMortarOperatorsInitialized = true;
    }

    // The slave geometry
    GeometryType& slave_geometry = this->GetGeometry();
    const array_1d<double, 3>& normal_slave = this->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables rVariables;

    // Create the current contact data
    DerivativeDataType rDerivativeData;
    rDerivativeData.Initialize(slave_geometry, rCurrentProcessInfo);

    // Create the mortar operators
    MortarConditionMatrices rThisMortarConditionMatrices;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo[DISTANCE_THRESHOLD];
    IntegrationUtility integration_utility = IntegrationUtility (BaseType::mIntegrationOrder, distance_threshold);

    // If we consider the normal variation
    const NormalDerivativesComputation consider_normal_variation = static_cast<NormalDerivativesComputation>(rCurrentProcessInfo[CONSIDER_NORMAL_VARIATION]);

    // The master geometry
    GeometryType& master_geometry = this->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& normal_master = this->GetValue(PAIRED_NORMAL);

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, normal_slave, master_geometry, normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(slave_geometry, conditions_points_slave, integration_area);

    const double geometry_area = slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-3 * geometry_area)) {
        IntegrationMethod this_integration_method = this->GetIntegrationMethod();

        // Initialize general variables for the current master element
        rVariables.Initialize();

        // Update slave element info
        rDerivativeData.UpdateMasterPair(master_geometry, rCurrentProcessInfo);

        // Initialize the mortar operators
        rThisMortarConditionMatrices.Initialize();

        const bool dual_LM = DerivativesUtilitiesType::CalculateAeAndDeltaAe(slave_geometry, normal_slave, master_geometry, rDerivativeData, rVariables, consider_normal_variation, conditions_points_slave, this_integration_method, this->GetAxisymmetricCoefficient(rVariables));

        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            array_1d<BelongType, TDim> belong_array;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                PointType global_point;
                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array[i_node] = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

            if (bad_shape == false) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates
                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    this->CalculateKinematics( rVariables, rDerivativeData, normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double integration_weight = integration_points_slave[point_number].Weight() * this->GetAxisymmetricCoefficient(rVariables);

                    rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);
                }
            }
        }

        // Setting the weighted gap
        // Mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rThisMortarConditionMatrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperator = rThisMortarConditionMatrices.MOperator;

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry);
        const BoundedMatrix<double, TNumNodes, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry);

        const BoundedMatrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(DOperator, x1) - prod(MOperator, x2);

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& normal = slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
            const array_1d<double, TDim> aux_array = row(D_x1_M_x2, i_node);

            double& weighted_gap = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);

            #pragma omp atomic
            weighted_gap += inner_prod(aux_array, - subrange(normal, 0, TDim));
        }

        // Setting the weighted slip
        // The increment of time
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        // Auxiliar slip derivative
        BoundedMatrix<double, TNumNodes, TDim> slip_time_derivative;

        // Delta mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TDim> x1old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry, false, 1);
        const BoundedMatrix<double, TNumNodes, TDim> x2old = MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry, false, 1);
        const BoundedMatrix<double, TNumNodes, TNumNodes> DeltaDOperator = mCurrentMortarOperators.DOperator - mPreviousMortarOperators.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodes> DeltaMOperator = mCurrentMortarOperators.MOperator - mPreviousMortarOperators.MOperator;

//         const double tolerance = std::numeric_limits<double>::epsilon();
//         if (norm_frobenius(DeltaDOperator) > tolerance && norm_frobenius(DeltaMOperator) > tolerance) {  // Frame indifferent
            const BoundedMatrix<double, TNumNodes, TDim> delta_D_x1_delta_M_x2 = prod(DeltaDOperator, x1old) - prod(DeltaMOperator, x2old);
//             const BoundedMatrix<double, TNumNodes, TDim> delta_D_x1_delta_M_x2 = prod(DeltaDOperator, x1) - prod(DeltaMOperator, x2);

            // The estimation of the slip time derivative
            slip_time_derivative = - delta_D_x1_delta_M_x2/delta_time;

//             const BoundedMatrix<double, TNumNodes, TDim> u1 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(slave_geometry, DISPLACEMENT, 0) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(slave_geometry, DISPLACEMENT, 1);
//             const BoundedMatrix<double, TNumNodes, TDim> u2 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(master_geometry, DISPLACEMENT, 0) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(master_geometry, DISPLACEMENT, 1);
//             const BoundedMatrix<double, TNumNodes, TDim> delta_D_u1_delta_M_u2 = prod(DeltaDOperator, u1) - prod(DeltaMOperator, u2);
//
//             // The estimation of the slip time derivative
//             slip_time_derivative = delta_D_u1_delta_M_u2/delta_time;
//
//         } else { // Standard definition
//             // Old coordinates
//             const BoundedMatrix<double, TNumNodes, TDim> delta_x1 = x1- MortarUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry, false, 1);
//             const BoundedMatrix<double, TNumNodes, TDim> delta_x2 = x2- MortarUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry, false, 1);
//
//             const BoundedMatrix<double, TNumNodes, TDim> D_delta_x1_M_delta_x2 = prod(DOperator, x1) - prod(MOperator, x2);
//             const BoundedMatrix<double, TNumNodes, TDim> D_delta_x1_M_delta_x2 = prod(DOperator, delta_x1) - prod(MOperator, delta_x2);
// //             const BoundedMatrix<double, TNumNodes, TDim> D_delta_x1_M_delta_x2 = prod(mCurrentMortarOperators.DOperator, delta_x1) - prod(mCurrentMortarOperators.MOperator, delta_x2);
//
// //             slip_time_derivative += D_delta_x1_M_delta_x2/delta_time;
// // //             const BoundedMatrix<double, TNumNodes, TDim> D_delta_x1_M_delta_x2 = prod(DOperator, delta_x1) - prod(MOperator, delta_x2);
// //
//             slip_time_derivative = D_delta_x1_M_delta_x2/delta_time;
// //         }

//         const BoundedMatrix<double, TNumNodes, TDim> slip_time_derivative = D_delta_x1_M_delta_x2/delta_time + delta_D_x1_delta_M_x2/delta_time;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            // We get the normal
            const array_1d<double, TDim>& normal = subrange(slave_geometry[i_node].FastGetSolutionStepValue(NORMAL), 0, TDim);

            // We compute the slip
            const array_1d<double, TDim>& slip_time_derivative_node = row(slip_time_derivative, i_node);
            const array_1d<double, TDim> slip_node = - delta_time * (slip_time_derivative_node - inner_prod(normal, slip_time_derivative_node) * normal);

            // The weighted slip
            array_1d<double, 3>& weighted_slip = slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP);

            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                #pragma omp atomic
                weighted_slip[i_dim] += slip_node[i_dim];
            }
        }

        // We reset the flag
        this->Set(ISOLATED, false);
    } else {
        // We set the flag
        this->Set(ISOLATED, true);
    }

    KRATOS_CATCH( "" );
}

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/


/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2, false>::CalculateLocalLHS(
    Matrix& rLocalLHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Initialize
    for (std::size_t i = 0; i < 12; ++i)
        for (std::size_t j = 0; j < 12; ++j)
            rLocalLHS(i, j) = 0.0;

    // The geometry of the condition
    GeometryType& geometry = this->GetGeometry();

    // Initialize values
    const BoundedMatrix<double, 2, 2>& u1 = rDerivativeData.u1;
    const BoundedMatrix<double, 2, 2>& u1old = rDerivativeData.u1old;
    const BoundedMatrix<double, 2, 2>& u2 = rDerivativeData.u2;
    const BoundedMatrix<double, 2, 2>& u2old = rDerivativeData.u2old;
    const BoundedMatrix<double, 2, 2>& X1 = rDerivativeData.X1;
    const BoundedMatrix<double, 2, 2>& X2 = rDerivativeData.X2;
    
    const BoundedMatrix<double, 2, 2> LM = MortarUtilities::GetVariableMatrix<2,2>(geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);
    
    // The normal and tangent vectors
    const BoundedMatrix<double, 2, 2>& NormalSlave = rDerivativeData.NormalSlave;
    const BoundedMatrix<double, 2, 2> TangentSlave = MortarUtilities::ComputeTangentMatrix<2,2>(geometry);
    const BoundedMatrix<double, 2, 2> TangentSlaveXi = MortarUtilities::GetVariableMatrix<2,2>(geometry, TANGENT_XI, 0);
    const BoundedMatrix<double, 2, 2> TangentSlaveEta = MortarUtilities::GetVariableMatrix<2,2>(geometry, TANGENT_ETA, 0);

    // The ALM parameters
    const array_1d<double, 2> DynamicFactor = MortarUtilities::GetVariableVector<2>(geometry, DYNAMIC_FACTOR);
    const double ScaleFactor = rDerivativeData.ScaleFactor;
    const array_1d<double, 2>& PenaltyParameter = rDerivativeData.PenaltyParameter;
    const double TangentFactor = rDerivativeData.TangentFactor;
    
    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;
    const BoundedMatrix<double, 2, 2>& StandardMOperator = mCurrentMortarOperators.MOperator;
    const BoundedMatrix<double, 2, 2>& StandardDOperator = mCurrentMortarOperators.DOperator;
    const BoundedMatrix<double, 2, 2>& StandardMOperatorold = mPreviousMortarOperators.MOperator;
    const BoundedMatrix<double, 2, 2>& StandardDOperatorold = mPreviousMortarOperators.DOperator;

    // Mortar operators derivatives
    const array_1d<BoundedMatrix<double, 2, 2>, 8>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;
    const array_1d<BoundedMatrix<double, 2, 2>, 8>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;

    // We get the friction coefficient
    const array_1d<double, 2> mu = GetFrictionCoefficient();

//    // The delta time
//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    
    // NODE 0
    if (geometry[0].IsNot(ACTIVE)) { // INACTIVE
        const double clhs0 =     std::pow(ScaleFactor, 2)/PenaltyParameter[0];
    
        rLocalLHS(8,8)+=NormalSlave(0,0)*clhs0;
        rLocalLHS(8,9)+=NormalSlave(0,1)*clhs0;
        rLocalLHS(9,8)+=TangentSlaveXi(0,0)*clhs0;
        rLocalLHS(9,9)+=TangentSlaveXi(0,1)*clhs0;
    } else if (geometry[0].Is(SLIP)) { // ACTIVE-SLIP
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     X1(0,1) + u1(0,1);
        const double clhs2 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs3 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs4 =     X1(1,1) + u1(1,1);
        const double clhs5 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X2(0,1) + u2(0,1);
        const double clhs8 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs9 =     X2(1,1) + u2(1,1);
        const double clhs10 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs11 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     NormalSlave(0,1)*(clhs1*clhs3 - clhs11*clhs9 + clhs4*clhs6 - clhs7*clhs8);
        const double clhs13 =     -clhs0;
        const double clhs14 =     X1(0,0) + u1(0,0);
        const double clhs15 =     clhs14*clhs3;
        const double clhs16 =     X1(1,0) + u1(1,0);
        const double clhs17 =     clhs16*clhs6;
        const double clhs18 =     X2(0,0) + u2(0,0);
        const double clhs19 =     clhs18*clhs8;
        const double clhs20 =     X2(1,0) + u2(1,0);
        const double clhs21 =     clhs11*clhs20;
        const double clhs22 =     NormalSlave(0,0)*(clhs13 + clhs15 + clhs17 - clhs19 - clhs21) + clhs12;
        const double clhs23 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs22;
        const double clhs24 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double clhs25 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double clhs26 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double clhs27 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double clhs28 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(clhs24*(X1(0,0) + u1old(0,0)) + clhs25*(X1(1,0) + u1old(1,0)) - clhs26*(X2(0,0) + u2old(0,0)) - clhs27*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(0,1)*(clhs24*(X1(0,1) + u1old(0,1)) + clhs25*(X1(1,1) + u1old(1,1)) - clhs26*(X2(0,1) + u2old(0,1)) - clhs27*(X2(1,1) + u2old(1,1))));
        const double clhs29 =     PenaltyParameter[0]*(NormalSlave(0,0)*(-clhs0*clhs18 - clhs10*clhs20 + clhs14*clhs2 + clhs16*clhs5) + NormalSlave(0,1)*(-clhs0*clhs7 + clhs1*clhs2 - clhs10*clhs9 + clhs4*clhs5));
        const double clhs30 =     LM(0,0)*ScaleFactor - NormalSlave(0,0)*clhs29 - TangentSlaveXi(0,0)*clhs28;
        const double clhs31 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs32 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs33 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs34 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs35 =     NormalSlave(0,0)*(clhs14*clhs31 + clhs16*clhs32 - clhs18*clhs33 - clhs20*clhs34);
        const double clhs36 =     clhs1*clhs31;
        const double clhs37 =     clhs32*clhs4;
        const double clhs38 =     clhs33*clhs7;
        const double clhs39 =     clhs34*clhs9;
        const double clhs40 =     NormalSlave(0,1)*(clhs13 + clhs36 + clhs37 - clhs38 - clhs39) + clhs35;
        const double clhs41 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs40;
        const double clhs42 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs43 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs44 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs45 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs46 =     NormalSlave(0,1)*(clhs1*clhs42 + clhs4*clhs43 - clhs44*clhs7 - clhs45*clhs9);
        const double clhs47 =     -clhs10;
        const double clhs48 =     clhs14*clhs42;
        const double clhs49 =     clhs16*clhs43;
        const double clhs50 =     clhs18*clhs44;
        const double clhs51 =     clhs20*clhs45;
        const double clhs52 =     NormalSlave(0,0)*(clhs47 + clhs48 + clhs49 - clhs50 - clhs51) + clhs46;
        const double clhs53 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs52;
        const double clhs54 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs55 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs56 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs57 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs58 =     NormalSlave(0,0)*(clhs14*clhs54 + clhs16*clhs55 - clhs18*clhs56 - clhs20*clhs57);
        const double clhs59 =     clhs1*clhs54;
        const double clhs60 =     clhs4*clhs55;
        const double clhs61 =     clhs56*clhs7;
        const double clhs62 =     clhs57*clhs9;
        const double clhs63 =     NormalSlave(0,1)*(clhs47 + clhs59 + clhs60 - clhs61 - clhs62) + clhs58;
        const double clhs64 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs63;
        const double clhs65 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs66 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs67 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs68 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs69 =     NormalSlave(0,0)*(clhs14*clhs65 + clhs16*clhs66 - clhs18*clhs67 + clhs2 - clhs20*clhs68) + NormalSlave(0,1)*(clhs1*clhs65 + clhs4*clhs66 - clhs67*clhs7 - clhs68*clhs9);
        const double clhs70 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs69;
        const double clhs71 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs72 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs73 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs74 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs75 =     NormalSlave(0,0)*(clhs14*clhs71 + clhs16*clhs72 - clhs18*clhs73 - clhs20*clhs74) + NormalSlave(0,1)*(clhs1*clhs71 + clhs2 + clhs4*clhs72 - clhs7*clhs73 - clhs74*clhs9);
        const double clhs76 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs75;
        const double clhs77 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs78 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs79 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs80 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs81 =     NormalSlave(0,0)*(clhs14*clhs77 + clhs16*clhs78 - clhs18*clhs79 - clhs20*clhs80 + clhs5) + NormalSlave(0,1)*(clhs1*clhs77 + clhs4*clhs78 - clhs7*clhs79 - clhs80*clhs9);
        const double clhs82 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs81;
        const double clhs83 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs84 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs85 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs86 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs87 =     NormalSlave(0,0)*(clhs14*clhs83 + clhs16*clhs84 - clhs18*clhs85 - clhs20*clhs86) + NormalSlave(0,1)*(clhs1*clhs83 + clhs4*clhs84 + clhs5 - clhs7*clhs85 - clhs86*clhs9);
        const double clhs88 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs87;
        const double clhs89 =     DynamicFactor[0]*ScaleFactor;
        const double clhs90 =     clhs0*clhs89;
        const double clhs91 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs22;
        const double clhs92 =     LM(0,1)*ScaleFactor - NormalSlave(0,1)*clhs29 - TangentSlaveXi(0,1)*clhs28;
        const double clhs93 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs40;
        const double clhs94 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs52;
        const double clhs95 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs63;
        const double clhs96 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs69;
        const double clhs97 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs75;
        const double clhs98 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs81;
        const double clhs99 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs87;
        const double clhs100 =     clhs10*clhs89;
        const double clhs101 =     -clhs2*clhs89;
        const double clhs102 =     -clhs5*clhs89;
    
        rLocalLHS(0,0)+=DynamicFactor[0]*(-clhs0*clhs23 + clhs30*clhs8);
        rLocalLHS(0,1)+=DynamicFactor[0]*(-clhs0*clhs41 + clhs30*clhs33);
        rLocalLHS(0,2)+=DynamicFactor[0]*(-clhs0*clhs53 + clhs30*clhs44);
        rLocalLHS(0,3)+=DynamicFactor[0]*(-clhs0*clhs64 + clhs30*clhs56);
        rLocalLHS(0,4)+=DynamicFactor[0]*(-clhs0*clhs70 + clhs30*clhs67);
        rLocalLHS(0,5)+=DynamicFactor[0]*(-clhs0*clhs76 + clhs30*clhs73);
        rLocalLHS(0,6)+=DynamicFactor[0]*(-clhs0*clhs82 + clhs30*clhs79);
        rLocalLHS(0,7)+=DynamicFactor[0]*(-clhs0*clhs88 + clhs30*clhs85);
        rLocalLHS(0,8)+=clhs90;
        rLocalLHS(1,0)+=DynamicFactor[0]*(-clhs0*clhs91 + clhs8*clhs92);
        rLocalLHS(1,1)+=DynamicFactor[0]*(-clhs0*clhs93 + clhs33*clhs92);
        rLocalLHS(1,2)+=DynamicFactor[0]*(-clhs0*clhs94 + clhs44*clhs92);
        rLocalLHS(1,3)+=DynamicFactor[0]*(-clhs0*clhs95 + clhs56*clhs92);
        rLocalLHS(1,4)+=DynamicFactor[0]*(-clhs0*clhs96 + clhs67*clhs92);
        rLocalLHS(1,5)+=DynamicFactor[0]*(-clhs0*clhs97 + clhs73*clhs92);
        rLocalLHS(1,6)+=DynamicFactor[0]*(-clhs0*clhs98 + clhs79*clhs92);
        rLocalLHS(1,7)+=DynamicFactor[0]*(-clhs0*clhs99 + clhs85*clhs92);
        rLocalLHS(1,9)+=clhs90;
        rLocalLHS(2,0)+=DynamicFactor[0]*(-clhs10*clhs23 + clhs11*clhs30);
        rLocalLHS(2,1)+=DynamicFactor[0]*(-clhs10*clhs41 + clhs30*clhs34);
        rLocalLHS(2,2)+=DynamicFactor[0]*(-clhs10*clhs53 + clhs30*clhs45);
        rLocalLHS(2,3)+=DynamicFactor[0]*(-clhs10*clhs64 + clhs30*clhs57);
        rLocalLHS(2,4)+=DynamicFactor[0]*(-clhs10*clhs70 + clhs30*clhs68);
        rLocalLHS(2,5)+=DynamicFactor[0]*(-clhs10*clhs76 + clhs30*clhs74);
        rLocalLHS(2,6)+=DynamicFactor[0]*(-clhs10*clhs82 + clhs30*clhs80);
        rLocalLHS(2,7)+=DynamicFactor[0]*(-clhs10*clhs88 + clhs30*clhs86);
        rLocalLHS(2,8)+=clhs100;
        rLocalLHS(3,0)+=DynamicFactor[0]*(-clhs10*clhs91 + clhs11*clhs92);
        rLocalLHS(3,1)+=DynamicFactor[0]*(-clhs10*clhs93 + clhs34*clhs92);
        rLocalLHS(3,2)+=DynamicFactor[0]*(-clhs10*clhs94 + clhs45*clhs92);
        rLocalLHS(3,3)+=DynamicFactor[0]*(-clhs10*clhs95 + clhs57*clhs92);
        rLocalLHS(3,4)+=DynamicFactor[0]*(-clhs10*clhs96 + clhs68*clhs92);
        rLocalLHS(3,5)+=DynamicFactor[0]*(-clhs10*clhs97 + clhs74*clhs92);
        rLocalLHS(3,6)+=DynamicFactor[0]*(-clhs10*clhs98 + clhs80*clhs92);
        rLocalLHS(3,7)+=DynamicFactor[0]*(-clhs10*clhs99 + clhs86*clhs92);
        rLocalLHS(3,9)+=clhs100;
        rLocalLHS(4,0)+=DynamicFactor[0]*(clhs2*clhs23 - clhs3*clhs30);
        rLocalLHS(4,1)+=DynamicFactor[0]*(clhs2*clhs41 - clhs30*clhs31);
        rLocalLHS(4,2)+=DynamicFactor[0]*(clhs2*clhs53 - clhs30*clhs42);
        rLocalLHS(4,3)+=DynamicFactor[0]*(clhs2*clhs64 - clhs30*clhs54);
        rLocalLHS(4,4)+=DynamicFactor[0]*(clhs2*clhs70 - clhs30*clhs65);
        rLocalLHS(4,5)+=DynamicFactor[0]*(clhs2*clhs76 - clhs30*clhs71);
        rLocalLHS(4,6)+=DynamicFactor[0]*(clhs2*clhs82 - clhs30*clhs77);
        rLocalLHS(4,7)+=DynamicFactor[0]*(clhs2*clhs88 - clhs30*clhs83);
        rLocalLHS(4,8)+=clhs101;
        rLocalLHS(5,0)+=DynamicFactor[0]*(clhs2*clhs91 - clhs3*clhs92);
        rLocalLHS(5,1)+=DynamicFactor[0]*(clhs2*clhs93 - clhs31*clhs92);
        rLocalLHS(5,2)+=DynamicFactor[0]*(clhs2*clhs94 - clhs42*clhs92);
        rLocalLHS(5,3)+=DynamicFactor[0]*(clhs2*clhs95 - clhs54*clhs92);
        rLocalLHS(5,4)+=DynamicFactor[0]*(clhs2*clhs96 - clhs65*clhs92);
        rLocalLHS(5,5)+=DynamicFactor[0]*(clhs2*clhs97 - clhs71*clhs92);
        rLocalLHS(5,6)+=DynamicFactor[0]*(clhs2*clhs98 - clhs77*clhs92);
        rLocalLHS(5,7)+=DynamicFactor[0]*(clhs2*clhs99 - clhs83*clhs92);
        rLocalLHS(5,9)+=clhs101;
        rLocalLHS(6,0)+=DynamicFactor[0]*(clhs23*clhs5 - clhs30*clhs6);
        rLocalLHS(6,1)+=DynamicFactor[0]*(-clhs30*clhs32 + clhs41*clhs5);
        rLocalLHS(6,2)+=DynamicFactor[0]*(-clhs30*clhs43 + clhs5*clhs53);
        rLocalLHS(6,3)+=DynamicFactor[0]*(-clhs30*clhs55 + clhs5*clhs64);
        rLocalLHS(6,4)+=DynamicFactor[0]*(-clhs30*clhs66 + clhs5*clhs70);
        rLocalLHS(6,5)+=DynamicFactor[0]*(-clhs30*clhs72 + clhs5*clhs76);
        rLocalLHS(6,6)+=DynamicFactor[0]*(-clhs30*clhs78 + clhs5*clhs82);
        rLocalLHS(6,7)+=DynamicFactor[0]*(-clhs30*clhs84 + clhs5*clhs88);
        rLocalLHS(6,8)+=clhs102;
        rLocalLHS(7,0)+=DynamicFactor[0]*(clhs5*clhs91 - clhs6*clhs92);
        rLocalLHS(7,1)+=DynamicFactor[0]*(-clhs32*clhs92 + clhs5*clhs93);
        rLocalLHS(7,2)+=DynamicFactor[0]*(-clhs43*clhs92 + clhs5*clhs94);
        rLocalLHS(7,3)+=DynamicFactor[0]*(clhs5*clhs95 - clhs55*clhs92);
        rLocalLHS(7,4)+=DynamicFactor[0]*(clhs5*clhs96 - clhs66*clhs92);
        rLocalLHS(7,5)+=DynamicFactor[0]*(clhs5*clhs97 - clhs72*clhs92);
        rLocalLHS(7,6)+=DynamicFactor[0]*(clhs5*clhs98 - clhs78*clhs92);
        rLocalLHS(7,7)+=DynamicFactor[0]*(clhs5*clhs99 - clhs84*clhs92);
        rLocalLHS(7,9)+=clhs102;
        rLocalLHS(8,0)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs0 - clhs15 - clhs17 + clhs19 + clhs21) + clhs12);
        rLocalLHS(8,1)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs0 - clhs36 - clhs37 + clhs38 + clhs39) + clhs35);
        rLocalLHS(8,2)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs10 - clhs48 - clhs49 + clhs50 + clhs51) + clhs46);
        rLocalLHS(8,3)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs10 - clhs59 - clhs60 + clhs61 + clhs62) + clhs58);
        rLocalLHS(8,4)+=-ScaleFactor*clhs69;
        rLocalLHS(8,5)+=-ScaleFactor*clhs75;
        rLocalLHS(8,6)+=-ScaleFactor*clhs81;
        rLocalLHS(8,7)+=-ScaleFactor*clhs87;
        rLocalLHS(9,8)+=TangentSlaveXi(0,0);
        rLocalLHS(9,9)+=TangentSlaveXi(0,1);
    } else { // ACTIVE-STICK
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     X1(0,1) + u1(0,1);
        const double clhs2 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs3 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs4 =     X1(1,1) + u1(1,1);
        const double clhs5 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X2(0,1) + u2(0,1);
        const double clhs8 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs9 =     X2(1,1) + u2(1,1);
        const double clhs10 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs11 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     NormalSlave(0,1)*(clhs1*clhs3 - clhs11*clhs9 + clhs4*clhs6 - clhs7*clhs8);
        const double clhs13 =     -clhs0;
        const double clhs14 =     X1(0,0) + u1(0,0);
        const double clhs15 =     clhs14*clhs3;
        const double clhs16 =     X1(1,0) + u1(1,0);
        const double clhs17 =     clhs16*clhs6;
        const double clhs18 =     X2(0,0) + u2(0,0);
        const double clhs19 =     clhs18*clhs8;
        const double clhs20 =     X2(1,0) + u2(1,0);
        const double clhs21 =     clhs11*clhs20;
        const double clhs22 =     NormalSlave(0,0)*(clhs13 + clhs15 + clhs17 - clhs19 - clhs21) + clhs12;
        const double clhs23 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs22;
        const double clhs24 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double clhs25 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double clhs26 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double clhs27 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double clhs28 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(clhs24*(X1(0,0) + u1old(0,0)) + clhs25*(X1(1,0) + u1old(1,0)) - clhs26*(X2(0,0) + u2old(0,0)) - clhs27*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(0,1)*(clhs24*(X1(0,1) + u1old(0,1)) + clhs25*(X1(1,1) + u1old(1,1)) - clhs26*(X2(0,1) + u2old(0,1)) - clhs27*(X2(1,1) + u2old(1,1))));
        const double clhs29 =     PenaltyParameter[0]*(NormalSlave(0,0)*(-clhs0*clhs18 - clhs10*clhs20 + clhs14*clhs2 + clhs16*clhs5) + NormalSlave(0,1)*(-clhs0*clhs7 + clhs1*clhs2 - clhs10*clhs9 + clhs4*clhs5));
        const double clhs30 =     LM(0,0)*ScaleFactor - NormalSlave(0,0)*clhs29 - TangentSlaveXi(0,0)*clhs28;
        const double clhs31 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs32 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs33 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs34 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs35 =     NormalSlave(0,0)*(clhs14*clhs31 + clhs16*clhs32 - clhs18*clhs33 - clhs20*clhs34);
        const double clhs36 =     clhs1*clhs31;
        const double clhs37 =     clhs32*clhs4;
        const double clhs38 =     clhs33*clhs7;
        const double clhs39 =     clhs34*clhs9;
        const double clhs40 =     NormalSlave(0,1)*(clhs13 + clhs36 + clhs37 - clhs38 - clhs39) + clhs35;
        const double clhs41 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs40;
        const double clhs42 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs43 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs44 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs45 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs46 =     NormalSlave(0,1)*(clhs1*clhs42 + clhs4*clhs43 - clhs44*clhs7 - clhs45*clhs9);
        const double clhs47 =     -clhs10;
        const double clhs48 =     clhs14*clhs42;
        const double clhs49 =     clhs16*clhs43;
        const double clhs50 =     clhs18*clhs44;
        const double clhs51 =     clhs20*clhs45;
        const double clhs52 =     NormalSlave(0,0)*(clhs47 + clhs48 + clhs49 - clhs50 - clhs51) + clhs46;
        const double clhs53 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs52;
        const double clhs54 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs55 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs56 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs57 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs58 =     NormalSlave(0,0)*(clhs14*clhs54 + clhs16*clhs55 - clhs18*clhs56 - clhs20*clhs57);
        const double clhs59 =     clhs1*clhs54;
        const double clhs60 =     clhs4*clhs55;
        const double clhs61 =     clhs56*clhs7;
        const double clhs62 =     clhs57*clhs9;
        const double clhs63 =     NormalSlave(0,1)*(clhs47 + clhs59 + clhs60 - clhs61 - clhs62) + clhs58;
        const double clhs64 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs63;
        const double clhs65 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs66 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs67 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs68 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs69 =     NormalSlave(0,0)*(clhs14*clhs65 + clhs16*clhs66 - clhs18*clhs67 + clhs2 - clhs20*clhs68) + NormalSlave(0,1)*(clhs1*clhs65 + clhs4*clhs66 - clhs67*clhs7 - clhs68*clhs9);
        const double clhs70 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs69;
        const double clhs71 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs72 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs73 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs74 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs75 =     NormalSlave(0,0)*(clhs14*clhs71 + clhs16*clhs72 - clhs18*clhs73 - clhs20*clhs74) + NormalSlave(0,1)*(clhs1*clhs71 + clhs2 + clhs4*clhs72 - clhs7*clhs73 - clhs74*clhs9);
        const double clhs76 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs75;
        const double clhs77 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs78 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs79 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs80 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs81 =     NormalSlave(0,0)*(clhs14*clhs77 + clhs16*clhs78 - clhs18*clhs79 - clhs20*clhs80 + clhs5) + NormalSlave(0,1)*(clhs1*clhs77 + clhs4*clhs78 - clhs7*clhs79 - clhs80*clhs9);
        const double clhs82 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs81;
        const double clhs83 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs84 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs85 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs86 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs87 =     NormalSlave(0,0)*(clhs14*clhs83 + clhs16*clhs84 - clhs18*clhs85 - clhs20*clhs86) + NormalSlave(0,1)*(clhs1*clhs83 + clhs4*clhs84 + clhs5 - clhs7*clhs85 - clhs86*clhs9);
        const double clhs88 =     NormalSlave(0,0)*PenaltyParameter[0]*clhs87;
        const double clhs89 =     DynamicFactor[0]*ScaleFactor;
        const double clhs90 =     clhs0*clhs89;
        const double clhs91 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs22;
        const double clhs92 =     LM(0,1)*ScaleFactor - NormalSlave(0,1)*clhs29 - TangentSlaveXi(0,1)*clhs28;
        const double clhs93 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs40;
        const double clhs94 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs52;
        const double clhs95 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs63;
        const double clhs96 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs69;
        const double clhs97 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs75;
        const double clhs98 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs81;
        const double clhs99 =     NormalSlave(0,1)*PenaltyParameter[0]*clhs87;
        const double clhs100 =     clhs10*clhs89;
        const double clhs101 =     -clhs2*clhs89;
        const double clhs102 =     -clhs5*clhs89;
    
        rLocalLHS(0,0)+=DynamicFactor[0]*(-clhs0*clhs23 + clhs30*clhs8);
        rLocalLHS(0,1)+=DynamicFactor[0]*(-clhs0*clhs41 + clhs30*clhs33);
        rLocalLHS(0,2)+=DynamicFactor[0]*(-clhs0*clhs53 + clhs30*clhs44);
        rLocalLHS(0,3)+=DynamicFactor[0]*(-clhs0*clhs64 + clhs30*clhs56);
        rLocalLHS(0,4)+=DynamicFactor[0]*(-clhs0*clhs70 + clhs30*clhs67);
        rLocalLHS(0,5)+=DynamicFactor[0]*(-clhs0*clhs76 + clhs30*clhs73);
        rLocalLHS(0,6)+=DynamicFactor[0]*(-clhs0*clhs82 + clhs30*clhs79);
        rLocalLHS(0,7)+=DynamicFactor[0]*(-clhs0*clhs88 + clhs30*clhs85);
        rLocalLHS(0,8)+=clhs90;
        rLocalLHS(1,0)+=DynamicFactor[0]*(-clhs0*clhs91 + clhs8*clhs92);
        rLocalLHS(1,1)+=DynamicFactor[0]*(-clhs0*clhs93 + clhs33*clhs92);
        rLocalLHS(1,2)+=DynamicFactor[0]*(-clhs0*clhs94 + clhs44*clhs92);
        rLocalLHS(1,3)+=DynamicFactor[0]*(-clhs0*clhs95 + clhs56*clhs92);
        rLocalLHS(1,4)+=DynamicFactor[0]*(-clhs0*clhs96 + clhs67*clhs92);
        rLocalLHS(1,5)+=DynamicFactor[0]*(-clhs0*clhs97 + clhs73*clhs92);
        rLocalLHS(1,6)+=DynamicFactor[0]*(-clhs0*clhs98 + clhs79*clhs92);
        rLocalLHS(1,7)+=DynamicFactor[0]*(-clhs0*clhs99 + clhs85*clhs92);
        rLocalLHS(1,9)+=clhs90;
        rLocalLHS(2,0)+=DynamicFactor[0]*(-clhs10*clhs23 + clhs11*clhs30);
        rLocalLHS(2,1)+=DynamicFactor[0]*(-clhs10*clhs41 + clhs30*clhs34);
        rLocalLHS(2,2)+=DynamicFactor[0]*(-clhs10*clhs53 + clhs30*clhs45);
        rLocalLHS(2,3)+=DynamicFactor[0]*(-clhs10*clhs64 + clhs30*clhs57);
        rLocalLHS(2,4)+=DynamicFactor[0]*(-clhs10*clhs70 + clhs30*clhs68);
        rLocalLHS(2,5)+=DynamicFactor[0]*(-clhs10*clhs76 + clhs30*clhs74);
        rLocalLHS(2,6)+=DynamicFactor[0]*(-clhs10*clhs82 + clhs30*clhs80);
        rLocalLHS(2,7)+=DynamicFactor[0]*(-clhs10*clhs88 + clhs30*clhs86);
        rLocalLHS(2,8)+=clhs100;
        rLocalLHS(3,0)+=DynamicFactor[0]*(-clhs10*clhs91 + clhs11*clhs92);
        rLocalLHS(3,1)+=DynamicFactor[0]*(-clhs10*clhs93 + clhs34*clhs92);
        rLocalLHS(3,2)+=DynamicFactor[0]*(-clhs10*clhs94 + clhs45*clhs92);
        rLocalLHS(3,3)+=DynamicFactor[0]*(-clhs10*clhs95 + clhs57*clhs92);
        rLocalLHS(3,4)+=DynamicFactor[0]*(-clhs10*clhs96 + clhs68*clhs92);
        rLocalLHS(3,5)+=DynamicFactor[0]*(-clhs10*clhs97 + clhs74*clhs92);
        rLocalLHS(3,6)+=DynamicFactor[0]*(-clhs10*clhs98 + clhs80*clhs92);
        rLocalLHS(3,7)+=DynamicFactor[0]*(-clhs10*clhs99 + clhs86*clhs92);
        rLocalLHS(3,9)+=clhs100;
        rLocalLHS(4,0)+=DynamicFactor[0]*(clhs2*clhs23 - clhs3*clhs30);
        rLocalLHS(4,1)+=DynamicFactor[0]*(clhs2*clhs41 - clhs30*clhs31);
        rLocalLHS(4,2)+=DynamicFactor[0]*(clhs2*clhs53 - clhs30*clhs42);
        rLocalLHS(4,3)+=DynamicFactor[0]*(clhs2*clhs64 - clhs30*clhs54);
        rLocalLHS(4,4)+=DynamicFactor[0]*(clhs2*clhs70 - clhs30*clhs65);
        rLocalLHS(4,5)+=DynamicFactor[0]*(clhs2*clhs76 - clhs30*clhs71);
        rLocalLHS(4,6)+=DynamicFactor[0]*(clhs2*clhs82 - clhs30*clhs77);
        rLocalLHS(4,7)+=DynamicFactor[0]*(clhs2*clhs88 - clhs30*clhs83);
        rLocalLHS(4,8)+=clhs101;
        rLocalLHS(5,0)+=DynamicFactor[0]*(clhs2*clhs91 - clhs3*clhs92);
        rLocalLHS(5,1)+=DynamicFactor[0]*(clhs2*clhs93 - clhs31*clhs92);
        rLocalLHS(5,2)+=DynamicFactor[0]*(clhs2*clhs94 - clhs42*clhs92);
        rLocalLHS(5,3)+=DynamicFactor[0]*(clhs2*clhs95 - clhs54*clhs92);
        rLocalLHS(5,4)+=DynamicFactor[0]*(clhs2*clhs96 - clhs65*clhs92);
        rLocalLHS(5,5)+=DynamicFactor[0]*(clhs2*clhs97 - clhs71*clhs92);
        rLocalLHS(5,6)+=DynamicFactor[0]*(clhs2*clhs98 - clhs77*clhs92);
        rLocalLHS(5,7)+=DynamicFactor[0]*(clhs2*clhs99 - clhs83*clhs92);
        rLocalLHS(5,9)+=clhs101;
        rLocalLHS(6,0)+=DynamicFactor[0]*(clhs23*clhs5 - clhs30*clhs6);
        rLocalLHS(6,1)+=DynamicFactor[0]*(-clhs30*clhs32 + clhs41*clhs5);
        rLocalLHS(6,2)+=DynamicFactor[0]*(-clhs30*clhs43 + clhs5*clhs53);
        rLocalLHS(6,3)+=DynamicFactor[0]*(-clhs30*clhs55 + clhs5*clhs64);
        rLocalLHS(6,4)+=DynamicFactor[0]*(-clhs30*clhs66 + clhs5*clhs70);
        rLocalLHS(6,5)+=DynamicFactor[0]*(-clhs30*clhs72 + clhs5*clhs76);
        rLocalLHS(6,6)+=DynamicFactor[0]*(-clhs30*clhs78 + clhs5*clhs82);
        rLocalLHS(6,7)+=DynamicFactor[0]*(-clhs30*clhs84 + clhs5*clhs88);
        rLocalLHS(6,8)+=clhs102;
        rLocalLHS(7,0)+=DynamicFactor[0]*(clhs5*clhs91 - clhs6*clhs92);
        rLocalLHS(7,1)+=DynamicFactor[0]*(-clhs32*clhs92 + clhs5*clhs93);
        rLocalLHS(7,2)+=DynamicFactor[0]*(-clhs43*clhs92 + clhs5*clhs94);
        rLocalLHS(7,3)+=DynamicFactor[0]*(clhs5*clhs95 - clhs55*clhs92);
        rLocalLHS(7,4)+=DynamicFactor[0]*(clhs5*clhs96 - clhs66*clhs92);
        rLocalLHS(7,5)+=DynamicFactor[0]*(clhs5*clhs97 - clhs72*clhs92);
        rLocalLHS(7,6)+=DynamicFactor[0]*(clhs5*clhs98 - clhs78*clhs92);
        rLocalLHS(7,7)+=DynamicFactor[0]*(clhs5*clhs99 - clhs84*clhs92);
        rLocalLHS(7,9)+=clhs102;
        rLocalLHS(8,0)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs0 - clhs15 - clhs17 + clhs19 + clhs21) + clhs12);
        rLocalLHS(8,1)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs0 - clhs36 - clhs37 + clhs38 + clhs39) + clhs35);
        rLocalLHS(8,2)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs10 - clhs48 - clhs49 + clhs50 + clhs51) + clhs46);
        rLocalLHS(8,3)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs10 - clhs59 - clhs60 + clhs61 + clhs62) + clhs58);
        rLocalLHS(8,4)+=-ScaleFactor*clhs69;
        rLocalLHS(8,5)+=-ScaleFactor*clhs75;
        rLocalLHS(8,6)+=-ScaleFactor*clhs81;
        rLocalLHS(8,7)+=-ScaleFactor*clhs87;
        rLocalLHS(9,8)+=TangentSlaveXi(0,0);
        rLocalLHS(9,9)+=TangentSlaveXi(0,1);
    }
    
    // NODE 1
    if (geometry[1].IsNot(ACTIVE)) { // INACTIVE
        const double clhs0 =     std::pow(ScaleFactor, 2)/PenaltyParameter[1];
    
        rLocalLHS(10,10)+=NormalSlave(1,0)*clhs0;
        rLocalLHS(10,11)+=NormalSlave(1,1)*clhs0;
        rLocalLHS(11,10)+=TangentSlaveXi(1,0)*clhs0;
        rLocalLHS(11,11)+=TangentSlaveXi(1,1)*clhs0;
    } else if (geometry[1].Is(SLIP)) { // ACTIVE-SLIP
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     X1(0,1) + u1(0,1);
        const double clhs2 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs3 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs4 =     X1(1,1) + u1(1,1);
        const double clhs5 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X2(0,1) + u2(0,1);
        const double clhs8 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs9 =     X2(1,1) + u2(1,1);
        const double clhs10 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs11 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     NormalSlave(1,1)*(clhs1*clhs3 - clhs11*clhs9 + clhs4*clhs6 - clhs7*clhs8);
        const double clhs13 =     -clhs0;
        const double clhs14 =     X1(0,0) + u1(0,0);
        const double clhs15 =     clhs14*clhs3;
        const double clhs16 =     X1(1,0) + u1(1,0);
        const double clhs17 =     clhs16*clhs6;
        const double clhs18 =     X2(0,0) + u2(0,0);
        const double clhs19 =     clhs18*clhs8;
        const double clhs20 =     X2(1,0) + u2(1,0);
        const double clhs21 =     clhs11*clhs20;
        const double clhs22 =     NormalSlave(1,0)*(clhs13 + clhs15 + clhs17 - clhs19 - clhs21) + clhs12;
        const double clhs23 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs22;
        const double clhs24 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double clhs25 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double clhs26 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double clhs27 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double clhs28 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(clhs24*(X1(0,0) + u1old(0,0)) + clhs25*(X1(1,0) + u1old(1,0)) - clhs26*(X2(0,0) + u2old(0,0)) - clhs27*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(1,1)*(clhs24*(X1(0,1) + u1old(0,1)) + clhs25*(X1(1,1) + u1old(1,1)) - clhs26*(X2(0,1) + u2old(0,1)) - clhs27*(X2(1,1) + u2old(1,1))));
        const double clhs29 =     PenaltyParameter[1]*(NormalSlave(1,0)*(-clhs0*clhs18 - clhs10*clhs20 + clhs14*clhs2 + clhs16*clhs5) + NormalSlave(1,1)*(-clhs0*clhs7 + clhs1*clhs2 - clhs10*clhs9 + clhs4*clhs5));
        const double clhs30 =     LM(1,0)*ScaleFactor - NormalSlave(1,0)*clhs29 - TangentSlaveXi(1,0)*clhs28;
        const double clhs31 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs32 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs33 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs34 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs35 =     NormalSlave(1,0)*(clhs14*clhs31 + clhs16*clhs32 - clhs18*clhs33 - clhs20*clhs34);
        const double clhs36 =     clhs1*clhs31;
        const double clhs37 =     clhs32*clhs4;
        const double clhs38 =     clhs33*clhs7;
        const double clhs39 =     clhs34*clhs9;
        const double clhs40 =     NormalSlave(1,1)*(clhs13 + clhs36 + clhs37 - clhs38 - clhs39) + clhs35;
        const double clhs41 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs40;
        const double clhs42 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs43 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs44 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs45 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs46 =     NormalSlave(1,1)*(clhs1*clhs42 + clhs4*clhs43 - clhs44*clhs7 - clhs45*clhs9);
        const double clhs47 =     -clhs10;
        const double clhs48 =     clhs14*clhs42;
        const double clhs49 =     clhs16*clhs43;
        const double clhs50 =     clhs18*clhs44;
        const double clhs51 =     clhs20*clhs45;
        const double clhs52 =     NormalSlave(1,0)*(clhs47 + clhs48 + clhs49 - clhs50 - clhs51) + clhs46;
        const double clhs53 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs52;
        const double clhs54 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs55 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs56 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs57 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs58 =     NormalSlave(1,0)*(clhs14*clhs54 + clhs16*clhs55 - clhs18*clhs56 - clhs20*clhs57);
        const double clhs59 =     clhs1*clhs54;
        const double clhs60 =     clhs4*clhs55;
        const double clhs61 =     clhs56*clhs7;
        const double clhs62 =     clhs57*clhs9;
        const double clhs63 =     NormalSlave(1,1)*(clhs47 + clhs59 + clhs60 - clhs61 - clhs62) + clhs58;
        const double clhs64 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs63;
        const double clhs65 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs66 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs67 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs68 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs69 =     NormalSlave(1,0)*(clhs14*clhs65 + clhs16*clhs66 - clhs18*clhs67 + clhs2 - clhs20*clhs68) + NormalSlave(1,1)*(clhs1*clhs65 + clhs4*clhs66 - clhs67*clhs7 - clhs68*clhs9);
        const double clhs70 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs69;
        const double clhs71 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs72 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs73 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs74 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs75 =     NormalSlave(1,0)*(clhs14*clhs71 + clhs16*clhs72 - clhs18*clhs73 - clhs20*clhs74) + NormalSlave(1,1)*(clhs1*clhs71 + clhs2 + clhs4*clhs72 - clhs7*clhs73 - clhs74*clhs9);
        const double clhs76 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs75;
        const double clhs77 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs78 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs79 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs80 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs81 =     NormalSlave(1,0)*(clhs14*clhs77 + clhs16*clhs78 - clhs18*clhs79 - clhs20*clhs80 + clhs5) + NormalSlave(1,1)*(clhs1*clhs77 + clhs4*clhs78 - clhs7*clhs79 - clhs80*clhs9);
        const double clhs82 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs81;
        const double clhs83 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs84 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs85 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs86 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs87 =     NormalSlave(1,0)*(clhs14*clhs83 + clhs16*clhs84 - clhs18*clhs85 - clhs20*clhs86) + NormalSlave(1,1)*(clhs1*clhs83 + clhs4*clhs84 + clhs5 - clhs7*clhs85 - clhs86*clhs9);
        const double clhs88 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs87;
        const double clhs89 =     DynamicFactor[1]*ScaleFactor;
        const double clhs90 =     clhs0*clhs89;
        const double clhs91 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs22;
        const double clhs92 =     LM(1,1)*ScaleFactor - NormalSlave(1,1)*clhs29 - TangentSlaveXi(1,1)*clhs28;
        const double clhs93 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs40;
        const double clhs94 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs52;
        const double clhs95 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs63;
        const double clhs96 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs69;
        const double clhs97 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs75;
        const double clhs98 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs81;
        const double clhs99 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs87;
        const double clhs100 =     clhs10*clhs89;
        const double clhs101 =     -clhs2*clhs89;
        const double clhs102 =     -clhs5*clhs89;
    
        rLocalLHS(0,0)+=DynamicFactor[1]*(-clhs0*clhs23 + clhs30*clhs8);
        rLocalLHS(0,1)+=DynamicFactor[1]*(-clhs0*clhs41 + clhs30*clhs33);
        rLocalLHS(0,2)+=DynamicFactor[1]*(-clhs0*clhs53 + clhs30*clhs44);
        rLocalLHS(0,3)+=DynamicFactor[1]*(-clhs0*clhs64 + clhs30*clhs56);
        rLocalLHS(0,4)+=DynamicFactor[1]*(-clhs0*clhs70 + clhs30*clhs67);
        rLocalLHS(0,5)+=DynamicFactor[1]*(-clhs0*clhs76 + clhs30*clhs73);
        rLocalLHS(0,6)+=DynamicFactor[1]*(-clhs0*clhs82 + clhs30*clhs79);
        rLocalLHS(0,7)+=DynamicFactor[1]*(-clhs0*clhs88 + clhs30*clhs85);
        rLocalLHS(0,10)+=clhs90;
        rLocalLHS(1,0)+=DynamicFactor[1]*(-clhs0*clhs91 + clhs8*clhs92);
        rLocalLHS(1,1)+=DynamicFactor[1]*(-clhs0*clhs93 + clhs33*clhs92);
        rLocalLHS(1,2)+=DynamicFactor[1]*(-clhs0*clhs94 + clhs44*clhs92);
        rLocalLHS(1,3)+=DynamicFactor[1]*(-clhs0*clhs95 + clhs56*clhs92);
        rLocalLHS(1,4)+=DynamicFactor[1]*(-clhs0*clhs96 + clhs67*clhs92);
        rLocalLHS(1,5)+=DynamicFactor[1]*(-clhs0*clhs97 + clhs73*clhs92);
        rLocalLHS(1,6)+=DynamicFactor[1]*(-clhs0*clhs98 + clhs79*clhs92);
        rLocalLHS(1,7)+=DynamicFactor[1]*(-clhs0*clhs99 + clhs85*clhs92);
        rLocalLHS(1,11)+=clhs90;
        rLocalLHS(2,0)+=DynamicFactor[1]*(-clhs10*clhs23 + clhs11*clhs30);
        rLocalLHS(2,1)+=DynamicFactor[1]*(-clhs10*clhs41 + clhs30*clhs34);
        rLocalLHS(2,2)+=DynamicFactor[1]*(-clhs10*clhs53 + clhs30*clhs45);
        rLocalLHS(2,3)+=DynamicFactor[1]*(-clhs10*clhs64 + clhs30*clhs57);
        rLocalLHS(2,4)+=DynamicFactor[1]*(-clhs10*clhs70 + clhs30*clhs68);
        rLocalLHS(2,5)+=DynamicFactor[1]*(-clhs10*clhs76 + clhs30*clhs74);
        rLocalLHS(2,6)+=DynamicFactor[1]*(-clhs10*clhs82 + clhs30*clhs80);
        rLocalLHS(2,7)+=DynamicFactor[1]*(-clhs10*clhs88 + clhs30*clhs86);
        rLocalLHS(2,10)+=clhs100;
        rLocalLHS(3,0)+=DynamicFactor[1]*(-clhs10*clhs91 + clhs11*clhs92);
        rLocalLHS(3,1)+=DynamicFactor[1]*(-clhs10*clhs93 + clhs34*clhs92);
        rLocalLHS(3,2)+=DynamicFactor[1]*(-clhs10*clhs94 + clhs45*clhs92);
        rLocalLHS(3,3)+=DynamicFactor[1]*(-clhs10*clhs95 + clhs57*clhs92);
        rLocalLHS(3,4)+=DynamicFactor[1]*(-clhs10*clhs96 + clhs68*clhs92);
        rLocalLHS(3,5)+=DynamicFactor[1]*(-clhs10*clhs97 + clhs74*clhs92);
        rLocalLHS(3,6)+=DynamicFactor[1]*(-clhs10*clhs98 + clhs80*clhs92);
        rLocalLHS(3,7)+=DynamicFactor[1]*(-clhs10*clhs99 + clhs86*clhs92);
        rLocalLHS(3,11)+=clhs100;
        rLocalLHS(4,0)+=DynamicFactor[1]*(clhs2*clhs23 - clhs3*clhs30);
        rLocalLHS(4,1)+=DynamicFactor[1]*(clhs2*clhs41 - clhs30*clhs31);
        rLocalLHS(4,2)+=DynamicFactor[1]*(clhs2*clhs53 - clhs30*clhs42);
        rLocalLHS(4,3)+=DynamicFactor[1]*(clhs2*clhs64 - clhs30*clhs54);
        rLocalLHS(4,4)+=DynamicFactor[1]*(clhs2*clhs70 - clhs30*clhs65);
        rLocalLHS(4,5)+=DynamicFactor[1]*(clhs2*clhs76 - clhs30*clhs71);
        rLocalLHS(4,6)+=DynamicFactor[1]*(clhs2*clhs82 - clhs30*clhs77);
        rLocalLHS(4,7)+=DynamicFactor[1]*(clhs2*clhs88 - clhs30*clhs83);
        rLocalLHS(4,10)+=clhs101;
        rLocalLHS(5,0)+=DynamicFactor[1]*(clhs2*clhs91 - clhs3*clhs92);
        rLocalLHS(5,1)+=DynamicFactor[1]*(clhs2*clhs93 - clhs31*clhs92);
        rLocalLHS(5,2)+=DynamicFactor[1]*(clhs2*clhs94 - clhs42*clhs92);
        rLocalLHS(5,3)+=DynamicFactor[1]*(clhs2*clhs95 - clhs54*clhs92);
        rLocalLHS(5,4)+=DynamicFactor[1]*(clhs2*clhs96 - clhs65*clhs92);
        rLocalLHS(5,5)+=DynamicFactor[1]*(clhs2*clhs97 - clhs71*clhs92);
        rLocalLHS(5,6)+=DynamicFactor[1]*(clhs2*clhs98 - clhs77*clhs92);
        rLocalLHS(5,7)+=DynamicFactor[1]*(clhs2*clhs99 - clhs83*clhs92);
        rLocalLHS(5,11)+=clhs101;
        rLocalLHS(6,0)+=DynamicFactor[1]*(clhs23*clhs5 - clhs30*clhs6);
        rLocalLHS(6,1)+=DynamicFactor[1]*(-clhs30*clhs32 + clhs41*clhs5);
        rLocalLHS(6,2)+=DynamicFactor[1]*(-clhs30*clhs43 + clhs5*clhs53);
        rLocalLHS(6,3)+=DynamicFactor[1]*(-clhs30*clhs55 + clhs5*clhs64);
        rLocalLHS(6,4)+=DynamicFactor[1]*(-clhs30*clhs66 + clhs5*clhs70);
        rLocalLHS(6,5)+=DynamicFactor[1]*(-clhs30*clhs72 + clhs5*clhs76);
        rLocalLHS(6,6)+=DynamicFactor[1]*(-clhs30*clhs78 + clhs5*clhs82);
        rLocalLHS(6,7)+=DynamicFactor[1]*(-clhs30*clhs84 + clhs5*clhs88);
        rLocalLHS(6,10)+=clhs102;
        rLocalLHS(7,0)+=DynamicFactor[1]*(clhs5*clhs91 - clhs6*clhs92);
        rLocalLHS(7,1)+=DynamicFactor[1]*(-clhs32*clhs92 + clhs5*clhs93);
        rLocalLHS(7,2)+=DynamicFactor[1]*(-clhs43*clhs92 + clhs5*clhs94);
        rLocalLHS(7,3)+=DynamicFactor[1]*(clhs5*clhs95 - clhs55*clhs92);
        rLocalLHS(7,4)+=DynamicFactor[1]*(clhs5*clhs96 - clhs66*clhs92);
        rLocalLHS(7,5)+=DynamicFactor[1]*(clhs5*clhs97 - clhs72*clhs92);
        rLocalLHS(7,6)+=DynamicFactor[1]*(clhs5*clhs98 - clhs78*clhs92);
        rLocalLHS(7,7)+=DynamicFactor[1]*(clhs5*clhs99 - clhs84*clhs92);
        rLocalLHS(7,11)+=clhs102;
        rLocalLHS(10,0)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs0 - clhs15 - clhs17 + clhs19 + clhs21) + clhs12);
        rLocalLHS(10,1)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs0 - clhs36 - clhs37 + clhs38 + clhs39) + clhs35);
        rLocalLHS(10,2)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs10 - clhs48 - clhs49 + clhs50 + clhs51) + clhs46);
        rLocalLHS(10,3)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs10 - clhs59 - clhs60 + clhs61 + clhs62) + clhs58);
        rLocalLHS(10,4)+=-ScaleFactor*clhs69;
        rLocalLHS(10,5)+=-ScaleFactor*clhs75;
        rLocalLHS(10,6)+=-ScaleFactor*clhs81;
        rLocalLHS(10,7)+=-ScaleFactor*clhs87;
        rLocalLHS(11,10)+=TangentSlaveXi(0,0);
        rLocalLHS(11,11)+=TangentSlaveXi(0,1);
    } else { // ACTIVE-STICK
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     X1(0,1) + u1(0,1);
        const double clhs2 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs3 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs4 =     X1(1,1) + u1(1,1);
        const double clhs5 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X2(0,1) + u2(0,1);
        const double clhs8 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs9 =     X2(1,1) + u2(1,1);
        const double clhs10 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs11 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     NormalSlave(1,1)*(clhs1*clhs3 - clhs11*clhs9 + clhs4*clhs6 - clhs7*clhs8);
        const double clhs13 =     -clhs0;
        const double clhs14 =     X1(0,0) + u1(0,0);
        const double clhs15 =     clhs14*clhs3;
        const double clhs16 =     X1(1,0) + u1(1,0);
        const double clhs17 =     clhs16*clhs6;
        const double clhs18 =     X2(0,0) + u2(0,0);
        const double clhs19 =     clhs18*clhs8;
        const double clhs20 =     X2(1,0) + u2(1,0);
        const double clhs21 =     clhs11*clhs20;
        const double clhs22 =     NormalSlave(1,0)*(clhs13 + clhs15 + clhs17 - clhs19 - clhs21) + clhs12;
        const double clhs23 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs22;
        const double clhs24 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double clhs25 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double clhs26 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double clhs27 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double clhs28 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(clhs24*(X1(0,0) + u1old(0,0)) + clhs25*(X1(1,0) + u1old(1,0)) - clhs26*(X2(0,0) + u2old(0,0)) - clhs27*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(1,1)*(clhs24*(X1(0,1) + u1old(0,1)) + clhs25*(X1(1,1) + u1old(1,1)) - clhs26*(X2(0,1) + u2old(0,1)) - clhs27*(X2(1,1) + u2old(1,1))));
        const double clhs29 =     PenaltyParameter[1]*(NormalSlave(1,0)*(-clhs0*clhs18 - clhs10*clhs20 + clhs14*clhs2 + clhs16*clhs5) + NormalSlave(1,1)*(-clhs0*clhs7 + clhs1*clhs2 - clhs10*clhs9 + clhs4*clhs5));
        const double clhs30 =     LM(1,0)*ScaleFactor - NormalSlave(1,0)*clhs29 - TangentSlaveXi(1,0)*clhs28;
        const double clhs31 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs32 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs33 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs34 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs35 =     NormalSlave(1,0)*(clhs14*clhs31 + clhs16*clhs32 - clhs18*clhs33 - clhs20*clhs34);
        const double clhs36 =     clhs1*clhs31;
        const double clhs37 =     clhs32*clhs4;
        const double clhs38 =     clhs33*clhs7;
        const double clhs39 =     clhs34*clhs9;
        const double clhs40 =     NormalSlave(1,1)*(clhs13 + clhs36 + clhs37 - clhs38 - clhs39) + clhs35;
        const double clhs41 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs40;
        const double clhs42 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs43 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs44 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs45 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs46 =     NormalSlave(1,1)*(clhs1*clhs42 + clhs4*clhs43 - clhs44*clhs7 - clhs45*clhs9);
        const double clhs47 =     -clhs10;
        const double clhs48 =     clhs14*clhs42;
        const double clhs49 =     clhs16*clhs43;
        const double clhs50 =     clhs18*clhs44;
        const double clhs51 =     clhs20*clhs45;
        const double clhs52 =     NormalSlave(1,0)*(clhs47 + clhs48 + clhs49 - clhs50 - clhs51) + clhs46;
        const double clhs53 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs52;
        const double clhs54 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs55 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs56 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs57 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs58 =     NormalSlave(1,0)*(clhs14*clhs54 + clhs16*clhs55 - clhs18*clhs56 - clhs20*clhs57);
        const double clhs59 =     clhs1*clhs54;
        const double clhs60 =     clhs4*clhs55;
        const double clhs61 =     clhs56*clhs7;
        const double clhs62 =     clhs57*clhs9;
        const double clhs63 =     NormalSlave(1,1)*(clhs47 + clhs59 + clhs60 - clhs61 - clhs62) + clhs58;
        const double clhs64 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs63;
        const double clhs65 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs66 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs67 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs68 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs69 =     NormalSlave(1,0)*(clhs14*clhs65 + clhs16*clhs66 - clhs18*clhs67 + clhs2 - clhs20*clhs68) + NormalSlave(1,1)*(clhs1*clhs65 + clhs4*clhs66 - clhs67*clhs7 - clhs68*clhs9);
        const double clhs70 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs69;
        const double clhs71 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs72 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs73 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs74 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs75 =     NormalSlave(1,0)*(clhs14*clhs71 + clhs16*clhs72 - clhs18*clhs73 - clhs20*clhs74) + NormalSlave(1,1)*(clhs1*clhs71 + clhs2 + clhs4*clhs72 - clhs7*clhs73 - clhs74*clhs9);
        const double clhs76 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs75;
        const double clhs77 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs78 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs79 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs80 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs81 =     NormalSlave(1,0)*(clhs14*clhs77 + clhs16*clhs78 - clhs18*clhs79 - clhs20*clhs80 + clhs5) + NormalSlave(1,1)*(clhs1*clhs77 + clhs4*clhs78 - clhs7*clhs79 - clhs80*clhs9);
        const double clhs82 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs81;
        const double clhs83 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs84 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs85 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs86 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs87 =     NormalSlave(1,0)*(clhs14*clhs83 + clhs16*clhs84 - clhs18*clhs85 - clhs20*clhs86) + NormalSlave(1,1)*(clhs1*clhs83 + clhs4*clhs84 + clhs5 - clhs7*clhs85 - clhs86*clhs9);
        const double clhs88 =     NormalSlave(1,0)*PenaltyParameter[1]*clhs87;
        const double clhs89 =     DynamicFactor[1]*ScaleFactor;
        const double clhs90 =     clhs0*clhs89;
        const double clhs91 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs22;
        const double clhs92 =     LM(1,1)*ScaleFactor - NormalSlave(1,1)*clhs29 - TangentSlaveXi(1,1)*clhs28;
        const double clhs93 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs40;
        const double clhs94 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs52;
        const double clhs95 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs63;
        const double clhs96 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs69;
        const double clhs97 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs75;
        const double clhs98 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs81;
        const double clhs99 =     NormalSlave(1,1)*PenaltyParameter[1]*clhs87;
        const double clhs100 =     clhs10*clhs89;
        const double clhs101 =     -clhs2*clhs89;
        const double clhs102 =     -clhs5*clhs89;
    
        rLocalLHS(0,0)+=DynamicFactor[1]*(-clhs0*clhs23 + clhs30*clhs8);
        rLocalLHS(0,1)+=DynamicFactor[1]*(-clhs0*clhs41 + clhs30*clhs33);
        rLocalLHS(0,2)+=DynamicFactor[1]*(-clhs0*clhs53 + clhs30*clhs44);
        rLocalLHS(0,3)+=DynamicFactor[1]*(-clhs0*clhs64 + clhs30*clhs56);
        rLocalLHS(0,4)+=DynamicFactor[1]*(-clhs0*clhs70 + clhs30*clhs67);
        rLocalLHS(0,5)+=DynamicFactor[1]*(-clhs0*clhs76 + clhs30*clhs73);
        rLocalLHS(0,6)+=DynamicFactor[1]*(-clhs0*clhs82 + clhs30*clhs79);
        rLocalLHS(0,7)+=DynamicFactor[1]*(-clhs0*clhs88 + clhs30*clhs85);
        rLocalLHS(0,10)+=clhs90;
        rLocalLHS(1,0)+=DynamicFactor[1]*(-clhs0*clhs91 + clhs8*clhs92);
        rLocalLHS(1,1)+=DynamicFactor[1]*(-clhs0*clhs93 + clhs33*clhs92);
        rLocalLHS(1,2)+=DynamicFactor[1]*(-clhs0*clhs94 + clhs44*clhs92);
        rLocalLHS(1,3)+=DynamicFactor[1]*(-clhs0*clhs95 + clhs56*clhs92);
        rLocalLHS(1,4)+=DynamicFactor[1]*(-clhs0*clhs96 + clhs67*clhs92);
        rLocalLHS(1,5)+=DynamicFactor[1]*(-clhs0*clhs97 + clhs73*clhs92);
        rLocalLHS(1,6)+=DynamicFactor[1]*(-clhs0*clhs98 + clhs79*clhs92);
        rLocalLHS(1,7)+=DynamicFactor[1]*(-clhs0*clhs99 + clhs85*clhs92);
        rLocalLHS(1,11)+=clhs90;
        rLocalLHS(2,0)+=DynamicFactor[1]*(-clhs10*clhs23 + clhs11*clhs30);
        rLocalLHS(2,1)+=DynamicFactor[1]*(-clhs10*clhs41 + clhs30*clhs34);
        rLocalLHS(2,2)+=DynamicFactor[1]*(-clhs10*clhs53 + clhs30*clhs45);
        rLocalLHS(2,3)+=DynamicFactor[1]*(-clhs10*clhs64 + clhs30*clhs57);
        rLocalLHS(2,4)+=DynamicFactor[1]*(-clhs10*clhs70 + clhs30*clhs68);
        rLocalLHS(2,5)+=DynamicFactor[1]*(-clhs10*clhs76 + clhs30*clhs74);
        rLocalLHS(2,6)+=DynamicFactor[1]*(-clhs10*clhs82 + clhs30*clhs80);
        rLocalLHS(2,7)+=DynamicFactor[1]*(-clhs10*clhs88 + clhs30*clhs86);
        rLocalLHS(2,10)+=clhs100;
        rLocalLHS(3,0)+=DynamicFactor[1]*(-clhs10*clhs91 + clhs11*clhs92);
        rLocalLHS(3,1)+=DynamicFactor[1]*(-clhs10*clhs93 + clhs34*clhs92);
        rLocalLHS(3,2)+=DynamicFactor[1]*(-clhs10*clhs94 + clhs45*clhs92);
        rLocalLHS(3,3)+=DynamicFactor[1]*(-clhs10*clhs95 + clhs57*clhs92);
        rLocalLHS(3,4)+=DynamicFactor[1]*(-clhs10*clhs96 + clhs68*clhs92);
        rLocalLHS(3,5)+=DynamicFactor[1]*(-clhs10*clhs97 + clhs74*clhs92);
        rLocalLHS(3,6)+=DynamicFactor[1]*(-clhs10*clhs98 + clhs80*clhs92);
        rLocalLHS(3,7)+=DynamicFactor[1]*(-clhs10*clhs99 + clhs86*clhs92);
        rLocalLHS(3,11)+=clhs100;
        rLocalLHS(4,0)+=DynamicFactor[1]*(clhs2*clhs23 - clhs3*clhs30);
        rLocalLHS(4,1)+=DynamicFactor[1]*(clhs2*clhs41 - clhs30*clhs31);
        rLocalLHS(4,2)+=DynamicFactor[1]*(clhs2*clhs53 - clhs30*clhs42);
        rLocalLHS(4,3)+=DynamicFactor[1]*(clhs2*clhs64 - clhs30*clhs54);
        rLocalLHS(4,4)+=DynamicFactor[1]*(clhs2*clhs70 - clhs30*clhs65);
        rLocalLHS(4,5)+=DynamicFactor[1]*(clhs2*clhs76 - clhs30*clhs71);
        rLocalLHS(4,6)+=DynamicFactor[1]*(clhs2*clhs82 - clhs30*clhs77);
        rLocalLHS(4,7)+=DynamicFactor[1]*(clhs2*clhs88 - clhs30*clhs83);
        rLocalLHS(4,10)+=clhs101;
        rLocalLHS(5,0)+=DynamicFactor[1]*(clhs2*clhs91 - clhs3*clhs92);
        rLocalLHS(5,1)+=DynamicFactor[1]*(clhs2*clhs93 - clhs31*clhs92);
        rLocalLHS(5,2)+=DynamicFactor[1]*(clhs2*clhs94 - clhs42*clhs92);
        rLocalLHS(5,3)+=DynamicFactor[1]*(clhs2*clhs95 - clhs54*clhs92);
        rLocalLHS(5,4)+=DynamicFactor[1]*(clhs2*clhs96 - clhs65*clhs92);
        rLocalLHS(5,5)+=DynamicFactor[1]*(clhs2*clhs97 - clhs71*clhs92);
        rLocalLHS(5,6)+=DynamicFactor[1]*(clhs2*clhs98 - clhs77*clhs92);
        rLocalLHS(5,7)+=DynamicFactor[1]*(clhs2*clhs99 - clhs83*clhs92);
        rLocalLHS(5,11)+=clhs101;
        rLocalLHS(6,0)+=DynamicFactor[1]*(clhs23*clhs5 - clhs30*clhs6);
        rLocalLHS(6,1)+=DynamicFactor[1]*(-clhs30*clhs32 + clhs41*clhs5);
        rLocalLHS(6,2)+=DynamicFactor[1]*(-clhs30*clhs43 + clhs5*clhs53);
        rLocalLHS(6,3)+=DynamicFactor[1]*(-clhs30*clhs55 + clhs5*clhs64);
        rLocalLHS(6,4)+=DynamicFactor[1]*(-clhs30*clhs66 + clhs5*clhs70);
        rLocalLHS(6,5)+=DynamicFactor[1]*(-clhs30*clhs72 + clhs5*clhs76);
        rLocalLHS(6,6)+=DynamicFactor[1]*(-clhs30*clhs78 + clhs5*clhs82);
        rLocalLHS(6,7)+=DynamicFactor[1]*(-clhs30*clhs84 + clhs5*clhs88);
        rLocalLHS(6,10)+=clhs102;
        rLocalLHS(7,0)+=DynamicFactor[1]*(clhs5*clhs91 - clhs6*clhs92);
        rLocalLHS(7,1)+=DynamicFactor[1]*(-clhs32*clhs92 + clhs5*clhs93);
        rLocalLHS(7,2)+=DynamicFactor[1]*(-clhs43*clhs92 + clhs5*clhs94);
        rLocalLHS(7,3)+=DynamicFactor[1]*(clhs5*clhs95 - clhs55*clhs92);
        rLocalLHS(7,4)+=DynamicFactor[1]*(clhs5*clhs96 - clhs66*clhs92);
        rLocalLHS(7,5)+=DynamicFactor[1]*(clhs5*clhs97 - clhs72*clhs92);
        rLocalLHS(7,6)+=DynamicFactor[1]*(clhs5*clhs98 - clhs78*clhs92);
        rLocalLHS(7,7)+=DynamicFactor[1]*(clhs5*clhs99 - clhs84*clhs92);
        rLocalLHS(7,11)+=clhs102;
        rLocalLHS(10,0)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs0 - clhs15 - clhs17 + clhs19 + clhs21) + clhs12);
        rLocalLHS(10,1)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs0 - clhs36 - clhs37 + clhs38 + clhs39) + clhs35);
        rLocalLHS(10,2)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs10 - clhs48 - clhs49 + clhs50 + clhs51) + clhs46);
        rLocalLHS(10,3)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs10 - clhs59 - clhs60 + clhs61 + clhs62) + clhs58);
        rLocalLHS(10,4)+=-ScaleFactor*clhs69;
        rLocalLHS(10,5)+=-ScaleFactor*clhs75;
        rLocalLHS(10,6)+=-ScaleFactor*clhs81;
        rLocalLHS(10,7)+=-ScaleFactor*clhs87;
        rLocalLHS(11,10)+=TangentSlaveXi(0,0);
        rLocalLHS(11,11)+=TangentSlaveXi(0,1);
    }
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

/***************************** BEGIN AD REPLACEMENT ********************************/
/***********************************************************************************/


/***********************************************************************************/
/***********************************************************************************/

template<>
void AugmentedLagrangianMethodFrictionalMortarContactCondition<2,2, false>::CalculateLocalRHS(
    Vector& rLocalRHS,
    const MortarConditionMatrices& rMortarConditionMatrices,
    const DerivativeDataType& rDerivativeData,
    const IndexType rActiveInactive,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Initialize
    for (std::size_t i = 0; i < 12; ++i)
        rLocalRHS[i] = 0.0;

    // The geometry of the condition
    GeometryType& geometry = this->GetGeometry();

    // Initialize values
    const BoundedMatrix<double, 2, 2>& u1 = rDerivativeData.u1;
    const BoundedMatrix<double, 2, 2>& u1old = rDerivativeData.u1old;
    const BoundedMatrix<double, 2, 2>& u2 = rDerivativeData.u2;
    const BoundedMatrix<double, 2, 2>& u2old = rDerivativeData.u2old;
    const BoundedMatrix<double, 2, 2>& X1 = rDerivativeData.X1;
    const BoundedMatrix<double, 2, 2>& X2 = rDerivativeData.X2;
    
    const BoundedMatrix<double, 2, 2> LM = MortarUtilities::GetVariableMatrix<2,2>(geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);
    
    // The normal and tangent vectors
    const BoundedMatrix<double, 2, 2>& NormalSlave = rDerivativeData.NormalSlave;
    const BoundedMatrix<double, 2, 2> TangentSlave = MortarUtilities::ComputeTangentMatrix<2,2>(geometry);
    const BoundedMatrix<double, 2, 2> TangentSlaveXi = MortarUtilities::GetVariableMatrix<2,2>(geometry, TANGENT_XI, 0);
    const BoundedMatrix<double, 2, 2> TangentSlaveEta = MortarUtilities::GetVariableMatrix<2,2>(geometry, TANGENT_ETA, 0);

    // The ALM parameters
    const array_1d<double, 2> DynamicFactor = MortarUtilities::GetVariableVector<2>(geometry, DYNAMIC_FACTOR);
    const double ScaleFactor = rDerivativeData.ScaleFactor;
    const array_1d<double, 2>& PenaltyParameter = rDerivativeData.PenaltyParameter;
    const double TangentFactor = rDerivativeData.TangentFactor;
    
    // Mortar operators
    const BoundedMatrix<double, 2, 2>& MOperator = rMortarConditionMatrices.MOperator;
    const BoundedMatrix<double, 2, 2>& DOperator = rMortarConditionMatrices.DOperator;
    const BoundedMatrix<double, 2, 2>& StandardMOperator = mCurrentMortarOperators.MOperator;
    const BoundedMatrix<double, 2, 2>& StandardDOperator = mCurrentMortarOperators.DOperator;
    const BoundedMatrix<double, 2, 2>& StandardMOperatorold = mPreviousMortarOperators.MOperator;
    const BoundedMatrix<double, 2, 2>& StandardDOperatorold = mPreviousMortarOperators.DOperator;

    // We get the friction coefficient
    const array_1d<double, 2> mu = GetFrictionCoefficient();

//    // The delta time
//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    
    // NODE 0
    if (geometry[0].IsNot(ACTIVE)) { // INACTIVE
        const double crhs0 =     std::pow(ScaleFactor, 2)/PenaltyParameter[0];
    
        rLocalRHS[8]+=-crhs0*(LM(0,0)*NormalSlave(0,0) + LM(0,1)*NormalSlave(0,1));
        rLocalRHS[9]+=-crhs0*(LM(0,0)*TangentSlaveXi(0,0) + LM(0,1)*TangentSlaveXi(0,1));
    } else if (geometry[0].Is(SLIP)) { // ACTIVE-SLIP
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double crhs2 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double crhs3 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double crhs4 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double crhs5 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(crhs1*(X1(0,0) + u1old(0,0)) + crhs2*(X1(1,0) + u1old(1,0)) - crhs3*(X2(0,0) + u2old(0,0)) - crhs4*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(0,1)*(crhs1*(X1(0,1) + u1old(0,1)) + crhs2*(X1(1,1) + u1old(1,1)) - crhs3*(X2(0,1) + u2old(0,1)) - crhs4*(X2(1,1) + u2old(1,1))));
        const double crhs6 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs7 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs8 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs9 =     NormalSlave(0,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs6*(X1(0,0) + u1(0,0)) + crhs7*(X1(1,0) + u1(1,0)) - crhs8*(X2(1,0) + u2(1,0))) + NormalSlave(0,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs6*(X1(0,1) + u1(0,1)) + crhs7*(X1(1,1) + u1(1,1)) - crhs8*(X2(1,1) + u2(1,1)));
        const double crhs10 =     PenaltyParameter[0]*crhs9;
        const double crhs11 =     DynamicFactor[0]*(-LM(0,0)*ScaleFactor + NormalSlave(0,0)*crhs10 + TangentSlaveXi(0,0)*crhs5);
        const double crhs12 =     DynamicFactor[0]*(-LM(0,1)*ScaleFactor + NormalSlave(0,1)*crhs10 + TangentSlaveXi(0,1)*crhs5);
    
        rLocalRHS[0]+=crhs0*crhs11;
        rLocalRHS[1]+=crhs0*crhs12;
        rLocalRHS[2]+=crhs11*crhs8;
        rLocalRHS[3]+=crhs12*crhs8;
        rLocalRHS[4]+=-crhs11*crhs6;
        rLocalRHS[5]+=-crhs12*crhs6;
        rLocalRHS[6]+=-crhs11*crhs7;
        rLocalRHS[7]+=-crhs12*crhs7;
        rLocalRHS[8]+=ScaleFactor*crhs9;
        rLocalRHS[9]+=-1.0*(ScaleFactor*(LM(0,0)*TangentSlaveXi(0,0) + LM(0,1)*TangentSlaveXi(0,1)) - crhs5)/(PenaltyParameter[0]*TangentFactor);
    } else { // ACTIVE-STICK
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double crhs2 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double crhs3 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double crhs4 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double crhs5 =     TangentSlaveXi(0,0)*(crhs1*(X1(0,0) + u1old(0,0)) + crhs2*(X1(1,0) + u1old(1,0)) - crhs3*(X2(0,0) + u2old(0,0)) - crhs4*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(0,1)*(crhs1*(X1(0,1) + u1old(0,1)) + crhs2*(X1(1,1) + u1old(1,1)) - crhs3*(X2(0,1) + u2old(0,1)) - crhs4*(X2(1,1) + u2old(1,1)));
        const double crhs6 =     PenaltyParameter[0]*TangentFactor*crhs5;
        const double crhs7 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs8 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs9 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs10 =     NormalSlave(0,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs7*(X1(0,0) + u1(0,0)) + crhs8*(X1(1,0) + u1(1,0)) - crhs9*(X2(1,0) + u2(1,0))) + NormalSlave(0,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs7*(X1(0,1) + u1(0,1)) + crhs8*(X1(1,1) + u1(1,1)) - crhs9*(X2(1,1) + u2(1,1)));
        const double crhs11 =     PenaltyParameter[0]*crhs10;
        const double crhs12 =     DynamicFactor[0]*(-LM(0,0)*ScaleFactor + NormalSlave(0,0)*crhs11 + TangentSlaveXi(0,0)*crhs6);
        const double crhs13 =     DynamicFactor[0]*(-LM(0,1)*ScaleFactor + NormalSlave(0,1)*crhs11 + TangentSlaveXi(0,1)*crhs6);
    
        rLocalRHS[0]+=crhs0*crhs12;
        rLocalRHS[1]+=crhs0*crhs13;
        rLocalRHS[2]+=crhs12*crhs9;
        rLocalRHS[3]+=crhs13*crhs9;
        rLocalRHS[4]+=-crhs12*crhs7;
        rLocalRHS[5]+=-crhs13*crhs7;
        rLocalRHS[6]+=-crhs12*crhs8;
        rLocalRHS[7]+=-crhs13*crhs8;
        rLocalRHS[8]+=ScaleFactor*crhs10;
        rLocalRHS[9]+=-ScaleFactor*crhs5;
    }
    
    // NODE 1
    if (geometry[1].IsNot(ACTIVE)) { // INACTIVE
        const double crhs0 =     std::pow(ScaleFactor, 2)/PenaltyParameter[1];
    
        rLocalRHS[10]+=-crhs0*(LM(1,0)*NormalSlave(1,0) + LM(1,1)*NormalSlave(1,1));
        rLocalRHS[11]+=-crhs0*(LM(1,0)*TangentSlaveXi(1,0) + LM(1,1)*TangentSlaveXi(1,1));
    } else if (geometry[1].Is(SLIP)) { // ACTIVE-SLIP
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double crhs2 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double crhs3 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double crhs4 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double crhs5 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(crhs1*(X1(0,0) + u1old(0,0)) + crhs2*(X1(1,0) + u1old(1,0)) - crhs3*(X2(0,0) + u2old(0,0)) - crhs4*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(1,1)*(crhs1*(X1(0,1) + u1old(0,1)) + crhs2*(X1(1,1) + u1old(1,1)) - crhs3*(X2(0,1) + u2old(0,1)) - crhs4*(X2(1,1) + u2old(1,1))));
        const double crhs6 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs7 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs8 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs9 =     NormalSlave(1,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs6*(X1(0,0) + u1(0,0)) + crhs7*(X1(1,0) + u1(1,0)) - crhs8*(X2(1,0) + u2(1,0))) + NormalSlave(1,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs6*(X1(0,1) + u1(0,1)) + crhs7*(X1(1,1) + u1(1,1)) - crhs8*(X2(1,1) + u2(1,1)));
        const double crhs10 =     PenaltyParameter[1]*crhs9;
        const double crhs11 =     DynamicFactor[1]*(-LM(1,0)*ScaleFactor + NormalSlave(1,0)*crhs10 + TangentSlaveXi(1,0)*crhs5);
        const double crhs12 =     DynamicFactor[1]*(-LM(1,1)*ScaleFactor + NormalSlave(1,1)*crhs10 + TangentSlaveXi(1,1)*crhs5);
    
        rLocalRHS[0]+=crhs0*crhs11;
        rLocalRHS[1]+=crhs0*crhs12;
        rLocalRHS[2]+=crhs11*crhs8;
        rLocalRHS[3]+=crhs12*crhs8;
        rLocalRHS[4]+=-crhs11*crhs6;
        rLocalRHS[5]+=-crhs12*crhs6;
        rLocalRHS[6]+=-crhs11*crhs7;
        rLocalRHS[7]+=-crhs12*crhs7;
        rLocalRHS[10]+=ScaleFactor*crhs9;
        rLocalRHS[11]+=-1.0*(ScaleFactor*(LM(1,0)*TangentSlaveXi(1,0) + LM(1,1)*TangentSlaveXi(1,1)) - crhs5)/(PenaltyParameter[1]*TangentFactor);
    } else { // ACTIVE-STICK
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double crhs2 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double crhs3 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double crhs4 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double crhs5 =     TangentSlaveXi(1,0)*(crhs1*(X1(0,0) + u1old(0,0)) + crhs2*(X1(1,0) + u1old(1,0)) - crhs3*(X2(0,0) + u2old(0,0)) - crhs4*(X2(1,0) + u2old(1,0))) + TangentSlaveXi(1,1)*(crhs1*(X1(0,1) + u1old(0,1)) + crhs2*(X1(1,1) + u1old(1,1)) - crhs3*(X2(0,1) + u2old(0,1)) - crhs4*(X2(1,1) + u2old(1,1)));
        const double crhs6 =     PenaltyParameter[1]*TangentFactor*crhs5;
        const double crhs7 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs8 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs9 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs10 =     NormalSlave(1,0)*(-crhs0*(X2(0,0) + u2(0,0)) + crhs7*(X1(0,0) + u1(0,0)) + crhs8*(X1(1,0) + u1(1,0)) - crhs9*(X2(1,0) + u2(1,0))) + NormalSlave(1,1)*(-crhs0*(X2(0,1) + u2(0,1)) + crhs7*(X1(0,1) + u1(0,1)) + crhs8*(X1(1,1) + u1(1,1)) - crhs9*(X2(1,1) + u2(1,1)));
        const double crhs11 =     PenaltyParameter[1]*crhs10;
        const double crhs12 =     DynamicFactor[1]*(-LM(1,0)*ScaleFactor + NormalSlave(1,0)*crhs11 + TangentSlaveXi(1,0)*crhs6);
        const double crhs13 =     DynamicFactor[1]*(-LM(1,1)*ScaleFactor + NormalSlave(1,1)*crhs11 + TangentSlaveXi(1,1)*crhs6);
    
        rLocalRHS[0]+=crhs0*crhs12;
        rLocalRHS[1]+=crhs0*crhs13;
        rLocalRHS[2]+=crhs12*crhs9;
        rLocalRHS[3]+=crhs13*crhs9;
        rLocalRHS[4]+=-crhs12*crhs7;
        rLocalRHS[5]+=-crhs13*crhs7;
        rLocalRHS[6]+=-crhs12*crhs8;
        rLocalRHS[7]+=-crhs13*crhs8;
        rLocalRHS[10]+=ScaleFactor*crhs10;
        rLocalRHS[11]+=-ScaleFactor*crhs5;
    }
}


/****************************** END AD REPLACEMENT *********************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& CurrentProcessInfo 
    )
{
    KRATOS_TRY;   
    
    if (rResult.size() != MatrixSize)
        rResult.resize( MatrixSize, false );
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;
    
    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ) { // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rResult[index++] = master_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = master_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = master_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( DISPLACEMENT_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( DISPLACEMENT_Z ).EquationId( );
    }

    // Slave Nodes  Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_X ).EquationId( );
        rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Y ).EquationId( );
        if (TDim == 3) rResult[index++] = slave_node.GetDof( VECTOR_LAGRANGE_MULTIPLIER_Z ).EquationId( );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
void AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo 
)
{
    KRATOS_TRY;
    
    if (rConditionalDofList.size() != MatrixSize)
        rConditionalDofList.resize( MatrixSize );
    
    IndexType index = 0;
    
    /* ORDER - [ MASTER, SLAVE, LAMBDA ] */
    GeometryType& current_master = this->GetPairedGeometry();;

    // Master Nodes Displacement Equation IDs
    for ( IndexType i_master = 0; i_master < TNumNodes; ++i_master ){ // NOTE: Assuming same number of nodes for master and slave
        NodeType& master_node = current_master[i_master];
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = master_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Displacement Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_X );
        rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( DISPLACEMENT_Z );
    }

    // Slave Nodes Lambda Equation IDs
    for ( IndexType i_slave = 0; i_slave < TNumNodes; ++i_slave ) {
        NodeType& slave_node = this->GetGeometry()[ i_slave ];
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_X );
        rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Y );
        if (TDim == 3) rConditionalDofList[index++] = slave_node.pGetDof( VECTOR_LAGRANGE_MULTIPLIER_Z );
    }
    
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
int AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes,TNormalVariation>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = BaseType::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NORMAL)
    KRATOS_CHECK_VARIABLE_KEY(VECTOR_LAGRANGE_MULTIPLIER)
    KRATOS_CHECK_VARIABLE_KEY(WEIGHTED_SLIP)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VECTOR_LAGRANGE_MULTIPLIER,rnode)

        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Z, rnode)
    }

    return ierr;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<2, 2, true>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 3, true>;
template class AugmentedLagrangianMethodFrictionalMortarContactCondition<3, 4, true>;

} // Namespace Kratos
