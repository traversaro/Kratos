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
        const double clhs1 =     TangentFactor*std::pow(TangentSlaveXi(0,0), 2);
        const double clhs2 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double clhs3 =     clhs1*clhs2;
        const double clhs4 =     X1(0,1) + u1(0,1);
        const double clhs5 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X1(1,1) + u1(1,1);
        const double clhs8 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs10 =     X2(0,1) + u2(0,1);
        const double clhs11 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs14 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     NormalSlave(0,1)*(-clhs10*clhs11 - clhs12*clhs14 + clhs4*clhs6 + clhs7*clhs9);
        const double clhs16 =     -clhs0;
        const double clhs17 =     X1(0,0) + u1(0,0);
        const double clhs18 =     clhs17*clhs6;
        const double clhs19 =     X1(1,0) + u1(1,0);
        const double clhs20 =     clhs19*clhs9;
        const double clhs21 =     X2(0,0) + u2(0,0);
        const double clhs22 =     clhs11*clhs21;
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     clhs14*clhs23;
        const double clhs25 =     NormalSlave(0,0)*(clhs16 + clhs18 + clhs20 - clhs22 - clhs24) + clhs15;
        const double clhs26 =     NormalSlave(0,0)*clhs25;
        const double clhs27 =     PenaltyParameter[0]*(clhs26 - clhs3);
        const double clhs28 =     LM(0,0)*ScaleFactor;
        const double clhs29 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double clhs30 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double clhs31 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double clhs32 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(clhs17*clhs29 + clhs19*clhs30 - clhs2*clhs21 - clhs23*clhs31) + TangentSlaveXi(0,1)*(-clhs10*clhs2 - clhs12*clhs31 + clhs29*clhs4 + clhs30*clhs7));
        const double clhs33 =     TangentSlaveXi(0,0)*clhs32;
        const double clhs34 =     PenaltyParameter[0]*(NormalSlave(0,0)*(-clhs0*clhs21 - clhs13*clhs23 + clhs17*clhs5 + clhs19*clhs8) + NormalSlave(0,1)*(-clhs0*clhs10 - clhs12*clhs13 + clhs4*clhs5 + clhs7*clhs8));
        const double clhs35 =     NormalSlave(0,0)*clhs34;
        const double clhs36 =     -clhs28 + clhs33 + clhs35;
        const double clhs37 =     TangentFactor*TangentSlaveXi(0,0)*TangentSlaveXi(0,1);
        const double clhs38 =     clhs2*clhs37;
        const double clhs39 =     -clhs38;
        const double clhs40 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs41 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs42 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs43 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs44 =     NormalSlave(0,0)*(clhs17*clhs40 + clhs19*clhs41 - clhs21*clhs42 - clhs23*clhs43);
        const double clhs45 =     clhs4*clhs40;
        const double clhs46 =     clhs41*clhs7;
        const double clhs47 =     clhs10*clhs42;
        const double clhs48 =     clhs12*clhs43;
        const double clhs49 =     NormalSlave(0,1)*(clhs16 + clhs45 + clhs46 - clhs47 - clhs48) + clhs44;
        const double clhs50 =     NormalSlave(0,0)*clhs49;
        const double clhs51 =     PenaltyParameter[0]*(clhs39 + clhs50);
        const double clhs52 =     clhs1*clhs31;
        const double clhs53 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs54 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs55 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs56 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs57 =     NormalSlave(0,1)*(-clhs10*clhs55 - clhs12*clhs56 + clhs4*clhs53 + clhs54*clhs7);
        const double clhs58 =     -clhs13;
        const double clhs59 =     clhs17*clhs53;
        const double clhs60 =     clhs19*clhs54;
        const double clhs61 =     clhs21*clhs55;
        const double clhs62 =     clhs23*clhs56;
        const double clhs63 =     NormalSlave(0,0)*(clhs58 + clhs59 + clhs60 - clhs61 - clhs62) + clhs57;
        const double clhs64 =     NormalSlave(0,0)*clhs63;
        const double clhs65 =     PenaltyParameter[0]*(-clhs52 + clhs64);
        const double clhs66 =     clhs31*clhs37;
        const double clhs67 =     -clhs66;
        const double clhs68 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs69 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs70 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs71 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs72 =     NormalSlave(0,0)*(clhs17*clhs68 + clhs19*clhs69 - clhs21*clhs70 - clhs23*clhs71);
        const double clhs73 =     clhs4*clhs68;
        const double clhs74 =     clhs69*clhs7;
        const double clhs75 =     clhs10*clhs70;
        const double clhs76 =     clhs12*clhs71;
        const double clhs77 =     NormalSlave(0,1)*(clhs58 + clhs73 + clhs74 - clhs75 - clhs76) + clhs72;
        const double clhs78 =     NormalSlave(0,0)*clhs77;
        const double clhs79 =     PenaltyParameter[0]*(clhs67 + clhs78);
        const double clhs80 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs81 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs82 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs83 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs84 =     NormalSlave(0,0)*(clhs17*clhs80 + clhs19*clhs81 - clhs21*clhs82 - clhs23*clhs83 + clhs5) + NormalSlave(0,1)*(-clhs10*clhs82 - clhs12*clhs83 + clhs4*clhs80 + clhs7*clhs81);
        const double clhs85 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs84 + clhs1*clhs29);
        const double clhs86 =     clhs28 - clhs33 - clhs35;
        const double clhs87 =     clhs29*clhs37;
        const double clhs88 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs89 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs90 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs91 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs92 =     NormalSlave(0,0)*(clhs17*clhs88 + clhs19*clhs89 - clhs21*clhs90 - clhs23*clhs91) + NormalSlave(0,1)*(-clhs10*clhs90 - clhs12*clhs91 + clhs4*clhs88 + clhs5 + clhs7*clhs89);
        const double clhs93 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs92 + clhs87);
        const double clhs94 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs95 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs96 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs97 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs98 =     NormalSlave(0,0)*(clhs17*clhs94 + clhs19*clhs95 - clhs21*clhs96 - clhs23*clhs97 + clhs8) + NormalSlave(0,1)*(-clhs10*clhs96 - clhs12*clhs97 + clhs4*clhs94 + clhs7*clhs95);
        const double clhs99 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs98 + clhs1*clhs30);
        const double clhs100 =     clhs30*clhs37;
        const double clhs101 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs102 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs103 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs104 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs105 =     NormalSlave(0,0)*(clhs101*clhs17 + clhs102*clhs19 - clhs103*clhs21 - clhs104*clhs23) + NormalSlave(0,1)*(-clhs10*clhs103 + clhs101*clhs4 + clhs102*clhs7 - clhs104*clhs12 + clhs8);
        const double clhs106 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs105 + clhs100);
        const double clhs107 =     DynamicFactor[0]*ScaleFactor;
        const double clhs108 =     clhs0*clhs107;
        const double clhs109 =     NormalSlave(0,1)*clhs25;
        const double clhs110 =     PenaltyParameter[0]*(clhs109 + clhs39);
        const double clhs111 =     LM(0,1)*ScaleFactor;
        const double clhs112 =     TangentSlaveXi(0,1)*clhs32;
        const double clhs113 =     NormalSlave(0,1)*clhs34;
        const double clhs114 =     -clhs111 + clhs112 + clhs113;
        const double clhs115 =     TangentFactor*std::pow(TangentSlaveXi(0,1), 2);
        const double clhs116 =     clhs115*clhs2;
        const double clhs117 =     NormalSlave(0,1)*clhs49;
        const double clhs118 =     PenaltyParameter[0]*(-clhs116 + clhs117);
        const double clhs119 =     NormalSlave(0,1)*clhs63;
        const double clhs120 =     PenaltyParameter[0]*(clhs119 + clhs67);
        const double clhs121 =     clhs115*clhs31;
        const double clhs122 =     NormalSlave(0,1)*clhs77;
        const double clhs123 =     PenaltyParameter[0]*(-clhs121 + clhs122);
        const double clhs124 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs84 + clhs87);
        const double clhs125 =     clhs111 - clhs112 - clhs113;
        const double clhs126 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs92 + clhs115*clhs29);
        const double clhs127 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs98 + clhs100);
        const double clhs128 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs105 + clhs115*clhs30);
        const double clhs129 =     clhs107*clhs13;
        const double clhs130 =     PenaltyParameter[0]*(-clhs26 + clhs3);
        const double clhs131 =     PenaltyParameter[0]*(clhs38 - clhs50);
        const double clhs132 =     PenaltyParameter[0]*(clhs52 - clhs64);
        const double clhs133 =     PenaltyParameter[0]*(clhs66 - clhs78);
        const double clhs134 =     -clhs107*clhs5;
        const double clhs135 =     PenaltyParameter[0]*(-clhs109 + clhs38);
        const double clhs136 =     PenaltyParameter[0]*(clhs116 - clhs117);
        const double clhs137 =     PenaltyParameter[0]*(-clhs119 + clhs66);
        const double clhs138 =     PenaltyParameter[0]*(clhs121 - clhs122);
        const double clhs139 =     -clhs107*clhs8;
        const double clhs140 =     1.0*TangentSlaveXi(0,0);
        const double clhs141 =     1.0*TangentSlaveXi(0,1);
    
        rLocalLHS(0,0)+=-DynamicFactor[0]*(clhs0*clhs27 + clhs11*clhs36);
        rLocalLHS(0,1)+=-DynamicFactor[0]*(clhs0*clhs51 + clhs36*clhs42);
        rLocalLHS(0,2)+=-DynamicFactor[0]*(clhs0*clhs65 + clhs36*clhs55);
        rLocalLHS(0,3)+=-DynamicFactor[0]*(clhs0*clhs79 + clhs36*clhs70);
        rLocalLHS(0,4)+=DynamicFactor[0]*(-clhs0*clhs85 + clhs82*clhs86);
        rLocalLHS(0,5)+=DynamicFactor[0]*(-clhs0*clhs93 + clhs86*clhs90);
        rLocalLHS(0,6)+=DynamicFactor[0]*(-clhs0*clhs99 + clhs86*clhs96);
        rLocalLHS(0,7)+=DynamicFactor[0]*(-clhs0*clhs106 + clhs103*clhs86);
        rLocalLHS(0,8)+=clhs108;
        rLocalLHS(1,0)+=-DynamicFactor[0]*(clhs0*clhs110 + clhs11*clhs114);
        rLocalLHS(1,1)+=-DynamicFactor[0]*(clhs0*clhs118 + clhs114*clhs42);
        rLocalLHS(1,2)+=-DynamicFactor[0]*(clhs0*clhs120 + clhs114*clhs55);
        rLocalLHS(1,3)+=-DynamicFactor[0]*(clhs0*clhs123 + clhs114*clhs70);
        rLocalLHS(1,4)+=DynamicFactor[0]*(-clhs0*clhs124 + clhs125*clhs82);
        rLocalLHS(1,5)+=DynamicFactor[0]*(-clhs0*clhs126 + clhs125*clhs90);
        rLocalLHS(1,6)+=DynamicFactor[0]*(-clhs0*clhs127 + clhs125*clhs96);
        rLocalLHS(1,7)+=DynamicFactor[0]*(-clhs0*clhs128 + clhs103*clhs125);
        rLocalLHS(1,9)+=clhs108;
        rLocalLHS(2,0)+=-DynamicFactor[0]*(clhs13*clhs27 + clhs14*clhs36);
        rLocalLHS(2,1)+=-DynamicFactor[0]*(clhs13*clhs51 + clhs36*clhs43);
        rLocalLHS(2,2)+=-DynamicFactor[0]*(clhs13*clhs65 + clhs36*clhs56);
        rLocalLHS(2,3)+=-DynamicFactor[0]*(clhs13*clhs79 + clhs36*clhs71);
        rLocalLHS(2,4)+=DynamicFactor[0]*(-clhs13*clhs85 + clhs83*clhs86);
        rLocalLHS(2,5)+=DynamicFactor[0]*(-clhs13*clhs93 + clhs86*clhs91);
        rLocalLHS(2,6)+=DynamicFactor[0]*(-clhs13*clhs99 + clhs86*clhs97);
        rLocalLHS(2,7)+=DynamicFactor[0]*(clhs104*clhs86 - clhs106*clhs13);
        rLocalLHS(2,8)+=clhs129;
        rLocalLHS(3,0)+=-DynamicFactor[0]*(clhs110*clhs13 + clhs114*clhs14);
        rLocalLHS(3,1)+=-DynamicFactor[0]*(clhs114*clhs43 + clhs118*clhs13);
        rLocalLHS(3,2)+=-DynamicFactor[0]*(clhs114*clhs56 + clhs120*clhs13);
        rLocalLHS(3,3)+=-DynamicFactor[0]*(clhs114*clhs71 + clhs123*clhs13);
        rLocalLHS(3,4)+=DynamicFactor[0]*(-clhs124*clhs13 + clhs125*clhs83);
        rLocalLHS(3,5)+=DynamicFactor[0]*(clhs125*clhs91 - clhs126*clhs13);
        rLocalLHS(3,6)+=DynamicFactor[0]*(clhs125*clhs97 - clhs127*clhs13);
        rLocalLHS(3,7)+=DynamicFactor[0]*(clhs104*clhs125 - clhs128*clhs13);
        rLocalLHS(3,9)+=clhs129;
        rLocalLHS(4,0)+=-DynamicFactor[0]*(clhs130*clhs5 + clhs6*clhs86);
        rLocalLHS(4,1)+=-DynamicFactor[0]*(clhs131*clhs5 + clhs40*clhs86);
        rLocalLHS(4,2)+=-DynamicFactor[0]*(clhs132*clhs5 + clhs53*clhs86);
        rLocalLHS(4,3)+=-DynamicFactor[0]*(clhs133*clhs5 + clhs68*clhs86);
        rLocalLHS(4,4)+=DynamicFactor[0]*(clhs5*clhs85 - clhs80*clhs86);
        rLocalLHS(4,5)+=DynamicFactor[0]*(clhs5*clhs93 - clhs86*clhs88);
        rLocalLHS(4,6)+=DynamicFactor[0]*(clhs5*clhs99 - clhs86*clhs94);
        rLocalLHS(4,7)+=DynamicFactor[0]*(-clhs101*clhs86 + clhs106*clhs5);
        rLocalLHS(4,8)+=clhs134;
        rLocalLHS(5,0)+=-DynamicFactor[0]*(clhs125*clhs6 + clhs135*clhs5);
        rLocalLHS(5,1)+=-DynamicFactor[0]*(clhs125*clhs40 + clhs136*clhs5);
        rLocalLHS(5,2)+=-DynamicFactor[0]*(clhs125*clhs53 + clhs137*clhs5);
        rLocalLHS(5,3)+=-DynamicFactor[0]*(clhs125*clhs68 + clhs138*clhs5);
        rLocalLHS(5,4)+=DynamicFactor[0]*(clhs124*clhs5 - clhs125*clhs80);
        rLocalLHS(5,5)+=DynamicFactor[0]*(-clhs125*clhs88 + clhs126*clhs5);
        rLocalLHS(5,6)+=DynamicFactor[0]*(-clhs125*clhs94 + clhs127*clhs5);
        rLocalLHS(5,7)+=DynamicFactor[0]*(-clhs101*clhs125 + clhs128*clhs5);
        rLocalLHS(5,9)+=clhs134;
        rLocalLHS(6,0)+=-DynamicFactor[0]*(clhs130*clhs8 + clhs86*clhs9);
        rLocalLHS(6,1)+=-DynamicFactor[0]*(clhs131*clhs8 + clhs41*clhs86);
        rLocalLHS(6,2)+=-DynamicFactor[0]*(clhs132*clhs8 + clhs54*clhs86);
        rLocalLHS(6,3)+=-DynamicFactor[0]*(clhs133*clhs8 + clhs69*clhs86);
        rLocalLHS(6,4)+=DynamicFactor[0]*(clhs8*clhs85 - clhs81*clhs86);
        rLocalLHS(6,5)+=DynamicFactor[0]*(clhs8*clhs93 - clhs86*clhs89);
        rLocalLHS(6,6)+=DynamicFactor[0]*(clhs8*clhs99 - clhs86*clhs95);
        rLocalLHS(6,7)+=DynamicFactor[0]*(-clhs102*clhs86 + clhs106*clhs8);
        rLocalLHS(6,8)+=clhs139;
        rLocalLHS(7,0)+=-DynamicFactor[0]*(clhs125*clhs9 + clhs135*clhs8);
        rLocalLHS(7,1)+=-DynamicFactor[0]*(clhs125*clhs41 + clhs136*clhs8);
        rLocalLHS(7,2)+=-DynamicFactor[0]*(clhs125*clhs54 + clhs137*clhs8);
        rLocalLHS(7,3)+=-DynamicFactor[0]*(clhs125*clhs69 + clhs138*clhs8);
        rLocalLHS(7,4)+=DynamicFactor[0]*(clhs124*clhs8 - clhs125*clhs81);
        rLocalLHS(7,5)+=DynamicFactor[0]*(-clhs125*clhs89 + clhs126*clhs8);
        rLocalLHS(7,6)+=DynamicFactor[0]*(-clhs125*clhs95 + clhs127*clhs8);
        rLocalLHS(7,7)+=DynamicFactor[0]*(-clhs102*clhs125 + clhs128*clhs8);
        rLocalLHS(7,9)+=clhs139;
        rLocalLHS(8,0)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs0 - clhs18 - clhs20 + clhs22 + clhs24) + clhs15);
        rLocalLHS(8,1)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs0 - clhs45 - clhs46 + clhs47 + clhs48) + clhs44);
        rLocalLHS(8,2)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs13 - clhs59 - clhs60 + clhs61 + clhs62) + clhs57);
        rLocalLHS(8,3)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs13 - clhs73 - clhs74 + clhs75 + clhs76) + clhs72);
        rLocalLHS(8,4)+=-ScaleFactor*clhs84;
        rLocalLHS(8,5)+=-ScaleFactor*clhs92;
        rLocalLHS(8,6)+=-ScaleFactor*clhs98;
        rLocalLHS(8,7)+=-ScaleFactor*clhs105;
        rLocalLHS(9,0)+=clhs140*clhs2;
        rLocalLHS(9,1)+=clhs141*clhs2;
        rLocalLHS(9,2)+=clhs140*clhs31;
        rLocalLHS(9,3)+=clhs141*clhs31;
        rLocalLHS(9,4)+=-clhs140*clhs29;
        rLocalLHS(9,5)+=-clhs141*clhs29;
        rLocalLHS(9,6)+=-clhs140*clhs30;
        rLocalLHS(9,7)+=-clhs141*clhs30;
        rLocalLHS(9,8)+=TangentSlaveXi(0,0);
        rLocalLHS(9,9)+=TangentSlaveXi(0,1);
    } else { // ACTIVE-STICK
        const double clhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     TangentFactor*std::pow(TangentSlaveXi(0,0), 2);
        const double clhs2 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double clhs3 =     clhs1*clhs2;
        const double clhs4 =     X1(0,1) + u1(0,1);
        const double clhs5 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X1(1,1) + u1(1,1);
        const double clhs8 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     DeltaDOperator[4](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs10 =     X2(0,1) + u2(0,1);
        const double clhs11 =     DeltaMOperator[4](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs14 =     DeltaMOperator[4](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     NormalSlave(0,1)*(-clhs10*clhs11 - clhs12*clhs14 + clhs4*clhs6 + clhs7*clhs9);
        const double clhs16 =     -clhs0;
        const double clhs17 =     X1(0,0) + u1(0,0);
        const double clhs18 =     clhs17*clhs6;
        const double clhs19 =     X1(1,0) + u1(1,0);
        const double clhs20 =     clhs19*clhs9;
        const double clhs21 =     X2(0,0) + u2(0,0);
        const double clhs22 =     clhs11*clhs21;
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     clhs14*clhs23;
        const double clhs25 =     NormalSlave(0,0)*(clhs16 + clhs18 + clhs20 - clhs22 - clhs24) + clhs15;
        const double clhs26 =     NormalSlave(0,0)*clhs25;
        const double clhs27 =     PenaltyParameter[0]*(clhs26 - clhs3);
        const double clhs28 =     LM(0,0)*ScaleFactor;
        const double clhs29 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double clhs30 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double clhs31 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double clhs32 =     TangentSlaveXi(0,0)*(clhs17*clhs29 + clhs19*clhs30 - clhs2*clhs21 - clhs23*clhs31) + TangentSlaveXi(0,1)*(-clhs10*clhs2 - clhs12*clhs31 + clhs29*clhs4 + clhs30*clhs7);
        const double clhs33 =     PenaltyParameter[0]*TangentFactor*clhs32;
        const double clhs34 =     TangentSlaveXi(0,0)*clhs33;
        const double clhs35 =     PenaltyParameter[0]*(NormalSlave(0,0)*(-clhs0*clhs21 - clhs13*clhs23 + clhs17*clhs5 + clhs19*clhs8) + NormalSlave(0,1)*(-clhs0*clhs10 - clhs12*clhs13 + clhs4*clhs5 + clhs7*clhs8));
        const double clhs36 =     NormalSlave(0,0)*clhs35;
        const double clhs37 =     -clhs28 + clhs34 + clhs36;
        const double clhs38 =     TangentFactor*TangentSlaveXi(0,0)*TangentSlaveXi(0,1);
        const double clhs39 =     clhs2*clhs38;
        const double clhs40 =     -clhs39;
        const double clhs41 =     DeltaDOperator[5](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs42 =     DeltaDOperator[5](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs43 =     DeltaMOperator[5](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs44 =     DeltaMOperator[5](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs45 =     NormalSlave(0,0)*(clhs17*clhs41 + clhs19*clhs42 - clhs21*clhs43 - clhs23*clhs44);
        const double clhs46 =     clhs4*clhs41;
        const double clhs47 =     clhs42*clhs7;
        const double clhs48 =     clhs10*clhs43;
        const double clhs49 =     clhs12*clhs44;
        const double clhs50 =     NormalSlave(0,1)*(clhs16 + clhs46 + clhs47 - clhs48 - clhs49) + clhs45;
        const double clhs51 =     NormalSlave(0,0)*clhs50;
        const double clhs52 =     PenaltyParameter[0]*(clhs40 + clhs51);
        const double clhs53 =     clhs1*clhs31;
        const double clhs54 =     DeltaDOperator[6](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs55 =     DeltaDOperator[6](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs56 =     DeltaMOperator[6](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs57 =     DeltaMOperator[6](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs58 =     NormalSlave(0,1)*(-clhs10*clhs56 - clhs12*clhs57 + clhs4*clhs54 + clhs55*clhs7);
        const double clhs59 =     -clhs13;
        const double clhs60 =     clhs17*clhs54;
        const double clhs61 =     clhs19*clhs55;
        const double clhs62 =     clhs21*clhs56;
        const double clhs63 =     clhs23*clhs57;
        const double clhs64 =     NormalSlave(0,0)*(clhs59 + clhs60 + clhs61 - clhs62 - clhs63) + clhs58;
        const double clhs65 =     NormalSlave(0,0)*clhs64;
        const double clhs66 =     PenaltyParameter[0]*(-clhs53 + clhs65);
        const double clhs67 =     clhs31*clhs38;
        const double clhs68 =     -clhs67;
        const double clhs69 =     DeltaDOperator[7](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs70 =     DeltaDOperator[7](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs71 =     DeltaMOperator[7](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs72 =     DeltaMOperator[7](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs73 =     NormalSlave(0,0)*(clhs17*clhs69 + clhs19*clhs70 - clhs21*clhs71 - clhs23*clhs72);
        const double clhs74 =     clhs4*clhs69;
        const double clhs75 =     clhs7*clhs70;
        const double clhs76 =     clhs10*clhs71;
        const double clhs77 =     clhs12*clhs72;
        const double clhs78 =     NormalSlave(0,1)*(clhs59 + clhs74 + clhs75 - clhs76 - clhs77) + clhs73;
        const double clhs79 =     NormalSlave(0,0)*clhs78;
        const double clhs80 =     PenaltyParameter[0]*(clhs68 + clhs79);
        const double clhs81 =     DeltaDOperator[0](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs82 =     DeltaDOperator[0](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs83 =     DeltaMOperator[0](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs84 =     DeltaMOperator[0](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs85 =     NormalSlave(0,0)*(clhs17*clhs81 + clhs19*clhs82 - clhs21*clhs83 - clhs23*clhs84 + clhs5) + NormalSlave(0,1)*(-clhs10*clhs83 - clhs12*clhs84 + clhs4*clhs81 + clhs7*clhs82);
        const double clhs86 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs85 + clhs1*clhs29);
        const double clhs87 =     clhs28 - clhs34 - clhs36;
        const double clhs88 =     clhs29*clhs38;
        const double clhs89 =     DeltaDOperator[1](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs90 =     DeltaDOperator[1](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs91 =     DeltaMOperator[1](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs92 =     DeltaMOperator[1](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs93 =     NormalSlave(0,0)*(clhs17*clhs89 + clhs19*clhs90 - clhs21*clhs91 - clhs23*clhs92) + NormalSlave(0,1)*(-clhs10*clhs91 - clhs12*clhs92 + clhs4*clhs89 + clhs5 + clhs7*clhs90);
        const double clhs94 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs93 + clhs88);
        const double clhs95 =     DeltaDOperator[2](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs96 =     DeltaDOperator[2](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs97 =     DeltaMOperator[2](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs98 =     DeltaMOperator[2](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs99 =     NormalSlave(0,0)*(clhs17*clhs95 + clhs19*clhs96 - clhs21*clhs97 - clhs23*clhs98 + clhs8) + NormalSlave(0,1)*(-clhs10*clhs97 - clhs12*clhs98 + clhs4*clhs95 + clhs7*clhs96);
        const double clhs100 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs99 + clhs1*clhs30);
        const double clhs101 =     clhs30*clhs38;
        const double clhs102 =     DeltaDOperator[3](0,0); // DERIVATIVE(DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs103 =     DeltaDOperator[3](0,1); // DERIVATIVE(DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs104 =     DeltaMOperator[3](0,0); // DERIVATIVE(MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs105 =     DeltaMOperator[3](0,1); // DERIVATIVE(MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs106 =     NormalSlave(0,0)*(clhs102*clhs17 + clhs103*clhs19 - clhs104*clhs21 - clhs105*clhs23) + NormalSlave(0,1)*(-clhs10*clhs104 + clhs102*clhs4 + clhs103*clhs7 - clhs105*clhs12 + clhs8);
        const double clhs107 =     PenaltyParameter[0]*(NormalSlave(0,0)*clhs106 + clhs101);
        const double clhs108 =     DynamicFactor[0]*ScaleFactor;
        const double clhs109 =     clhs0*clhs108;
        const double clhs110 =     NormalSlave(0,1)*clhs25;
        const double clhs111 =     PenaltyParameter[0]*(clhs110 + clhs40);
        const double clhs112 =     LM(0,1)*ScaleFactor;
        const double clhs113 =     TangentSlaveXi(0,1)*clhs33;
        const double clhs114 =     NormalSlave(0,1)*clhs35;
        const double clhs115 =     -clhs112 + clhs113 + clhs114;
        const double clhs116 =     TangentFactor*std::pow(TangentSlaveXi(0,1), 2);
        const double clhs117 =     clhs116*clhs2;
        const double clhs118 =     NormalSlave(0,1)*clhs50;
        const double clhs119 =     PenaltyParameter[0]*(-clhs117 + clhs118);
        const double clhs120 =     NormalSlave(0,1)*clhs64;
        const double clhs121 =     PenaltyParameter[0]*(clhs120 + clhs68);
        const double clhs122 =     clhs116*clhs31;
        const double clhs123 =     NormalSlave(0,1)*clhs78;
        const double clhs124 =     PenaltyParameter[0]*(-clhs122 + clhs123);
        const double clhs125 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs85 + clhs88);
        const double clhs126 =     clhs112 - clhs113 - clhs114;
        const double clhs127 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs93 + clhs116*clhs29);
        const double clhs128 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs99 + clhs101);
        const double clhs129 =     PenaltyParameter[0]*(NormalSlave(0,1)*clhs106 + clhs116*clhs30);
        const double clhs130 =     clhs108*clhs13;
        const double clhs131 =     PenaltyParameter[0]*(-clhs26 + clhs3);
        const double clhs132 =     PenaltyParameter[0]*(clhs39 - clhs51);
        const double clhs133 =     PenaltyParameter[0]*(clhs53 - clhs65);
        const double clhs134 =     PenaltyParameter[0]*(clhs67 - clhs79);
        const double clhs135 =     -clhs108*clhs5;
        const double clhs136 =     PenaltyParameter[0]*(-clhs110 + clhs39);
        const double clhs137 =     PenaltyParameter[0]*(clhs117 - clhs118);
        const double clhs138 =     PenaltyParameter[0]*(-clhs120 + clhs67);
        const double clhs139 =     PenaltyParameter[0]*(clhs122 - clhs123);
        const double clhs140 =     -clhs108*clhs8;
        const double clhs141 =     PenaltyParameter[0]*TangentFactor*mu[0];
        const double clhs142 =     ScaleFactor*(LM(0,0)*NormalSlave(0,0) + LM(0,1)*NormalSlave(0,1));
        const double clhs143 =     clhs142 - clhs35;
        const double clhs144 =     TangentSlaveXi(0,0)*clhs143;
        const double clhs145 =     PenaltyParameter[0]*clhs32;
        const double clhs146 =     TangentSlaveXi(0,1)*clhs143;
        const double clhs147 =     -clhs142 + clhs35;
        const double clhs148 =     TangentSlaveXi(0,0)*clhs147;
        const double clhs149 =     TangentSlaveXi(0,1)*clhs147;
    
        rLocalLHS(0,0)+=-DynamicFactor[0]*(clhs0*clhs27 + clhs11*clhs37);
        rLocalLHS(0,1)+=-DynamicFactor[0]*(clhs0*clhs52 + clhs37*clhs43);
        rLocalLHS(0,2)+=-DynamicFactor[0]*(clhs0*clhs66 + clhs37*clhs56);
        rLocalLHS(0,3)+=-DynamicFactor[0]*(clhs0*clhs80 + clhs37*clhs71);
        rLocalLHS(0,4)+=DynamicFactor[0]*(-clhs0*clhs86 + clhs83*clhs87);
        rLocalLHS(0,5)+=DynamicFactor[0]*(-clhs0*clhs94 + clhs87*clhs91);
        rLocalLHS(0,6)+=DynamicFactor[0]*(-clhs0*clhs100 + clhs87*clhs97);
        rLocalLHS(0,7)+=DynamicFactor[0]*(-clhs0*clhs107 + clhs104*clhs87);
        rLocalLHS(0,8)+=clhs109;
        rLocalLHS(1,0)+=-DynamicFactor[0]*(clhs0*clhs111 + clhs11*clhs115);
        rLocalLHS(1,1)+=-DynamicFactor[0]*(clhs0*clhs119 + clhs115*clhs43);
        rLocalLHS(1,2)+=-DynamicFactor[0]*(clhs0*clhs121 + clhs115*clhs56);
        rLocalLHS(1,3)+=-DynamicFactor[0]*(clhs0*clhs124 + clhs115*clhs71);
        rLocalLHS(1,4)+=DynamicFactor[0]*(-clhs0*clhs125 + clhs126*clhs83);
        rLocalLHS(1,5)+=DynamicFactor[0]*(-clhs0*clhs127 + clhs126*clhs91);
        rLocalLHS(1,6)+=DynamicFactor[0]*(-clhs0*clhs128 + clhs126*clhs97);
        rLocalLHS(1,7)+=DynamicFactor[0]*(-clhs0*clhs129 + clhs104*clhs126);
        rLocalLHS(1,9)+=clhs109;
        rLocalLHS(2,0)+=-DynamicFactor[0]*(clhs13*clhs27 + clhs14*clhs37);
        rLocalLHS(2,1)+=-DynamicFactor[0]*(clhs13*clhs52 + clhs37*clhs44);
        rLocalLHS(2,2)+=-DynamicFactor[0]*(clhs13*clhs66 + clhs37*clhs57);
        rLocalLHS(2,3)+=-DynamicFactor[0]*(clhs13*clhs80 + clhs37*clhs72);
        rLocalLHS(2,4)+=DynamicFactor[0]*(-clhs13*clhs86 + clhs84*clhs87);
        rLocalLHS(2,5)+=DynamicFactor[0]*(-clhs13*clhs94 + clhs87*clhs92);
        rLocalLHS(2,6)+=DynamicFactor[0]*(-clhs100*clhs13 + clhs87*clhs98);
        rLocalLHS(2,7)+=DynamicFactor[0]*(clhs105*clhs87 - clhs107*clhs13);
        rLocalLHS(2,8)+=clhs130;
        rLocalLHS(3,0)+=-DynamicFactor[0]*(clhs111*clhs13 + clhs115*clhs14);
        rLocalLHS(3,1)+=-DynamicFactor[0]*(clhs115*clhs44 + clhs119*clhs13);
        rLocalLHS(3,2)+=-DynamicFactor[0]*(clhs115*clhs57 + clhs121*clhs13);
        rLocalLHS(3,3)+=-DynamicFactor[0]*(clhs115*clhs72 + clhs124*clhs13);
        rLocalLHS(3,4)+=DynamicFactor[0]*(-clhs125*clhs13 + clhs126*clhs84);
        rLocalLHS(3,5)+=DynamicFactor[0]*(clhs126*clhs92 - clhs127*clhs13);
        rLocalLHS(3,6)+=DynamicFactor[0]*(clhs126*clhs98 - clhs128*clhs13);
        rLocalLHS(3,7)+=DynamicFactor[0]*(clhs105*clhs126 - clhs129*clhs13);
        rLocalLHS(3,9)+=clhs130;
        rLocalLHS(4,0)+=-DynamicFactor[0]*(clhs131*clhs5 + clhs6*clhs87);
        rLocalLHS(4,1)+=-DynamicFactor[0]*(clhs132*clhs5 + clhs41*clhs87);
        rLocalLHS(4,2)+=-DynamicFactor[0]*(clhs133*clhs5 + clhs54*clhs87);
        rLocalLHS(4,3)+=-DynamicFactor[0]*(clhs134*clhs5 + clhs69*clhs87);
        rLocalLHS(4,4)+=DynamicFactor[0]*(clhs5*clhs86 - clhs81*clhs87);
        rLocalLHS(4,5)+=DynamicFactor[0]*(clhs5*clhs94 - clhs87*clhs89);
        rLocalLHS(4,6)+=DynamicFactor[0]*(clhs100*clhs5 - clhs87*clhs95);
        rLocalLHS(4,7)+=DynamicFactor[0]*(-clhs102*clhs87 + clhs107*clhs5);
        rLocalLHS(4,8)+=clhs135;
        rLocalLHS(5,0)+=-DynamicFactor[0]*(clhs126*clhs6 + clhs136*clhs5);
        rLocalLHS(5,1)+=-DynamicFactor[0]*(clhs126*clhs41 + clhs137*clhs5);
        rLocalLHS(5,2)+=-DynamicFactor[0]*(clhs126*clhs54 + clhs138*clhs5);
        rLocalLHS(5,3)+=-DynamicFactor[0]*(clhs126*clhs69 + clhs139*clhs5);
        rLocalLHS(5,4)+=DynamicFactor[0]*(clhs125*clhs5 - clhs126*clhs81);
        rLocalLHS(5,5)+=DynamicFactor[0]*(-clhs126*clhs89 + clhs127*clhs5);
        rLocalLHS(5,6)+=DynamicFactor[0]*(-clhs126*clhs95 + clhs128*clhs5);
        rLocalLHS(5,7)+=DynamicFactor[0]*(-clhs102*clhs126 + clhs129*clhs5);
        rLocalLHS(5,9)+=clhs135;
        rLocalLHS(6,0)+=-DynamicFactor[0]*(clhs131*clhs8 + clhs87*clhs9);
        rLocalLHS(6,1)+=-DynamicFactor[0]*(clhs132*clhs8 + clhs42*clhs87);
        rLocalLHS(6,2)+=-DynamicFactor[0]*(clhs133*clhs8 + clhs55*clhs87);
        rLocalLHS(6,3)+=-DynamicFactor[0]*(clhs134*clhs8 + clhs70*clhs87);
        rLocalLHS(6,4)+=DynamicFactor[0]*(clhs8*clhs86 - clhs82*clhs87);
        rLocalLHS(6,5)+=DynamicFactor[0]*(clhs8*clhs94 - clhs87*clhs90);
        rLocalLHS(6,6)+=DynamicFactor[0]*(clhs100*clhs8 - clhs87*clhs96);
        rLocalLHS(6,7)+=DynamicFactor[0]*(-clhs103*clhs87 + clhs107*clhs8);
        rLocalLHS(6,8)+=clhs140;
        rLocalLHS(7,0)+=-DynamicFactor[0]*(clhs126*clhs9 + clhs136*clhs8);
        rLocalLHS(7,1)+=-DynamicFactor[0]*(clhs126*clhs42 + clhs137*clhs8);
        rLocalLHS(7,2)+=-DynamicFactor[0]*(clhs126*clhs55 + clhs138*clhs8);
        rLocalLHS(7,3)+=-DynamicFactor[0]*(clhs126*clhs70 + clhs139*clhs8);
        rLocalLHS(7,4)+=DynamicFactor[0]*(clhs125*clhs8 - clhs126*clhs82);
        rLocalLHS(7,5)+=DynamicFactor[0]*(-clhs126*clhs90 + clhs127*clhs8);
        rLocalLHS(7,6)+=DynamicFactor[0]*(-clhs126*clhs96 + clhs128*clhs8);
        rLocalLHS(7,7)+=DynamicFactor[0]*(-clhs103*clhs126 + clhs129*clhs8);
        rLocalLHS(7,9)+=clhs140;
        rLocalLHS(8,0)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs0 - clhs18 - clhs20 + clhs22 + clhs24) + clhs15);
        rLocalLHS(8,1)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs0 - clhs46 - clhs47 + clhs48 + clhs49) + clhs45);
        rLocalLHS(8,2)+=-ScaleFactor*(-NormalSlave(0,0)*(clhs13 - clhs60 - clhs61 + clhs62 + clhs63) + clhs58);
        rLocalLHS(8,3)+=-ScaleFactor*(-NormalSlave(0,1)*(clhs13 - clhs74 - clhs75 + clhs76 + clhs77) + clhs73);
        rLocalLHS(8,4)+=-ScaleFactor*clhs85;
        rLocalLHS(8,5)+=-ScaleFactor*clhs93;
        rLocalLHS(8,6)+=-ScaleFactor*clhs99;
        rLocalLHS(8,7)+=-ScaleFactor*clhs106;
        rLocalLHS(9,0)+=-clhs141*(clhs144*clhs2 + clhs145*clhs25);
        rLocalLHS(9,1)+=-clhs141*(clhs145*clhs50 + clhs146*clhs2);
        rLocalLHS(9,2)+=-clhs141*(clhs144*clhs31 + clhs145*clhs64);
        rLocalLHS(9,3)+=-clhs141*(clhs145*clhs78 + clhs146*clhs31);
        rLocalLHS(9,4)+=-clhs141*(clhs145*clhs85 + clhs148*clhs29);
        rLocalLHS(9,5)+=-clhs141*(clhs145*clhs93 + clhs149*clhs29);
        rLocalLHS(9,6)+=-clhs141*(clhs145*clhs99 + clhs148*clhs30);
        rLocalLHS(9,7)+=-clhs141*(clhs106*clhs145 + clhs149*clhs30);
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
        const double clhs1 =     TangentFactor*std::pow(TangentSlaveXi(1,0), 2);
        const double clhs2 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double clhs3 =     clhs1*clhs2;
        const double clhs4 =     X1(0,1) + u1(0,1);
        const double clhs5 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X1(1,1) + u1(1,1);
        const double clhs8 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs10 =     X2(0,1) + u2(0,1);
        const double clhs11 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs14 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     NormalSlave(1,1)*(-clhs10*clhs11 - clhs12*clhs14 + clhs4*clhs6 + clhs7*clhs9);
        const double clhs16 =     -clhs0;
        const double clhs17 =     X1(0,0) + u1(0,0);
        const double clhs18 =     clhs17*clhs6;
        const double clhs19 =     X1(1,0) + u1(1,0);
        const double clhs20 =     clhs19*clhs9;
        const double clhs21 =     X2(0,0) + u2(0,0);
        const double clhs22 =     clhs11*clhs21;
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     clhs14*clhs23;
        const double clhs25 =     NormalSlave(1,0)*(clhs16 + clhs18 + clhs20 - clhs22 - clhs24) + clhs15;
        const double clhs26 =     NormalSlave(1,0)*clhs25;
        const double clhs27 =     PenaltyParameter[1]*(clhs26 - clhs3);
        const double clhs28 =     LM(1,0)*ScaleFactor;
        const double clhs29 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double clhs30 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double clhs31 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double clhs32 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(clhs17*clhs29 + clhs19*clhs30 - clhs2*clhs21 - clhs23*clhs31) + TangentSlaveXi(1,1)*(-clhs10*clhs2 - clhs12*clhs31 + clhs29*clhs4 + clhs30*clhs7));
        const double clhs33 =     TangentSlaveXi(1,0)*clhs32;
        const double clhs34 =     PenaltyParameter[1]*(NormalSlave(1,0)*(-clhs0*clhs21 - clhs13*clhs23 + clhs17*clhs5 + clhs19*clhs8) + NormalSlave(1,1)*(-clhs0*clhs10 - clhs12*clhs13 + clhs4*clhs5 + clhs7*clhs8));
        const double clhs35 =     NormalSlave(1,0)*clhs34;
        const double clhs36 =     -clhs28 + clhs33 + clhs35;
        const double clhs37 =     TangentFactor*TangentSlaveXi(1,0)*TangentSlaveXi(1,1);
        const double clhs38 =     clhs2*clhs37;
        const double clhs39 =     -clhs38;
        const double clhs40 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs41 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs42 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs43 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs44 =     NormalSlave(1,0)*(clhs17*clhs40 + clhs19*clhs41 - clhs21*clhs42 - clhs23*clhs43);
        const double clhs45 =     clhs4*clhs40;
        const double clhs46 =     clhs41*clhs7;
        const double clhs47 =     clhs10*clhs42;
        const double clhs48 =     clhs12*clhs43;
        const double clhs49 =     NormalSlave(1,1)*(clhs16 + clhs45 + clhs46 - clhs47 - clhs48) + clhs44;
        const double clhs50 =     NormalSlave(1,0)*clhs49;
        const double clhs51 =     PenaltyParameter[1]*(clhs39 + clhs50);
        const double clhs52 =     clhs1*clhs31;
        const double clhs53 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs54 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs55 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs56 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs57 =     NormalSlave(1,1)*(-clhs10*clhs55 - clhs12*clhs56 + clhs4*clhs53 + clhs54*clhs7);
        const double clhs58 =     -clhs13;
        const double clhs59 =     clhs17*clhs53;
        const double clhs60 =     clhs19*clhs54;
        const double clhs61 =     clhs21*clhs55;
        const double clhs62 =     clhs23*clhs56;
        const double clhs63 =     NormalSlave(1,0)*(clhs58 + clhs59 + clhs60 - clhs61 - clhs62) + clhs57;
        const double clhs64 =     NormalSlave(1,0)*clhs63;
        const double clhs65 =     PenaltyParameter[1]*(-clhs52 + clhs64);
        const double clhs66 =     clhs31*clhs37;
        const double clhs67 =     -clhs66;
        const double clhs68 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs69 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs70 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs71 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs72 =     NormalSlave(1,0)*(clhs17*clhs68 + clhs19*clhs69 - clhs21*clhs70 - clhs23*clhs71);
        const double clhs73 =     clhs4*clhs68;
        const double clhs74 =     clhs69*clhs7;
        const double clhs75 =     clhs10*clhs70;
        const double clhs76 =     clhs12*clhs71;
        const double clhs77 =     NormalSlave(1,1)*(clhs58 + clhs73 + clhs74 - clhs75 - clhs76) + clhs72;
        const double clhs78 =     NormalSlave(1,0)*clhs77;
        const double clhs79 =     PenaltyParameter[1]*(clhs67 + clhs78);
        const double clhs80 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs81 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs82 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs83 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs84 =     NormalSlave(1,0)*(clhs17*clhs80 + clhs19*clhs81 - clhs21*clhs82 - clhs23*clhs83 + clhs5) + NormalSlave(1,1)*(-clhs10*clhs82 - clhs12*clhs83 + clhs4*clhs80 + clhs7*clhs81);
        const double clhs85 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs84 + clhs1*clhs29);
        const double clhs86 =     clhs28 - clhs33 - clhs35;
        const double clhs87 =     clhs29*clhs37;
        const double clhs88 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs89 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs90 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs91 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs92 =     NormalSlave(1,0)*(clhs17*clhs88 + clhs19*clhs89 - clhs21*clhs90 - clhs23*clhs91) + NormalSlave(1,1)*(-clhs10*clhs90 - clhs12*clhs91 + clhs4*clhs88 + clhs5 + clhs7*clhs89);
        const double clhs93 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs92 + clhs87);
        const double clhs94 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs95 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs96 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs97 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs98 =     NormalSlave(1,0)*(clhs17*clhs94 + clhs19*clhs95 - clhs21*clhs96 - clhs23*clhs97 + clhs8) + NormalSlave(1,1)*(-clhs10*clhs96 - clhs12*clhs97 + clhs4*clhs94 + clhs7*clhs95);
        const double clhs99 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs98 + clhs1*clhs30);
        const double clhs100 =     clhs30*clhs37;
        const double clhs101 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs102 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs103 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs104 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs105 =     NormalSlave(1,0)*(clhs101*clhs17 + clhs102*clhs19 - clhs103*clhs21 - clhs104*clhs23) + NormalSlave(1,1)*(-clhs10*clhs103 + clhs101*clhs4 + clhs102*clhs7 - clhs104*clhs12 + clhs8);
        const double clhs106 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs105 + clhs100);
        const double clhs107 =     DynamicFactor[1]*ScaleFactor;
        const double clhs108 =     clhs0*clhs107;
        const double clhs109 =     NormalSlave(1,1)*clhs25;
        const double clhs110 =     PenaltyParameter[1]*(clhs109 + clhs39);
        const double clhs111 =     LM(1,1)*ScaleFactor;
        const double clhs112 =     TangentSlaveXi(1,1)*clhs32;
        const double clhs113 =     NormalSlave(1,1)*clhs34;
        const double clhs114 =     -clhs111 + clhs112 + clhs113;
        const double clhs115 =     TangentFactor*std::pow(TangentSlaveXi(1,1), 2);
        const double clhs116 =     clhs115*clhs2;
        const double clhs117 =     NormalSlave(1,1)*clhs49;
        const double clhs118 =     PenaltyParameter[1]*(-clhs116 + clhs117);
        const double clhs119 =     NormalSlave(1,1)*clhs63;
        const double clhs120 =     PenaltyParameter[1]*(clhs119 + clhs67);
        const double clhs121 =     clhs115*clhs31;
        const double clhs122 =     NormalSlave(1,1)*clhs77;
        const double clhs123 =     PenaltyParameter[1]*(-clhs121 + clhs122);
        const double clhs124 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs84 + clhs87);
        const double clhs125 =     clhs111 - clhs112 - clhs113;
        const double clhs126 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs92 + clhs115*clhs29);
        const double clhs127 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs98 + clhs100);
        const double clhs128 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs105 + clhs115*clhs30);
        const double clhs129 =     clhs107*clhs13;
        const double clhs130 =     PenaltyParameter[1]*(-clhs26 + clhs3);
        const double clhs131 =     PenaltyParameter[1]*(clhs38 - clhs50);
        const double clhs132 =     PenaltyParameter[1]*(clhs52 - clhs64);
        const double clhs133 =     PenaltyParameter[1]*(clhs66 - clhs78);
        const double clhs134 =     -clhs107*clhs5;
        const double clhs135 =     PenaltyParameter[1]*(-clhs109 + clhs38);
        const double clhs136 =     PenaltyParameter[1]*(clhs116 - clhs117);
        const double clhs137 =     PenaltyParameter[1]*(-clhs119 + clhs66);
        const double clhs138 =     PenaltyParameter[1]*(clhs121 - clhs122);
        const double clhs139 =     -clhs107*clhs8;
        const double clhs140 =     1.0*TangentSlaveXi(1,0);
        const double clhs141 =     1.0*TangentSlaveXi(1,1);
    
        rLocalLHS(0,0)+=-DynamicFactor[1]*(clhs0*clhs27 + clhs11*clhs36);
        rLocalLHS(0,1)+=-DynamicFactor[1]*(clhs0*clhs51 + clhs36*clhs42);
        rLocalLHS(0,2)+=-DynamicFactor[1]*(clhs0*clhs65 + clhs36*clhs55);
        rLocalLHS(0,3)+=-DynamicFactor[1]*(clhs0*clhs79 + clhs36*clhs70);
        rLocalLHS(0,4)+=DynamicFactor[1]*(-clhs0*clhs85 + clhs82*clhs86);
        rLocalLHS(0,5)+=DynamicFactor[1]*(-clhs0*clhs93 + clhs86*clhs90);
        rLocalLHS(0,6)+=DynamicFactor[1]*(-clhs0*clhs99 + clhs86*clhs96);
        rLocalLHS(0,7)+=DynamicFactor[1]*(-clhs0*clhs106 + clhs103*clhs86);
        rLocalLHS(0,10)+=clhs108;
        rLocalLHS(1,0)+=-DynamicFactor[1]*(clhs0*clhs110 + clhs11*clhs114);
        rLocalLHS(1,1)+=-DynamicFactor[1]*(clhs0*clhs118 + clhs114*clhs42);
        rLocalLHS(1,2)+=-DynamicFactor[1]*(clhs0*clhs120 + clhs114*clhs55);
        rLocalLHS(1,3)+=-DynamicFactor[1]*(clhs0*clhs123 + clhs114*clhs70);
        rLocalLHS(1,4)+=DynamicFactor[1]*(-clhs0*clhs124 + clhs125*clhs82);
        rLocalLHS(1,5)+=DynamicFactor[1]*(-clhs0*clhs126 + clhs125*clhs90);
        rLocalLHS(1,6)+=DynamicFactor[1]*(-clhs0*clhs127 + clhs125*clhs96);
        rLocalLHS(1,7)+=DynamicFactor[1]*(-clhs0*clhs128 + clhs103*clhs125);
        rLocalLHS(1,11)+=clhs108;
        rLocalLHS(2,0)+=-DynamicFactor[1]*(clhs13*clhs27 + clhs14*clhs36);
        rLocalLHS(2,1)+=-DynamicFactor[1]*(clhs13*clhs51 + clhs36*clhs43);
        rLocalLHS(2,2)+=-DynamicFactor[1]*(clhs13*clhs65 + clhs36*clhs56);
        rLocalLHS(2,3)+=-DynamicFactor[1]*(clhs13*clhs79 + clhs36*clhs71);
        rLocalLHS(2,4)+=DynamicFactor[1]*(-clhs13*clhs85 + clhs83*clhs86);
        rLocalLHS(2,5)+=DynamicFactor[1]*(-clhs13*clhs93 + clhs86*clhs91);
        rLocalLHS(2,6)+=DynamicFactor[1]*(-clhs13*clhs99 + clhs86*clhs97);
        rLocalLHS(2,7)+=DynamicFactor[1]*(clhs104*clhs86 - clhs106*clhs13);
        rLocalLHS(2,10)+=clhs129;
        rLocalLHS(3,0)+=-DynamicFactor[1]*(clhs110*clhs13 + clhs114*clhs14);
        rLocalLHS(3,1)+=-DynamicFactor[1]*(clhs114*clhs43 + clhs118*clhs13);
        rLocalLHS(3,2)+=-DynamicFactor[1]*(clhs114*clhs56 + clhs120*clhs13);
        rLocalLHS(3,3)+=-DynamicFactor[1]*(clhs114*clhs71 + clhs123*clhs13);
        rLocalLHS(3,4)+=DynamicFactor[1]*(-clhs124*clhs13 + clhs125*clhs83);
        rLocalLHS(3,5)+=DynamicFactor[1]*(clhs125*clhs91 - clhs126*clhs13);
        rLocalLHS(3,6)+=DynamicFactor[1]*(clhs125*clhs97 - clhs127*clhs13);
        rLocalLHS(3,7)+=DynamicFactor[1]*(clhs104*clhs125 - clhs128*clhs13);
        rLocalLHS(3,11)+=clhs129;
        rLocalLHS(4,0)+=-DynamicFactor[1]*(clhs130*clhs5 + clhs6*clhs86);
        rLocalLHS(4,1)+=-DynamicFactor[1]*(clhs131*clhs5 + clhs40*clhs86);
        rLocalLHS(4,2)+=-DynamicFactor[1]*(clhs132*clhs5 + clhs53*clhs86);
        rLocalLHS(4,3)+=-DynamicFactor[1]*(clhs133*clhs5 + clhs68*clhs86);
        rLocalLHS(4,4)+=DynamicFactor[1]*(clhs5*clhs85 - clhs80*clhs86);
        rLocalLHS(4,5)+=DynamicFactor[1]*(clhs5*clhs93 - clhs86*clhs88);
        rLocalLHS(4,6)+=DynamicFactor[1]*(clhs5*clhs99 - clhs86*clhs94);
        rLocalLHS(4,7)+=DynamicFactor[1]*(-clhs101*clhs86 + clhs106*clhs5);
        rLocalLHS(4,10)+=clhs134;
        rLocalLHS(5,0)+=-DynamicFactor[1]*(clhs125*clhs6 + clhs135*clhs5);
        rLocalLHS(5,1)+=-DynamicFactor[1]*(clhs125*clhs40 + clhs136*clhs5);
        rLocalLHS(5,2)+=-DynamicFactor[1]*(clhs125*clhs53 + clhs137*clhs5);
        rLocalLHS(5,3)+=-DynamicFactor[1]*(clhs125*clhs68 + clhs138*clhs5);
        rLocalLHS(5,4)+=DynamicFactor[1]*(clhs124*clhs5 - clhs125*clhs80);
        rLocalLHS(5,5)+=DynamicFactor[1]*(-clhs125*clhs88 + clhs126*clhs5);
        rLocalLHS(5,6)+=DynamicFactor[1]*(-clhs125*clhs94 + clhs127*clhs5);
        rLocalLHS(5,7)+=DynamicFactor[1]*(-clhs101*clhs125 + clhs128*clhs5);
        rLocalLHS(5,11)+=clhs134;
        rLocalLHS(6,0)+=-DynamicFactor[1]*(clhs130*clhs8 + clhs86*clhs9);
        rLocalLHS(6,1)+=-DynamicFactor[1]*(clhs131*clhs8 + clhs41*clhs86);
        rLocalLHS(6,2)+=-DynamicFactor[1]*(clhs132*clhs8 + clhs54*clhs86);
        rLocalLHS(6,3)+=-DynamicFactor[1]*(clhs133*clhs8 + clhs69*clhs86);
        rLocalLHS(6,4)+=DynamicFactor[1]*(clhs8*clhs85 - clhs81*clhs86);
        rLocalLHS(6,5)+=DynamicFactor[1]*(clhs8*clhs93 - clhs86*clhs89);
        rLocalLHS(6,6)+=DynamicFactor[1]*(clhs8*clhs99 - clhs86*clhs95);
        rLocalLHS(6,7)+=DynamicFactor[1]*(-clhs102*clhs86 + clhs106*clhs8);
        rLocalLHS(6,10)+=clhs139;
        rLocalLHS(7,0)+=-DynamicFactor[1]*(clhs125*clhs9 + clhs135*clhs8);
        rLocalLHS(7,1)+=-DynamicFactor[1]*(clhs125*clhs41 + clhs136*clhs8);
        rLocalLHS(7,2)+=-DynamicFactor[1]*(clhs125*clhs54 + clhs137*clhs8);
        rLocalLHS(7,3)+=-DynamicFactor[1]*(clhs125*clhs69 + clhs138*clhs8);
        rLocalLHS(7,4)+=DynamicFactor[1]*(clhs124*clhs8 - clhs125*clhs81);
        rLocalLHS(7,5)+=DynamicFactor[1]*(-clhs125*clhs89 + clhs126*clhs8);
        rLocalLHS(7,6)+=DynamicFactor[1]*(-clhs125*clhs95 + clhs127*clhs8);
        rLocalLHS(7,7)+=DynamicFactor[1]*(-clhs102*clhs125 + clhs128*clhs8);
        rLocalLHS(7,11)+=clhs139;
        rLocalLHS(10,0)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs0 - clhs18 - clhs20 + clhs22 + clhs24) + clhs15);
        rLocalLHS(10,1)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs0 - clhs45 - clhs46 + clhs47 + clhs48) + clhs44);
        rLocalLHS(10,2)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs13 - clhs59 - clhs60 + clhs61 + clhs62) + clhs57);
        rLocalLHS(10,3)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs13 - clhs73 - clhs74 + clhs75 + clhs76) + clhs72);
        rLocalLHS(10,4)+=-ScaleFactor*clhs84;
        rLocalLHS(10,5)+=-ScaleFactor*clhs92;
        rLocalLHS(10,6)+=-ScaleFactor*clhs98;
        rLocalLHS(10,7)+=-ScaleFactor*clhs105;
        rLocalLHS(11,0)+=clhs140*clhs2;
        rLocalLHS(11,1)+=clhs141*clhs2;
        rLocalLHS(11,2)+=clhs140*clhs31;
        rLocalLHS(11,3)+=clhs141*clhs31;
        rLocalLHS(11,4)+=-clhs140*clhs29;
        rLocalLHS(11,5)+=-clhs141*clhs29;
        rLocalLHS(11,6)+=-clhs140*clhs30;
        rLocalLHS(11,7)+=-clhs141*clhs30;
        rLocalLHS(11,10)+=TangentSlaveXi(0,0);
        rLocalLHS(11,11)+=TangentSlaveXi(0,1);
    } else { // ACTIVE-STICK
        const double clhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs1 =     TangentFactor*std::pow(TangentSlaveXi(1,0), 2);
        const double clhs2 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double clhs3 =     clhs1*clhs2;
        const double clhs4 =     X1(0,1) + u1(0,1);
        const double clhs5 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs6 =     DeltaDOperator[4](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs7 =     X1(1,1) + u1(1,1);
        const double clhs8 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs9 =     DeltaDOperator[4](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs10 =     X2(0,1) + u2(0,1);
        const double clhs11 =     DeltaMOperator[4](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs12 =     X2(1,1) + u2(1,1);
        const double clhs13 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double clhs14 =     DeltaMOperator[4](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
        const double clhs15 =     NormalSlave(1,1)*(-clhs10*clhs11 - clhs12*clhs14 + clhs4*clhs6 + clhs7*clhs9);
        const double clhs16 =     -clhs0;
        const double clhs17 =     X1(0,0) + u1(0,0);
        const double clhs18 =     clhs17*clhs6;
        const double clhs19 =     X1(1,0) + u1(1,0);
        const double clhs20 =     clhs19*clhs9;
        const double clhs21 =     X2(0,0) + u2(0,0);
        const double clhs22 =     clhs11*clhs21;
        const double clhs23 =     X2(1,0) + u2(1,0);
        const double clhs24 =     clhs14*clhs23;
        const double clhs25 =     NormalSlave(1,0)*(clhs16 + clhs18 + clhs20 - clhs22 - clhs24) + clhs15;
        const double clhs26 =     NormalSlave(1,0)*clhs25;
        const double clhs27 =     PenaltyParameter[1]*(clhs26 - clhs3);
        const double clhs28 =     LM(1,0)*ScaleFactor;
        const double clhs29 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double clhs30 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double clhs31 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double clhs32 =     TangentSlaveXi(1,0)*(clhs17*clhs29 + clhs19*clhs30 - clhs2*clhs21 - clhs23*clhs31) + TangentSlaveXi(1,1)*(-clhs10*clhs2 - clhs12*clhs31 + clhs29*clhs4 + clhs30*clhs7);
        const double clhs33 =     PenaltyParameter[1]*TangentFactor*clhs32;
        const double clhs34 =     TangentSlaveXi(1,0)*clhs33;
        const double clhs35 =     PenaltyParameter[1]*(NormalSlave(1,0)*(-clhs0*clhs21 - clhs13*clhs23 + clhs17*clhs5 + clhs19*clhs8) + NormalSlave(1,1)*(-clhs0*clhs10 - clhs12*clhs13 + clhs4*clhs5 + clhs7*clhs8));
        const double clhs36 =     NormalSlave(1,0)*clhs35;
        const double clhs37 =     -clhs28 + clhs34 + clhs36;
        const double clhs38 =     TangentFactor*TangentSlaveXi(1,0)*TangentSlaveXi(1,1);
        const double clhs39 =     clhs2*clhs38;
        const double clhs40 =     -clhs39;
        const double clhs41 =     DeltaDOperator[5](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs42 =     DeltaDOperator[5](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs43 =     DeltaMOperator[5](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs44 =     DeltaMOperator[5](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
        const double clhs45 =     NormalSlave(1,0)*(clhs17*clhs41 + clhs19*clhs42 - clhs21*clhs43 - clhs23*clhs44);
        const double clhs46 =     clhs4*clhs41;
        const double clhs47 =     clhs42*clhs7;
        const double clhs48 =     clhs10*clhs43;
        const double clhs49 =     clhs12*clhs44;
        const double clhs50 =     NormalSlave(1,1)*(clhs16 + clhs46 + clhs47 - clhs48 - clhs49) + clhs45;
        const double clhs51 =     NormalSlave(1,0)*clhs50;
        const double clhs52 =     PenaltyParameter[1]*(clhs40 + clhs51);
        const double clhs53 =     clhs1*clhs31;
        const double clhs54 =     DeltaDOperator[6](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs55 =     DeltaDOperator[6](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs56 =     DeltaMOperator[6](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs57 =     DeltaMOperator[6](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
        const double clhs58 =     NormalSlave(1,1)*(-clhs10*clhs56 - clhs12*clhs57 + clhs4*clhs54 + clhs55*clhs7);
        const double clhs59 =     -clhs13;
        const double clhs60 =     clhs17*clhs54;
        const double clhs61 =     clhs19*clhs55;
        const double clhs62 =     clhs21*clhs56;
        const double clhs63 =     clhs23*clhs57;
        const double clhs64 =     NormalSlave(1,0)*(clhs59 + clhs60 + clhs61 - clhs62 - clhs63) + clhs58;
        const double clhs65 =     NormalSlave(1,0)*clhs64;
        const double clhs66 =     PenaltyParameter[1]*(-clhs53 + clhs65);
        const double clhs67 =     clhs31*clhs38;
        const double clhs68 =     -clhs67;
        const double clhs69 =     DeltaDOperator[7](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs70 =     DeltaDOperator[7](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs71 =     DeltaMOperator[7](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs72 =     DeltaMOperator[7](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
        const double clhs73 =     NormalSlave(1,0)*(clhs17*clhs69 + clhs19*clhs70 - clhs21*clhs71 - clhs23*clhs72);
        const double clhs74 =     clhs4*clhs69;
        const double clhs75 =     clhs7*clhs70;
        const double clhs76 =     clhs10*clhs71;
        const double clhs77 =     clhs12*clhs72;
        const double clhs78 =     NormalSlave(1,1)*(clhs59 + clhs74 + clhs75 - clhs76 - clhs77) + clhs73;
        const double clhs79 =     NormalSlave(1,0)*clhs78;
        const double clhs80 =     PenaltyParameter[1]*(clhs68 + clhs79);
        const double clhs81 =     DeltaDOperator[0](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs82 =     DeltaDOperator[0](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs83 =     DeltaMOperator[0](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs84 =     DeltaMOperator[0](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
        const double clhs85 =     NormalSlave(1,0)*(clhs17*clhs81 + clhs19*clhs82 - clhs21*clhs83 - clhs23*clhs84 + clhs5) + NormalSlave(1,1)*(-clhs10*clhs83 - clhs12*clhs84 + clhs4*clhs81 + clhs7*clhs82);
        const double clhs86 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs85 + clhs1*clhs29);
        const double clhs87 =     clhs28 - clhs34 - clhs36;
        const double clhs88 =     clhs29*clhs38;
        const double clhs89 =     DeltaDOperator[1](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs90 =     DeltaDOperator[1](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs91 =     DeltaMOperator[1](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs92 =     DeltaMOperator[1](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
        const double clhs93 =     NormalSlave(1,0)*(clhs17*clhs89 + clhs19*clhs90 - clhs21*clhs91 - clhs23*clhs92) + NormalSlave(1,1)*(-clhs10*clhs91 - clhs12*clhs92 + clhs4*clhs89 + clhs5 + clhs7*clhs90);
        const double clhs94 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs93 + clhs88);
        const double clhs95 =     DeltaDOperator[2](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs96 =     DeltaDOperator[2](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs97 =     DeltaMOperator[2](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs98 =     DeltaMOperator[2](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
        const double clhs99 =     NormalSlave(1,0)*(clhs17*clhs95 + clhs19*clhs96 - clhs21*clhs97 - clhs23*clhs98 + clhs8) + NormalSlave(1,1)*(-clhs10*clhs97 - clhs12*clhs98 + clhs4*clhs95 + clhs7*clhs96);
        const double clhs100 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs99 + clhs1*clhs30);
        const double clhs101 =     clhs30*clhs38;
        const double clhs102 =     DeltaDOperator[3](1,0); // DERIVATIVE(DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs103 =     DeltaDOperator[3](1,1); // DERIVATIVE(DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs104 =     DeltaMOperator[3](1,0); // DERIVATIVE(MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs105 =     DeltaMOperator[3](1,1); // DERIVATIVE(MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
        const double clhs106 =     NormalSlave(1,0)*(clhs102*clhs17 + clhs103*clhs19 - clhs104*clhs21 - clhs105*clhs23) + NormalSlave(1,1)*(-clhs10*clhs104 + clhs102*clhs4 + clhs103*clhs7 - clhs105*clhs12 + clhs8);
        const double clhs107 =     PenaltyParameter[1]*(NormalSlave(1,0)*clhs106 + clhs101);
        const double clhs108 =     DynamicFactor[1]*ScaleFactor;
        const double clhs109 =     clhs0*clhs108;
        const double clhs110 =     NormalSlave(1,1)*clhs25;
        const double clhs111 =     PenaltyParameter[1]*(clhs110 + clhs40);
        const double clhs112 =     LM(1,1)*ScaleFactor;
        const double clhs113 =     TangentSlaveXi(1,1)*clhs33;
        const double clhs114 =     NormalSlave(1,1)*clhs35;
        const double clhs115 =     -clhs112 + clhs113 + clhs114;
        const double clhs116 =     TangentFactor*std::pow(TangentSlaveXi(1,1), 2);
        const double clhs117 =     clhs116*clhs2;
        const double clhs118 =     NormalSlave(1,1)*clhs50;
        const double clhs119 =     PenaltyParameter[1]*(-clhs117 + clhs118);
        const double clhs120 =     NormalSlave(1,1)*clhs64;
        const double clhs121 =     PenaltyParameter[1]*(clhs120 + clhs68);
        const double clhs122 =     clhs116*clhs31;
        const double clhs123 =     NormalSlave(1,1)*clhs78;
        const double clhs124 =     PenaltyParameter[1]*(-clhs122 + clhs123);
        const double clhs125 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs85 + clhs88);
        const double clhs126 =     clhs112 - clhs113 - clhs114;
        const double clhs127 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs93 + clhs116*clhs29);
        const double clhs128 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs99 + clhs101);
        const double clhs129 =     PenaltyParameter[1]*(NormalSlave(1,1)*clhs106 + clhs116*clhs30);
        const double clhs130 =     clhs108*clhs13;
        const double clhs131 =     PenaltyParameter[1]*(-clhs26 + clhs3);
        const double clhs132 =     PenaltyParameter[1]*(clhs39 - clhs51);
        const double clhs133 =     PenaltyParameter[1]*(clhs53 - clhs65);
        const double clhs134 =     PenaltyParameter[1]*(clhs67 - clhs79);
        const double clhs135 =     -clhs108*clhs5;
        const double clhs136 =     PenaltyParameter[1]*(-clhs110 + clhs39);
        const double clhs137 =     PenaltyParameter[1]*(clhs117 - clhs118);
        const double clhs138 =     PenaltyParameter[1]*(-clhs120 + clhs67);
        const double clhs139 =     PenaltyParameter[1]*(clhs122 - clhs123);
        const double clhs140 =     -clhs108*clhs8;
        const double clhs141 =     PenaltyParameter[1]*TangentFactor*mu[1];
        const double clhs142 =     ScaleFactor*(LM(1,0)*NormalSlave(1,0) + LM(1,1)*NormalSlave(1,1));
        const double clhs143 =     clhs142 - clhs35;
        const double clhs144 =     TangentSlaveXi(1,0)*clhs143;
        const double clhs145 =     PenaltyParameter[1]*clhs32;
        const double clhs146 =     TangentSlaveXi(1,1)*clhs143;
        const double clhs147 =     -clhs142 + clhs35;
        const double clhs148 =     TangentSlaveXi(1,0)*clhs147;
        const double clhs149 =     TangentSlaveXi(1,1)*clhs147;
    
        rLocalLHS(0,0)+=-DynamicFactor[1]*(clhs0*clhs27 + clhs11*clhs37);
        rLocalLHS(0,1)+=-DynamicFactor[1]*(clhs0*clhs52 + clhs37*clhs43);
        rLocalLHS(0,2)+=-DynamicFactor[1]*(clhs0*clhs66 + clhs37*clhs56);
        rLocalLHS(0,3)+=-DynamicFactor[1]*(clhs0*clhs80 + clhs37*clhs71);
        rLocalLHS(0,4)+=DynamicFactor[1]*(-clhs0*clhs86 + clhs83*clhs87);
        rLocalLHS(0,5)+=DynamicFactor[1]*(-clhs0*clhs94 + clhs87*clhs91);
        rLocalLHS(0,6)+=DynamicFactor[1]*(-clhs0*clhs100 + clhs87*clhs97);
        rLocalLHS(0,7)+=DynamicFactor[1]*(-clhs0*clhs107 + clhs104*clhs87);
        rLocalLHS(0,10)+=clhs109;
        rLocalLHS(1,0)+=-DynamicFactor[1]*(clhs0*clhs111 + clhs11*clhs115);
        rLocalLHS(1,1)+=-DynamicFactor[1]*(clhs0*clhs119 + clhs115*clhs43);
        rLocalLHS(1,2)+=-DynamicFactor[1]*(clhs0*clhs121 + clhs115*clhs56);
        rLocalLHS(1,3)+=-DynamicFactor[1]*(clhs0*clhs124 + clhs115*clhs71);
        rLocalLHS(1,4)+=DynamicFactor[1]*(-clhs0*clhs125 + clhs126*clhs83);
        rLocalLHS(1,5)+=DynamicFactor[1]*(-clhs0*clhs127 + clhs126*clhs91);
        rLocalLHS(1,6)+=DynamicFactor[1]*(-clhs0*clhs128 + clhs126*clhs97);
        rLocalLHS(1,7)+=DynamicFactor[1]*(-clhs0*clhs129 + clhs104*clhs126);
        rLocalLHS(1,11)+=clhs109;
        rLocalLHS(2,0)+=-DynamicFactor[1]*(clhs13*clhs27 + clhs14*clhs37);
        rLocalLHS(2,1)+=-DynamicFactor[1]*(clhs13*clhs52 + clhs37*clhs44);
        rLocalLHS(2,2)+=-DynamicFactor[1]*(clhs13*clhs66 + clhs37*clhs57);
        rLocalLHS(2,3)+=-DynamicFactor[1]*(clhs13*clhs80 + clhs37*clhs72);
        rLocalLHS(2,4)+=DynamicFactor[1]*(-clhs13*clhs86 + clhs84*clhs87);
        rLocalLHS(2,5)+=DynamicFactor[1]*(-clhs13*clhs94 + clhs87*clhs92);
        rLocalLHS(2,6)+=DynamicFactor[1]*(-clhs100*clhs13 + clhs87*clhs98);
        rLocalLHS(2,7)+=DynamicFactor[1]*(clhs105*clhs87 - clhs107*clhs13);
        rLocalLHS(2,10)+=clhs130;
        rLocalLHS(3,0)+=-DynamicFactor[1]*(clhs111*clhs13 + clhs115*clhs14);
        rLocalLHS(3,1)+=-DynamicFactor[1]*(clhs115*clhs44 + clhs119*clhs13);
        rLocalLHS(3,2)+=-DynamicFactor[1]*(clhs115*clhs57 + clhs121*clhs13);
        rLocalLHS(3,3)+=-DynamicFactor[1]*(clhs115*clhs72 + clhs124*clhs13);
        rLocalLHS(3,4)+=DynamicFactor[1]*(-clhs125*clhs13 + clhs126*clhs84);
        rLocalLHS(3,5)+=DynamicFactor[1]*(clhs126*clhs92 - clhs127*clhs13);
        rLocalLHS(3,6)+=DynamicFactor[1]*(clhs126*clhs98 - clhs128*clhs13);
        rLocalLHS(3,7)+=DynamicFactor[1]*(clhs105*clhs126 - clhs129*clhs13);
        rLocalLHS(3,11)+=clhs130;
        rLocalLHS(4,0)+=-DynamicFactor[1]*(clhs131*clhs5 + clhs6*clhs87);
        rLocalLHS(4,1)+=-DynamicFactor[1]*(clhs132*clhs5 + clhs41*clhs87);
        rLocalLHS(4,2)+=-DynamicFactor[1]*(clhs133*clhs5 + clhs54*clhs87);
        rLocalLHS(4,3)+=-DynamicFactor[1]*(clhs134*clhs5 + clhs69*clhs87);
        rLocalLHS(4,4)+=DynamicFactor[1]*(clhs5*clhs86 - clhs81*clhs87);
        rLocalLHS(4,5)+=DynamicFactor[1]*(clhs5*clhs94 - clhs87*clhs89);
        rLocalLHS(4,6)+=DynamicFactor[1]*(clhs100*clhs5 - clhs87*clhs95);
        rLocalLHS(4,7)+=DynamicFactor[1]*(-clhs102*clhs87 + clhs107*clhs5);
        rLocalLHS(4,10)+=clhs135;
        rLocalLHS(5,0)+=-DynamicFactor[1]*(clhs126*clhs6 + clhs136*clhs5);
        rLocalLHS(5,1)+=-DynamicFactor[1]*(clhs126*clhs41 + clhs137*clhs5);
        rLocalLHS(5,2)+=-DynamicFactor[1]*(clhs126*clhs54 + clhs138*clhs5);
        rLocalLHS(5,3)+=-DynamicFactor[1]*(clhs126*clhs69 + clhs139*clhs5);
        rLocalLHS(5,4)+=DynamicFactor[1]*(clhs125*clhs5 - clhs126*clhs81);
        rLocalLHS(5,5)+=DynamicFactor[1]*(-clhs126*clhs89 + clhs127*clhs5);
        rLocalLHS(5,6)+=DynamicFactor[1]*(-clhs126*clhs95 + clhs128*clhs5);
        rLocalLHS(5,7)+=DynamicFactor[1]*(-clhs102*clhs126 + clhs129*clhs5);
        rLocalLHS(5,11)+=clhs135;
        rLocalLHS(6,0)+=-DynamicFactor[1]*(clhs131*clhs8 + clhs87*clhs9);
        rLocalLHS(6,1)+=-DynamicFactor[1]*(clhs132*clhs8 + clhs42*clhs87);
        rLocalLHS(6,2)+=-DynamicFactor[1]*(clhs133*clhs8 + clhs55*clhs87);
        rLocalLHS(6,3)+=-DynamicFactor[1]*(clhs134*clhs8 + clhs70*clhs87);
        rLocalLHS(6,4)+=DynamicFactor[1]*(clhs8*clhs86 - clhs82*clhs87);
        rLocalLHS(6,5)+=DynamicFactor[1]*(clhs8*clhs94 - clhs87*clhs90);
        rLocalLHS(6,6)+=DynamicFactor[1]*(clhs100*clhs8 - clhs87*clhs96);
        rLocalLHS(6,7)+=DynamicFactor[1]*(-clhs103*clhs87 + clhs107*clhs8);
        rLocalLHS(6,10)+=clhs140;
        rLocalLHS(7,0)+=-DynamicFactor[1]*(clhs126*clhs9 + clhs136*clhs8);
        rLocalLHS(7,1)+=-DynamicFactor[1]*(clhs126*clhs42 + clhs137*clhs8);
        rLocalLHS(7,2)+=-DynamicFactor[1]*(clhs126*clhs55 + clhs138*clhs8);
        rLocalLHS(7,3)+=-DynamicFactor[1]*(clhs126*clhs70 + clhs139*clhs8);
        rLocalLHS(7,4)+=DynamicFactor[1]*(clhs125*clhs8 - clhs126*clhs82);
        rLocalLHS(7,5)+=DynamicFactor[1]*(-clhs126*clhs90 + clhs127*clhs8);
        rLocalLHS(7,6)+=DynamicFactor[1]*(-clhs126*clhs96 + clhs128*clhs8);
        rLocalLHS(7,7)+=DynamicFactor[1]*(-clhs103*clhs126 + clhs129*clhs8);
        rLocalLHS(7,11)+=clhs140;
        rLocalLHS(10,0)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs0 - clhs18 - clhs20 + clhs22 + clhs24) + clhs15);
        rLocalLHS(10,1)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs0 - clhs46 - clhs47 + clhs48 + clhs49) + clhs45);
        rLocalLHS(10,2)+=-ScaleFactor*(-NormalSlave(1,0)*(clhs13 - clhs60 - clhs61 + clhs62 + clhs63) + clhs58);
        rLocalLHS(10,3)+=-ScaleFactor*(-NormalSlave(1,1)*(clhs13 - clhs74 - clhs75 + clhs76 + clhs77) + clhs73);
        rLocalLHS(10,4)+=-ScaleFactor*clhs85;
        rLocalLHS(10,5)+=-ScaleFactor*clhs93;
        rLocalLHS(10,6)+=-ScaleFactor*clhs99;
        rLocalLHS(10,7)+=-ScaleFactor*clhs106;
        rLocalLHS(11,0)+=-clhs141*(clhs144*clhs2 + clhs145*clhs25);
        rLocalLHS(11,1)+=-clhs141*(clhs145*clhs50 + clhs146*clhs2);
        rLocalLHS(11,2)+=-clhs141*(clhs144*clhs31 + clhs145*clhs64);
        rLocalLHS(11,3)+=-clhs141*(clhs145*clhs78 + clhs146*clhs31);
        rLocalLHS(11,4)+=-clhs141*(clhs145*clhs85 + clhs148*clhs29);
        rLocalLHS(11,5)+=-clhs141*(clhs145*clhs93 + clhs149*clhs29);
        rLocalLHS(11,6)+=-clhs141*(clhs145*clhs99 + clhs148*clhs30);
        rLocalLHS(11,7)+=-clhs141*(clhs106*clhs145 + clhs149*clhs30);
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
        const double crhs1 =     X1(0,0) + u1(0,0);
        const double crhs2 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double crhs3 =     X1(1,0) + u1(1,0);
        const double crhs4 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double crhs5 =     X2(0,0) + u2(0,0);
        const double crhs6 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double crhs7 =     X2(1,0) + u2(1,0);
        const double crhs8 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double crhs9 =     X1(0,1) + u1(0,1);
        const double crhs10 =     X1(1,1) + u1(1,1);
        const double crhs11 =     X2(0,1) + u2(0,1);
        const double crhs12 =     X2(1,1) + u2(1,1);
        const double crhs13 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(crhs1*crhs2 + crhs3*crhs4 - crhs5*crhs6 - crhs7*crhs8) + TangentSlaveXi(0,1)*(crhs10*crhs4 - crhs11*crhs6 - crhs12*crhs8 + crhs2*crhs9));
        const double crhs14 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs15 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs16 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs17 =     NormalSlave(0,0)*(-crhs0*crhs5 + crhs1*crhs14 + crhs15*crhs3 - crhs16*crhs7) + NormalSlave(0,1)*(-crhs0*crhs11 + crhs10*crhs15 - crhs12*crhs16 + crhs14*crhs9);
        const double crhs18 =     PenaltyParameter[0]*crhs17;
        const double crhs19 =     DynamicFactor[0]*(-LM(0,0)*ScaleFactor + NormalSlave(0,0)*crhs18 + TangentSlaveXi(0,0)*crhs13);
        const double crhs20 =     DynamicFactor[0]*(-LM(0,1)*ScaleFactor + NormalSlave(0,1)*crhs18 + TangentSlaveXi(0,1)*crhs13);
    
        rLocalRHS[0]+=crhs0*crhs19;
        rLocalRHS[1]+=crhs0*crhs20;
        rLocalRHS[2]+=crhs16*crhs19;
        rLocalRHS[3]+=crhs16*crhs20;
        rLocalRHS[4]+=-crhs14*crhs19;
        rLocalRHS[5]+=-crhs14*crhs20;
        rLocalRHS[6]+=-crhs15*crhs19;
        rLocalRHS[7]+=-crhs15*crhs20;
        rLocalRHS[8]+=ScaleFactor*crhs17;
        rLocalRHS[9]+=-1.0*(ScaleFactor*(LM(0,0)*TangentSlaveXi(0,0) + LM(0,1)*TangentSlaveXi(0,1)) - crhs13)/(PenaltyParameter[0]*TangentFactor);
    } else { // ACTIVE-STICK
        const double crhs0 =     MOperator(0,0); // MOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     X1(0,0) + u1(0,0);
        const double crhs2 =     StandardDOperator(0,0) - StandardDOperatorold(0,0);
        const double crhs3 =     X1(1,0) + u1(1,0);
        const double crhs4 =     StandardDOperator(0,1) - StandardDOperatorold(0,1);
        const double crhs5 =     X2(0,0) + u2(0,0);
        const double crhs6 =     StandardMOperator(0,0) - StandardMOperatorold(0,0);
        const double crhs7 =     X2(1,0) + u2(1,0);
        const double crhs8 =     StandardMOperator(0,1) - StandardMOperatorold(0,1);
        const double crhs9 =     X1(0,1) + u1(0,1);
        const double crhs10 =     X1(1,1) + u1(1,1);
        const double crhs11 =     X2(0,1) + u2(0,1);
        const double crhs12 =     X2(1,1) + u2(1,1);
        const double crhs13 =     PenaltyParameter[0]*TangentFactor*(TangentSlaveXi(0,0)*(crhs1*crhs2 + crhs3*crhs4 - crhs5*crhs6 - crhs7*crhs8) + TangentSlaveXi(0,1)*(crhs10*crhs4 - crhs11*crhs6 - crhs12*crhs8 + crhs2*crhs9));
        const double crhs14 =     DOperator(0,0); // DOPERATOR(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs15 =     DOperator(0,1); // DOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs16 =     MOperator(0,1); // MOPERATOR(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs17 =     NormalSlave(0,0)*(-crhs0*crhs5 + crhs1*crhs14 + crhs15*crhs3 - crhs16*crhs7) + NormalSlave(0,1)*(-crhs0*crhs11 + crhs10*crhs15 - crhs12*crhs16 + crhs14*crhs9);
        const double crhs18 =     PenaltyParameter[0]*crhs17;
        const double crhs19 =     DynamicFactor[0]*(-LM(0,0)*ScaleFactor + NormalSlave(0,0)*crhs18 + TangentSlaveXi(0,0)*crhs13);
        const double crhs20 =     DynamicFactor[0]*(-LM(0,1)*ScaleFactor + NormalSlave(0,1)*crhs18 + TangentSlaveXi(0,1)*crhs13);
    
        rLocalRHS[0]+=crhs0*crhs19;
        rLocalRHS[1]+=crhs0*crhs20;
        rLocalRHS[2]+=crhs16*crhs19;
        rLocalRHS[3]+=crhs16*crhs20;
        rLocalRHS[4]+=-crhs14*crhs19;
        rLocalRHS[5]+=-crhs14*crhs20;
        rLocalRHS[6]+=-crhs15*crhs19;
        rLocalRHS[7]+=-crhs15*crhs20;
        rLocalRHS[8]+=ScaleFactor*crhs17;
        rLocalRHS[9]+=crhs13*mu[0]*(-ScaleFactor*(LM(0,0)*NormalSlave(0,0) + LM(0,1)*NormalSlave(0,1)) + crhs18);
    }
    
    // NODE 1
    if (geometry[1].IsNot(ACTIVE)) { // INACTIVE
        const double crhs0 =     std::pow(ScaleFactor, 2)/PenaltyParameter[1];
    
        rLocalRHS[10]+=-crhs0*(LM(1,0)*NormalSlave(1,0) + LM(1,1)*NormalSlave(1,1));
        rLocalRHS[11]+=-crhs0*(LM(1,0)*TangentSlaveXi(1,0) + LM(1,1)*TangentSlaveXi(1,1));
    } else if (geometry[1].Is(SLIP)) { // ACTIVE-SLIP
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     X1(0,0) + u1(0,0);
        const double crhs2 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double crhs3 =     X1(1,0) + u1(1,0);
        const double crhs4 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double crhs5 =     X2(0,0) + u2(0,0);
        const double crhs6 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double crhs7 =     X2(1,0) + u2(1,0);
        const double crhs8 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double crhs9 =     X1(0,1) + u1(0,1);
        const double crhs10 =     X1(1,1) + u1(1,1);
        const double crhs11 =     X2(0,1) + u2(0,1);
        const double crhs12 =     X2(1,1) + u2(1,1);
        const double crhs13 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(crhs1*crhs2 + crhs3*crhs4 - crhs5*crhs6 - crhs7*crhs8) + TangentSlaveXi(1,1)*(crhs10*crhs4 - crhs11*crhs6 - crhs12*crhs8 + crhs2*crhs9));
        const double crhs14 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs15 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs16 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs17 =     NormalSlave(1,0)*(-crhs0*crhs5 + crhs1*crhs14 + crhs15*crhs3 - crhs16*crhs7) + NormalSlave(1,1)*(-crhs0*crhs11 + crhs10*crhs15 - crhs12*crhs16 + crhs14*crhs9);
        const double crhs18 =     PenaltyParameter[1]*crhs17;
        const double crhs19 =     DynamicFactor[1]*(-LM(1,0)*ScaleFactor + NormalSlave(1,0)*crhs18 + TangentSlaveXi(1,0)*crhs13);
        const double crhs20 =     DynamicFactor[1]*(-LM(1,1)*ScaleFactor + NormalSlave(1,1)*crhs18 + TangentSlaveXi(1,1)*crhs13);
    
        rLocalRHS[0]+=crhs0*crhs19;
        rLocalRHS[1]+=crhs0*crhs20;
        rLocalRHS[2]+=crhs16*crhs19;
        rLocalRHS[3]+=crhs16*crhs20;
        rLocalRHS[4]+=-crhs14*crhs19;
        rLocalRHS[5]+=-crhs14*crhs20;
        rLocalRHS[6]+=-crhs15*crhs19;
        rLocalRHS[7]+=-crhs15*crhs20;
        rLocalRHS[10]+=ScaleFactor*crhs17;
        rLocalRHS[11]+=-1.0*(ScaleFactor*(LM(1,0)*TangentSlaveXi(1,0) + LM(1,1)*TangentSlaveXi(1,1)) - crhs13)/(PenaltyParameter[1]*TangentFactor);
    } else { // ACTIVE-STICK
        const double crhs0 =     MOperator(1,0); // MOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs1 =     X1(0,0) + u1(0,0);
        const double crhs2 =     StandardDOperator(1,0) - StandardDOperatorold(1,0);
        const double crhs3 =     X1(1,0) + u1(1,0);
        const double crhs4 =     StandardDOperator(1,1) - StandardDOperatorold(1,1);
        const double crhs5 =     X2(0,0) + u2(0,0);
        const double crhs6 =     StandardMOperator(1,0) - StandardMOperatorold(1,0);
        const double crhs7 =     X2(1,0) + u2(1,0);
        const double crhs8 =     StandardMOperator(1,1) - StandardMOperatorold(1,1);
        const double crhs9 =     X1(0,1) + u1(0,1);
        const double crhs10 =     X1(1,1) + u1(1,1);
        const double crhs11 =     X2(0,1) + u2(0,1);
        const double crhs12 =     X2(1,1) + u2(1,1);
        const double crhs13 =     PenaltyParameter[1]*TangentFactor*(TangentSlaveXi(1,0)*(crhs1*crhs2 + crhs3*crhs4 - crhs5*crhs6 - crhs7*crhs8) + TangentSlaveXi(1,1)*(crhs10*crhs4 - crhs11*crhs6 - crhs12*crhs8 + crhs2*crhs9));
        const double crhs14 =     DOperator(1,0); // DOPERATOR(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs15 =     DOperator(1,1); // DOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs16 =     MOperator(1,1); // MOPERATOR(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
        const double crhs17 =     NormalSlave(1,0)*(-crhs0*crhs5 + crhs1*crhs14 + crhs15*crhs3 - crhs16*crhs7) + NormalSlave(1,1)*(-crhs0*crhs11 + crhs10*crhs15 - crhs12*crhs16 + crhs14*crhs9);
        const double crhs18 =     PenaltyParameter[1]*crhs17;
        const double crhs19 =     DynamicFactor[1]*(-LM(1,0)*ScaleFactor + NormalSlave(1,0)*crhs18 + TangentSlaveXi(1,0)*crhs13);
        const double crhs20 =     DynamicFactor[1]*(-LM(1,1)*ScaleFactor + NormalSlave(1,1)*crhs18 + TangentSlaveXi(1,1)*crhs13);
    
        rLocalRHS[0]+=crhs0*crhs19;
        rLocalRHS[1]+=crhs0*crhs20;
        rLocalRHS[2]+=crhs16*crhs19;
        rLocalRHS[3]+=crhs16*crhs20;
        rLocalRHS[4]+=-crhs14*crhs19;
        rLocalRHS[5]+=-crhs14*crhs20;
        rLocalRHS[6]+=-crhs15*crhs19;
        rLocalRHS[7]+=-crhs15*crhs20;
        rLocalRHS[10]+=ScaleFactor*crhs17;
        rLocalRHS[11]+=crhs13*mu[1]*(-ScaleFactor*(LM(1,0)*NormalSlave(1,0) + LM(1,1)*NormalSlave(1,1)) + crhs18);
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
