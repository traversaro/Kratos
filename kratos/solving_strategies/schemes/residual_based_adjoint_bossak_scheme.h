//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela, https://github.com/jcotela
//

#if !defined(KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <unordered_set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "utilities/indirect_scalar.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace Internals
{
    class AdjointBossakSchemeUtils
    {
        struct Hash
        {
            std::size_t operator()(const VariableData* const& p) const
            {
                return p->Key();
            }
        };

        struct Pred
        {
            bool operator()(const VariableData* const l, const VariableData* const r) const
            {
                return *l == *r;
            }
        };

        public:
            using TGetAdjointVector = Variable<Function<void(std::vector<VariableData const*>&)>>;

            template <typename TContainerType>
            std::vector<const VariableData*> Execute(const TGetAdjointVector& rFunc, const TContainerType& rElements)
            {
                KRATOS_TRY;
                const int number_of_elements = rElements.size();
                const int num_threads = OpenMPUtils::GetNumThreads();
                std::vector<const VariableData*> var_vec;
                std::vector<std::unordered_set<const VariableData*, Hash, Pred>> var_sets(num_threads);
#pragma omp parallel for private(var_vec)
                for (int i = 0; i < number_of_elements; ++i)
                {
                    auto& r_element = *(rElements.begin() + i);
                    r_element.GetValue(rFunc)(var_vec);
                    const int k = OpenMPUtils::ThisThread();
                    var_sets[k].insert(var_vec.begin(), var_vec.end());
                }
                std::unordered_set<const VariableData*, Hash, Pred> total_var_set;
                for (int i = 0; i < num_threads; ++i)
                {
                    total_var_set.insert(var_sets[i].begin(), var_sets[i].end());
                }
                return std::vector<const VariableData*>{total_var_set.begin(), total_var_set.end()};
                KRATOS_CATCH("");
            }

            void SetToZero_AdjointVars(const std::vector<const VariableData*>& rVariables, ModelPart::NodesContainerType& rNodes)
            {
                KRATOS_TRY;
                for (auto p_variable_data : rVariables)
                {
                    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(
                            p_variable_data->Name()))
                    {
                        const auto& r_variable =
                            KratosComponents<Variable<array_1d<double, 3>>>::Get(
                                p_variable_data->Name());
                        VariableUtils().SetToZero_VectorVar(r_variable, rNodes);
                    }
                    else if (KratosComponents<Variable<double>>::Has(
                                 p_variable_data->Name()))
                    {
                        const auto& r_variable = KratosComponents<Variable<double>>::Get(
                            p_variable_data->Name());
                        VariableUtils().SetToZero_ScalarVar(r_variable, rNodes);
                    }
                    else
                    {
                        KRATOS_ERROR << "Variable \"" << p_variable_data->Name()
                                     << "\" not found!\n";
                    }
                }
                KRATOS_CATCH("");
            }

            void Assemble_AdjointVars(const std::vector<const VariableData*>& rVariables,
                                      Communicator& rComm)
            {
                KRATOS_TRY;
                for (auto p_variable_data : rVariables)
                {
                    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(
                            p_variable_data->Name()))
                    {
                        const auto& r_variable =
                            KratosComponents<Variable<array_1d<double, 3>>>::Get(
                                p_variable_data->Name());
                        rComm.AssembleCurrentData(r_variable);
                    }
                    else if (KratosComponents<Variable<double>>::Has(
                                 p_variable_data->Name()))
                    {
                        const auto& r_variable = KratosComponents<Variable<double>>::Get(
                            p_variable_data->Name());
                        rComm.AssembleCurrentData(r_variable);
                    }
                    else
                    {
                        KRATOS_ERROR << "Variable \"" << p_variable_data->Name()
                                     << "\" not found!\n";
                    }
                }
                KRATOS_CATCH("");
            }
    };
}

/// A scheme for dynamic adjoint equations, using Bossak time integration.
/**
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedAdjointBossakScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedAdjointBossakScheme);

    typedef ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::SystemVectorType SystemVectorType;
    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResidualBasedAdjointBossakScheme(Parameters& rParameters, AdjointResponseFunction::Pointer pResponseFunction):
        ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction)
    {

        Parameters default_parameters(R"({
            "scheme_type": "bossak",
            "alpha_bossak": -0.3
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();
        mBetaNewmark = 0.25 * (1.0 - mAlphaBossak) * (1.0 - mAlphaBossak);
        mGammaNewmark = 0.5 - mAlphaBossak;
    }

    /// Destructor.
    ~ResidualBasedAdjointBossakScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        int num_threads = OpenMPUtils::GetNumThreads();
        mLeftHandSide.resize(num_threads);
        mResponseGradient.resize(num_threads);
        mFirstDerivsLHS.resize(num_threads);
        mFirstDerivsResponseGradient.resize(num_threads);
        mSecondDerivsLHS.resize(num_threads);
        mSecondDerivsResponseGradient.resize(num_threads);
        mAdjointValuesVector.resize(num_threads);
        mAdjointIndirectVector2.resize(num_threads);
        mAdjointIndirectVector3.resize(num_threads);
        mAuxAdjointIndirectVector1.resize(num_threads);

        InitializeNodeNeighbourCount(rModelPart.Nodes());

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(
                    ModelPart& rModelPart,
                    SystemMatrixType& rA,
                    SystemVectorType& rDx,
                    SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Get current time step.
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const ProcessInfo& r_last_process_info = r_current_process_info.GetPreviousSolutionStepInfo(1);

        // Note: solution is backwards in time, but we still want a positive time step
        // (it is the time step in the "forward" Bossak scheme).
        mTimeStep = r_last_process_info.GetValue(TIME) - r_current_process_info.GetValue(TIME);
        KRATOS_ERROR_IF(mTimeStep <= 0.0) << "Backwards in time solution is not decreasing time from last step." << std::endl;

        mInverseDt = 1.0 / mTimeStep;

        CalculateNodeNeighbourCount(rModelPart);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        this->UpdateAuxiliaryVariable(rModelPart);

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Update degrees of freedom: adjoint variables associated to the residual of the physical problem.
        this->mpDofUpdater->UpdateDofs(rDofSet,rDx);

        // Update adjoint variables associated to time integration.
        this->UpdateTimeSchemeAdjoints(rModelPart);

        KRATOS_CATCH("");
    }


    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(rModelPart.GetCommunicator().TotalProcesses() != 1) << "MPI version not implemented yet." << std::endl;

        KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << "No elements found in the ModelPart." << std::endl;

        // Check that the element and space dimensions match
        const unsigned int working_space_dimension = rModelPart.ElementsBegin()->WorkingSpaceDimension();
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size = static_cast<unsigned int>(r_process_info.GetValue(DOMAIN_SIZE));

        KRATOS_ERROR_IF(domain_size != 2 && domain_size != 3) <<
            "invalid DOMAIN_SIZE: " << domain_size << "." << std::endl;
        KRATOS_ERROR_IF(domain_size != working_space_dimension) <<
            "DOMAIN_SIZE " << domain_size << " not equal to the element's Working Space Dimension " <<
            working_space_dimension << "." << std::endl;

        return BaseType::Check(rModelPart); // Check elements and conditions.

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        // Check and resize rLHS and rRHS
        this->CheckAndResizeLocalSystem(pCurrentElement, rLHS_Contribution, rRHS_Contribution);

        // Contribution from variable gradients
        this->CalculateGradientContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from first derivative gradients
        this->CalculateFirstDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from second derivative gradients
        this->CalculateSecondDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contributions from the previos time step
        this->CalculatePreviousTimeStepContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Make the local contribution residual
        this->CalculateResidualLocalContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        LocalSystemVectorType RHS_Contribution;

        CalculateSystemContributions(
            pCurrentElement, rLHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                Condition::EquationIdVectorType& rEquationId,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        pCurrentCondition->CalculateLocalSystem(
            rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              Condition::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    void CheckAndResizeLocalSystem(Element::Pointer pCurrentElement,
                                   LocalSystemMatrixType& rLHS_Contribution,
                                   LocalSystemVectorType& rRHS_Contribution)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_residual_adjoint = mAdjointValuesVector[k];
        pCurrentElement->GetValuesVector(r_residual_adjoint);

        if (rRHS_Contribution.size() != r_residual_adjoint.size()) {
            rRHS_Contribution.resize(r_residual_adjoint.size(),false);
        }

        if (rLHS_Contribution.size1() != r_residual_adjoint.size()) {
            rLHS_Contribution.resize(r_residual_adjoint.size(),r_residual_adjoint.size(),false);
        }

        this->CheckAndResizeThreadStorage(r_residual_adjoint.size());
    }

    virtual void CalculateGradientContributions(Element::Pointer pCurrentElement,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        auto& r_lhs = mLeftHandSide[k];
        auto& r_response_gradient = mResponseGradient[k];
        pCurrentElement->CalculateLeftHandSide(r_lhs,rCurrentProcessInfo);
        r_response_function.CalculateGradient(*pCurrentElement, r_lhs, r_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) = r_lhs;
        noalias(rRHS_Contribution) = -1. * r_response_gradient;
    }

    virtual void CalculateFirstDerivativeContributions(Element::Pointer pCurrentElement,
                                                       LocalSystemMatrixType& rLHS_Contribution,
                                                       LocalSystemVectorType& rRHS_Contribution,
                                                       ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        const double first_deriv_coeff = mGammaNewmark / (mBetaNewmark*mTimeStep);
        auto& r_first_lhs = mFirstDerivsLHS[k];
        auto& r_first_response_gradient = mFirstDerivsResponseGradient[k];
        pCurrentElement->CalculateFirstDerivativesLHS(r_first_lhs,rCurrentProcessInfo);
        r_response_function.CalculateFirstDerivativesGradient(*pCurrentElement,r_first_lhs, r_first_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) += first_deriv_coeff * r_first_lhs;
        noalias(rRHS_Contribution) -= first_deriv_coeff * r_first_response_gradient;
    }

    virtual void CalculateSecondDerivativeContributions(Element::Pointer pCurrentElement,
                                                        LocalSystemMatrixType& rLHS_Contribution,
                                                        LocalSystemVectorType& rRHS_Contribution,
                                                        ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        const double second_deriv_coeff = mInverseDt * mInverseDt / mBetaNewmark;
        auto& r_second_lhs = mSecondDerivsLHS[k];
        auto& r_second_response_gradient = mSecondDerivsResponseGradient[k];
        pCurrentElement->CalculateSecondDerivativesLHS(r_second_lhs,rCurrentProcessInfo);
        r_second_lhs *= (1.0 - mAlphaBossak);
        r_response_function.CalculateSecondDerivativesGradient(*pCurrentElement,r_second_lhs, r_second_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) += second_deriv_coeff * r_second_lhs;
        noalias(rRHS_Contribution) -= second_deriv_coeff * r_second_response_gradient;
    }

    virtual void CalculatePreviousTimeStepContributions(Element::Pointer pCurrentElement,
                                                        LocalSystemMatrixType& rLHS_Contribution,
                                                        LocalSystemVectorType& rRHS_Contribution,
                                                        ProcessInfo& rCurrentProcessInfo)
    {
        const double old_adjoint_velocity_coeff = mInverseDt * (mBetaNewmark - mGammaNewmark * (mGammaNewmark + 0.5)) / (mBetaNewmark*mBetaNewmark);
        const double old_adjoint_acceleration_coeff = -1.0 * (mInverseDt*mInverseDt) * (mGammaNewmark + 0.5) / (mBetaNewmark*mBetaNewmark);
        const double second_deriv_coeff = mInverseDt * mInverseDt / mBetaNewmark;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int num_nodes = r_geometry.PointsNumber();
        const int k = OpenMPUtils::ThisThread();

        unsigned int local_index = 0;
        for (unsigned int i_node = 0; i_node < num_nodes; ++i_node) {
            auto& r_node = r_geometry[i_node];
            pCurrentElement->GetValue(GetFirstDerivativesIndirectVector)(
                i_node, mAdjointIndirectVector2[k], 1);
            pCurrentElement->GetValue(GetSecondDerivativesIndirectVector)(
                i_node, mAdjointIndirectVector3[k], 1);
            pCurrentElement->GetValue(GetAuxAdjointIndirectVector)(
                i_node, mAuxAdjointIndirectVector1[k], 1);
            const double weight = 1.0 / r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);

            for (unsigned int d = 0; d < mAdjointIndirectVector2[k].size(); d++) {
                rRHS_Contribution[local_index] += weight * second_deriv_coeff * mAuxAdjointIndirectVector1[k][d];
                rRHS_Contribution[local_index] += weight * old_adjoint_velocity_coeff * mAdjointIndirectVector2[k][d];
                rRHS_Contribution[local_index] += weight * old_adjoint_acceleration_coeff * mAdjointIndirectVector3[k][d];
                ++local_index;
            }
        }
    }

    virtual void CalculateResidualLocalContributions(Element::Pointer pCurrentElement,
                                                     LocalSystemMatrixType& rLHS_Contribution,
                                                     LocalSystemVectorType& rRHS_Contribution,
                                                     ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_residual_adjoint = mAdjointValuesVector[k];
        pCurrentElement->GetValuesVector(r_residual_adjoint);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, r_residual_adjoint);
    }

    virtual void InitializeNodeNeighbourCount(ModelPart::NodesContainerType& rNodes)
    {
        // This loop should not be omp parallel
        // The operation is not threadsafe if the value is uninitialized
        for (auto& r_node : rNodes)
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0.0);
    }

    virtual void CalculateNodeNeighbourCount(ModelPart& rModelPart)
    {
        // Calculate number of neighbour elements for each node.
        const int num_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < num_nodes; i++) {
            Node<3>& r_node = *(rModelPart.Nodes().begin()+i);
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0.0);
        }

        const int num_elements = rModelPart.NumberOfElements();
        #pragma omp parallel for
        for (int i = 0; i < num_elements; i++) {
            Element& r_element = *(rModelPart.Elements().begin()+i);
            Geometry<Node<3>>& r_geometry = r_element.GetGeometry();
            for (unsigned int j = 0; j < r_geometry.PointsNumber(); j++) {
                double& r_num_neighbour = r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                #pragma omp atomic
                r_num_neighbour += 1.0;
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);
    }

    virtual void UpdateTimeSchemeAdjoints(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        auto lambda2_vars = Internals::AdjointBossakSchemeUtils().Execute(
            GetFirstDerivativesVariables, rModelPart.Elements());
        auto lambda3_vars = Internals::AdjointBossakSchemeUtils().Execute(
            GetSecondDerivativesVariables, rModelPart.Elements());
        Internals::AdjointBossakSchemeUtils().SetToZero_AdjointVars(lambda2_vars, rModelPart.Nodes());
        Internals::AdjointBossakSchemeUtils().SetToZero_AdjointVars(lambda3_vars, rModelPart.Nodes());
        auto& r_response_function = *(this->mpResponseFunction);

        const double a22 = 1.0 - mGammaNewmark/mBetaNewmark;
        const double a23 = -1.0 / (mBetaNewmark*mTimeStep);
        const double a32 = (1.0 - 0.5*mGammaNewmark/mBetaNewmark)*mTimeStep;
        const double a33 = (1.0 - 0.5/mBetaNewmark);

        const int number_of_elements = rModelPart.NumberOfElements();
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        Vector aux_adjoint2, aux_adjoint3;
        std::vector<IndirectScalar<double>> adjoint2, adjoint3, aux_adjoint,
            adjoint2_old, adjoint3_old;
        #pragma omp parallel for private(aux_adjoint2, aux_adjoint3, adjoint2, adjoint3, aux_adjoint, adjoint2_old, adjoint3_old)
        for (int i = 0; i < number_of_elements; i++) {
            Element& r_element = *(rModelPart.ElementsBegin()+i);
            const int k = OpenMPUtils::ThisThread();

            r_element.GetValuesVector(mAdjointValuesVector[k]);
            this->CheckAndResizeThreadStorage(mAdjointValuesVector[k].size());

            r_element.CalculateFirstDerivativesLHS(mFirstDerivsLHS[k],r_process_info);
            r_response_function.CalculateFirstDerivativesGradient(r_element,mFirstDerivsLHS[k],mFirstDerivsResponseGradient[k],r_process_info);

            r_element.CalculateSecondDerivativesLHS(mSecondDerivsLHS[k],r_process_info);
            mSecondDerivsLHS[k] *= (1.0 - mAlphaBossak);
            r_response_function.CalculateSecondDerivativesGradient(r_element,mSecondDerivsLHS[k],mSecondDerivsResponseGradient[k],r_process_info);

            if (aux_adjoint2.size() != mFirstDerivsResponseGradient[k].size())
                aux_adjoint2.resize(mFirstDerivsResponseGradient[k].size(), false);
            noalias(aux_adjoint2) = -mFirstDerivsResponseGradient[k] - prod(mFirstDerivsLHS[k], mAdjointValuesVector[k]);
            if (aux_adjoint3.size() != mSecondDerivsResponseGradient[k].size())
                aux_adjoint3.resize(mSecondDerivsResponseGradient[k].size(), false);
            noalias(aux_adjoint3) = - mSecondDerivsResponseGradient[k] - prod(mSecondDerivsLHS[k],mAdjointValuesVector[k]);

            // Assemble the contributions to the corresponding nodal unknowns.
            unsigned int local_index = 0;
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
                r_element.GetValue(GetFirstDerivativesIndirectVector)(i_node, adjoint2, 0);
                r_element.GetValue(GetSecondDerivativesIndirectVector)(i_node, adjoint3, 0);
                r_element.GetValue(GetFirstDerivativesIndirectVector)(i_node, adjoint2_old, 1);
                r_element.GetValue(GetSecondDerivativesIndirectVector)(i_node, adjoint3_old, 1);
                r_element.GetValue(GetAuxAdjointIndirectVector)(i_node, aux_adjoint, 1);
                Node<3>& r_node = r_geometry[i_node];
                const double weight = 1.0 / r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                r_node.SetLock();
                for (unsigned int d = 0; d < adjoint2.size(); ++d) {

                    adjoint2[d] += aux_adjoint2[local_index];
                    adjoint2[d] += a22 * weight * adjoint2_old[d];
                    adjoint2[d] += a23 * weight * adjoint3_old[d];
                    adjoint3[d] += aux_adjoint3[local_index];
                    adjoint3[d] += a32 * weight * adjoint2_old[d];
                    adjoint3[d] += a33 * weight * adjoint3_old[d];
                    adjoint3[d] += weight * aux_adjoint[d];
                    ++local_index;
                }
                r_node.UnSetLock();
            }
        }

        // Finalize global assembly
        Internals::AdjointBossakSchemeUtils().Assemble_AdjointVars(lambda2_vars, rModelPart.GetCommunicator());
        Internals::AdjointBossakSchemeUtils().Assemble_AdjointVars(lambda3_vars, rModelPart.GetCommunicator());
        KRATOS_CATCH("");
    }

    virtual void UpdateAuxiliaryVariable(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        auto aux_vars = Internals::AdjointBossakSchemeUtils().Execute(
            GetAuxAdjointVariables, rModelPart.Elements());
        Internals::AdjointBossakSchemeUtils().SetToZero_AdjointVars(
            aux_vars, rModelPart.Nodes());
        auto& r_response_function = *(this->mpResponseFunction);

        // Loop over elements to assemble the remaining terms
        const int number_of_elements = rModelPart.NumberOfElements();
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        std::vector<IndirectScalar<double>> aux_adjoint_on_node;
        Vector aux_adjoint_vector;
#pragma omp parallel for private(aux_adjoint_on_node, aux_adjoint_vector)
        for (int i = 0; i < number_of_elements; i++) {
            Element& r_element = *(rModelPart.ElementsBegin()+i);
            const int k = OpenMPUtils::ThisThread();

            r_element.GetValuesVector(mAdjointValuesVector[k]);
            this->CheckAndResizeThreadStorage(mAdjointValuesVector[k].size());

            r_element.CalculateSecondDerivativesLHS(mSecondDerivsLHS[k],r_process_info);
            mSecondDerivsLHS[k] *= mAlphaBossak;
            r_response_function.CalculateSecondDerivativesGradient(r_element,mSecondDerivsLHS[k],mSecondDerivsResponseGradient[k],r_process_info);

            if (aux_adjoint_vector.size() != mSecondDerivsLHS[k].size1())
                aux_adjoint_vector.resize(mSecondDerivsLHS[k].size1(), false);
            noalias(aux_adjoint_vector) = -prod(mSecondDerivsLHS[k], mAdjointValuesVector[k]) - mSecondDerivsResponseGradient[k];

            // Assemble the contributions to the corresponding nodal unknowns.
            unsigned int local_index = 0;
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {

                Node<3>& r_node = r_geometry[i_node];
                r_element.GetValue(GetAuxAdjointIndirectVector)(i_node, aux_adjoint_on_node, 0);

                r_node.SetLock();
                for (unsigned int d = 0; d < aux_adjoint_on_node.size(); ++d) {

                    aux_adjoint_on_node[d] += aux_adjoint_vector[local_index];
                    ++local_index;
                }
                r_node.UnSetLock();
            }
        }

        // Finalize global assembly
        Internals::AdjointBossakSchemeUtils().Assemble_AdjointVars(aux_vars, rModelPart.GetCommunicator());
        KRATOS_CATCH("");
    }

    /// Free memory allocated by this class.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mAlphaBossak;
    double mBetaNewmark;
    double mGammaNewmark;
    double mTimeStep;
    double mInverseDt;

    std::vector< LocalSystemMatrixType > mLeftHandSide;
    std::vector< LocalSystemVectorType > mResponseGradient;
    std::vector< LocalSystemMatrixType > mFirstDerivsLHS;
    std::vector< LocalSystemVectorType > mFirstDerivsResponseGradient;
    std::vector< LocalSystemMatrixType > mSecondDerivsLHS;
    std::vector< LocalSystemVectorType > mSecondDerivsResponseGradient;
    std::vector< LocalSystemVectorType > mAdjointValuesVector;
    std::vector<std::vector<IndirectScalar<double>>> mAdjointIndirectVector2;
    std::vector<std::vector<IndirectScalar<double>>> mAdjointIndirectVector3;
    std::vector<std::vector<IndirectScalar<double>>> mAuxAdjointIndirectVector1;

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckAndResizeThreadStorage(unsigned int SystemSize) {
        const int k = OpenMPUtils::ThisThread();

        if (mLeftHandSide[k].size1() != SystemSize || mLeftHandSide[k].size2() != SystemSize) {
            mLeftHandSide[k].resize(SystemSize, SystemSize, false);
        }

        if (mFirstDerivsLHS[k].size1() != SystemSize || mFirstDerivsLHS[k].size2() != SystemSize) {
            mFirstDerivsLHS[k].resize(SystemSize, SystemSize, false);
        }

        if (mSecondDerivsLHS[k].size1() != SystemSize || mSecondDerivsLHS[k].size2() != SystemSize) {
            mSecondDerivsLHS[k].resize(SystemSize, SystemSize, false);
        }

        if (mResponseGradient[k].size() != SystemSize) {
            mResponseGradient[k].resize(SystemSize,false);
        }

        if (mFirstDerivsResponseGradient[k].size() != SystemSize) {
            mFirstDerivsResponseGradient[k].resize(SystemSize,false);
        }

        if (mSecondDerivsResponseGradient[k].size() != SystemSize) {
            mSecondDerivsResponseGradient[k].resize(SystemSize,false);
        }
    }

    Variable< array_1d<double,3> > GetVariableFromParameters(
        Parameters& rParameters,
        const std::string& rLabel) {
        KRATOS_TRY;

        std::string variable_name = rParameters[rLabel].GetString();
        bool have_variable = KratosComponents< Variable< array_1d<double,3> > >::Has(variable_name);

        KRATOS_ERROR_IF_NOT( have_variable ) << "Variable " << variable_name <<
        " passed as Parameters argument \"" << rLabel << "\" is not an array_1d<double,3> Variable defined in Kratos." << std::endl <<
        "If it is an application variable, you may need to imoirt the application that defines it." << std::endl;

        return KratosComponents< Variable< array_1d<double,3> > >::Get(variable_name);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedAdjointBossakScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED defined */
