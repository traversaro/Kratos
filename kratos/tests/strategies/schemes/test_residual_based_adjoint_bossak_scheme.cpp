#include<iostream>
#include<iomanip> // for debug only
#include<algorithm>
#include<exception>
#include <sstream>

#include "testing/testing.h"
#include "includes/define.h"
#include "includes/shared_pointers.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "containers/pointer_vector.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "utilities/indirect_scalar.h"
#include "response_functions/adjoint_response_function.h"
#include "response_functions/sensitivity_builder.h"


namespace Kratos
{
namespace Testing
{
struct ResultsData
{
    virtual void StoreCurrentResult(const ModelPart& rModelPart) = 0;
    virtual void LoadCurrentResult(ModelPart& rModelPart) = 0;
};

const double eps = 0./*1.e-4*/;
namespace NonLinearMassSpringDamper
{
/**
 * @class TwoMassSpringDamperSystem
 * @brief A system of two mass-spring-dampers for testing a second-order ode.
 * @details Taken from L.F. Fernandez, D.A. Tortorelli, Semi-analytical
 * sensitivity analysis for nonlinear transient problems.
 * 
 *  |                _____                 _____
 *  |---[ Damper ]--|  m1 |---[ Damper ]--|  m2 |
 *  |-----/\/\/\----|_____|-----/\/\/\----|_____|
 *  |
 * 
 * Spring force: fe = x + kd * x^3
 * Damper force: fc = kc * x'
 * 
 */
class TwoMassSpringDamperSystem : public Element
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TwoMassSpringDamperSystem);

    static Pointer Create(Node<3>::Pointer pNode1, Node<3>::Pointer pNode2)
    {
        auto nodes = PointerVector<Node<3>>{};
        nodes.push_back(pNode1);
        nodes.push_back(pNode2);
        return Kratos::make_shared<TwoMassSpringDamperSystem>(nodes);
    }

    TwoMassSpringDamperSystem(const NodesArrayType& ThisNodes)
        : Element(0, ThisNodes)
    {
    }

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        rResult.resize(2);
        rResult[0] = this->GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
        rResult[1] = this->GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
    }

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        rElementalDofList.resize(2);
        rElementalDofList[0] = this->GetGeometry()[0].pGetDof(DISPLACEMENT_X);
        rElementalDofList[1] = this->GetGeometry()[1].pGetDof(DISPLACEMENT_X);
    }

    void GetValuesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    }

    void GetFirstDerivativesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
    }

    void GetSecondDerivativesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION_X);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix.resize(2, 2, false);
        const double& x1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        const double& x2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
        rLeftHandSideMatrix(0, 0) = 2. + 3. * kd * x1 * x1 + 3. * kd * (x2 - x1) * (x2 - x1);
        rLeftHandSideMatrix(0, 1) = -1. - 3. * kd * (x2 - x1) * (x2 - x1);
        rLeftHandSideMatrix(1, 0) = -1. - 3. * kd * (x2 - x1) * (x2 - x1);
        rLeftHandSideMatrix(1, 1) = 1. + 3. * kd * (x2 - x1) * (x2 - x1);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector.resize(2, false);
        const double& x1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        const double& x2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
        const double x21 = x2 - x1;
        rRightHandSideVector(0) = -(2. * x1 - x2 + kd * x1 * x1 * x1 - kd * x21 * x21 * x21);
        rRightHandSideVector(1) = -(-x1 + x2 + kd * x21 * x21 * x21);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        rMassMatrix.resize(2, 2, false);
        rMassMatrix(0, 0) = m1;
        rMassMatrix(0, 1) = 0.;
        rMassMatrix(1, 1) = m2;
        rMassMatrix(1, 0) = 0.;
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        rDampingMatrix.resize(2, 2, false);
        rDampingMatrix(0, 0) = 2. * kc;
        rDampingMatrix(0, 1) = -kc;
        rDampingMatrix(1, 1) = kc;
        rDampingMatrix(1, 0) = -kc;
    }

private:
    const double m1 = 1.;
    const double m2 = 1.;
    const double kd = 1. + eps;
    const double kc = 0.1;
};

class TwoMassSpringDamperSystemAdjoint : public Element
{
    struct GetFirstDerivativesVectorImpl
    {
        Element* mpElement;
        void operator()(std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
        }
    };

    struct GetSecondDerivativesVectorImpl
    {
        Element* mpElement;
        void operator()(std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
        }
    };

    struct GetAuxAdjointVectorImpl
    {
        Element* mpElement;
        void operator()(std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
        }
    };

    struct GetFirstDerivativesVariablesImpl
    {
        Element* mpElement;
        void operator()(std::vector<VariableData const*>& rVariables)
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
        }
    };
    
    struct GetSecondDerivativesVariablesImpl
    {
        Element* mpElement;
        void operator()(std::vector<VariableData const*>& rVariables)
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
        }
    };
    
    struct GetAuxAdjointVariablesImpl
    {
        Element* mpElement;
        void operator()(std::vector<VariableData const*>& rVariables)
        {
            rVariables.resize(1);
            rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
        }
    };

public:
    KRATOS_CLASS_POINTER_DEFINITION(TwoMassSpringDamperSystemAdjoint);

    static Pointer Create(Node<3>::Pointer pNode1, Node<3>::Pointer pNode2)
    {
        auto nodes = PointerVector<Node<3>>{};
        nodes.push_back(pNode1);
        nodes.push_back(pNode2);
        return Kratos::make_shared<TwoMassSpringDamperSystemAdjoint>(nodes);
    }

    TwoMassSpringDamperSystemAdjoint(const NodesArrayType& ThisNodes)
        : Element(0, ThisNodes), mPrimalElement(ThisNodes)
    {
        SetValue(GetFirstDerivativesIndirectVector, GetFirstDerivativesVectorImpl{this});
        SetValue(GetSecondDerivativesIndirectVector, GetSecondDerivativesVectorImpl{this});
        SetValue(GetAuxAdjointIndirectVector, GetAuxAdjointVectorImpl{this});
        SetValue(GetFirstDerivativesVariables, GetFirstDerivativesVariablesImpl{this});
        SetValue(GetSecondDerivativesVariables, GetSecondDerivativesVariablesImpl{this});
        SetValue(GetAuxAdjointVariables, GetAuxAdjointVariablesImpl{this});
    }

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        rResult.resize(2);
        rResult[0] = this->GetGeometry()[0].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
        rResult[1] = this->GetGeometry()[1].GetDof(ADJOINT_FLUID_VECTOR_1_X).EquationId();
    }

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        rElementalDofList.resize(2);
        rElementalDofList[0] = this->GetGeometry()[0].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
        rElementalDofList[1] = this->GetGeometry()[1].pGetDof(ADJOINT_FLUID_VECTOR_1_X);
    }

    void GetValuesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X);
    }

    void GetFirstDerivativesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2_X);
    }

    void GetSecondDerivativesVector(Vector& values, int Step = 0) override
    {
        values.resize(2);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_X);
        values[1] = this->GetGeometry()[1].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_X);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalElement.CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
        KRATOS_CATCH("");
    }

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        if (rDesignVariable == NORMAL_SENSITIVITY)
        {
            rOutput.resize(1, 2, false);
            const double& x1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
            const double& x2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
            const double x21 = x2 - x1;
            rOutput(0, 0) = -(x1 * x1 * x1 - x21 * x21 * x21);
            rOutput(0, 1) = -(x21 * x21 * x21);
        }
        else
        {
            KRATOS_ERROR << "Invalid variable: " << rDesignVariable << std::endl;
        }
        KRATOS_CATCH("");
    }

    ///@}

private:
    TwoMassSpringDamperSystem mPrimalElement;
};

class TwoMassSpringDamperSystemResponseFunction : public AdjointResponseFunction
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(TwoMassSpringDamperSystemResponseFunction);
        
        TwoMassSpringDamperSystemResponseFunction(ModelPart& rModelPart)
            : mrModelPart(rModelPart)
        {
        }

        void CalculateGradient(const Element& rAdjointElement,
                               const Matrix& rResidualGradient,
                               Vector& rResponseGradient,
                               const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(2, false);
            rResponseGradient(0) = 2. * rAdjointElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
            rResponseGradient(1) = 2. * rAdjointElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
        }

        void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                               const Matrix& rResidualGradient,
                                               Vector& rResponseGradient,
                                               const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(2, false);
            rResponseGradient(0) = 2. * rAdjointElement.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
            rResponseGradient(1) = 2. * rAdjointElement.GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
        }

        void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                                const Matrix& rResidualGradient,
                                                Vector& rResponseGradient,
                                                const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(2, false);
            rResponseGradient(0) = 0.;
            rResponseGradient(1) = 0.;
        }

        void CalculatePartialSensitivity(Element& rAdjointElement,
                                         const Variable<double>& rVariable,
                                         const Matrix& rSensitivityMatrix,
                                         Vector& rSensitivityGradient,
                                         const ProcessInfo& rProcessInfo) override
        {
            rSensitivityGradient.resize(1, false);
            rSensitivityGradient(0) = 0.;
        }

        void CalculatePartialSensitivity(Element& rAdjointElement,
                                         const Variable<array_1d<double, 3>>& rVariable,
                                         const Matrix& rSensitivityMatrix,
                                         Vector& rSensitivityGradient,
                                         const ProcessInfo& rProcessInfo) override
        {
            rSensitivityGradient.resize(1, false);
            rSensitivityGradient(0) = 0.;
        }

        double CalculateValue() override
        {
            const double& x1 =
                mrModelPart.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT_X);
            const double& x2 =
                mrModelPart.GetNode(2).FastGetSolutionStepValue(DISPLACEMENT_X);
            const double& v1 = mrModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY_X);
            const double& v2 = mrModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY_X);
            return x1 * x1 + x2 * x2 + v1 * v1 + v2 * v2;
        }

    private:
        ModelPart& mrModelPart;
};

struct TwoMassSpringDamperSystemResultsData : ResultsData
{
    std::vector<double> time;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> v1;
    std::vector<double> v2;
    std::vector<double> a1;
    std::vector<double> a2;

    void StoreCurrentResult(const ModelPart& rModelPart) override
    {
        this->time.push_back(rModelPart.GetProcessInfo()[TIME]);
        auto& node1 = rModelPart.GetNode(1);
        auto& node2 = rModelPart.GetNode(2);
        this->x1.push_back(node1.FastGetSolutionStepValue(DISPLACEMENT_X));
        this->v1.push_back(node1.FastGetSolutionStepValue(VELOCITY_X));
        this->a1.push_back(node1.FastGetSolutionStepValue(ACCELERATION_X));
        this->x2.push_back(node2.FastGetSolutionStepValue(DISPLACEMENT_X));
        this->v2.push_back(node2.FastGetSolutionStepValue(VELOCITY_X));
        this->a2.push_back(node2.FastGetSolutionStepValue(ACCELERATION_X));
    }

    void LoadCurrentResult(ModelPart& rModelPart) override
    {
        const double current_time = rModelPart.GetProcessInfo()[TIME];
        if (current_time < 1e-8 || current_time > this->time.back() + 1e-8)
        {
            std::stringstream ss;
            ss << "Adjoint time = " << current_time
               << " outside of primal solution time range!\n";
            throw std::runtime_error{ss.str()};
        }
        auto it = std::find_if(
            this->time.begin(), this->time.end(),
            [current_time](const double& t) -> bool { return current_time <= t; });
        std::size_t pos = it - this->time.begin();
        if (pos == this->time.size())
        {
            --pos;
        }
        else if (pos > 0)
        {
            const auto t0 = this->time.at(pos - 1);
            const auto t1 = this->time.at(pos);
            if (std::abs(current_time - t0) < 1e-8)
            {
                pos = pos - 1;
            }
            else if (!(std::abs(current_time - t1) < 1e-8))
            {
                std::stringstream ss;
                ss << "Adjoint time = " << current_time
                   << " does not match primal solution time!\n";
                throw std::runtime_error{ss.str()};
            }
        }
        auto& node1 = rModelPart.GetNode(1);
        auto& node2 = rModelPart.GetNode(2);
        node1.FastGetSolutionStepValue(DISPLACEMENT_X) = this->x1.at(pos);
        node1.FastGetSolutionStepValue(VELOCITY_X) = this->v1.at(pos);
        node1.FastGetSolutionStepValue(ACCELERATION_X) = this->a1.at(pos);
        node2.FastGetSolutionStepValue(DISPLACEMENT_X) = this->x2.at(pos);
        node2.FastGetSolutionStepValue(VELOCITY_X) = this->v2.at(pos);
        node2.FastGetSolutionStepValue(ACCELERATION_X) = this->a2.at(pos);
    }
};

// double ResponseValue(TwoMassSpringDamperSystemResultsData& rd)
// {
//     double response_value = 0.;
//     for (std::size_t i = 0; i < rd.time.size(); ++i)
//     {
//         const double x1 = rd.x1[i];
//         const double x2 = rd.x2[i];
//         const double v1 = rd.v1[i];
//         const double v2 = rd.v2[i];
//         const double delta_time = rd.time[i] - rd.time[i - 1];
//         response_value += delta_time * (x1 * x1 + x2 * x2 + v1 * v1 + v2 * v2);
//     }
//     return response_value;
// }
}

namespace Solvers
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

    SolvingStrategyType::Pointer CreateResidualBasedNewtonRaphsonStrategy(ModelPart& rModelPart)
    {
        LinearSolverType::Pointer p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        SchemeType::Pointer p_scheme =
            Kratos::make_shared<ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>>(/*-0.3*/ 0.);
        ConvergenceCriteriaType::Pointer p_conv_criteria =
            Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(
                1e-26, 1e-27);
        return Kratos::make_shared<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            rModelPart, p_scheme, p_linear_solver, p_conv_criteria, 40, true, false, true);
    }

    class AdjointBossakSolver
    {
    public:
        AdjointBossakSolver(ModelPart& rModelPart,
                            Kratos::shared_ptr<ResultsData> pResultsData,
                            AdjointResponseFunction::Pointer pResponseFunction)
            : mrModelPart(rModelPart), mpResultsData(pResultsData)
        {
            LinearSolverType::Pointer p_linear_solver =
                Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
            auto scheme_settings = Parameters{R"({ "alpha_bossak": 0.0 })"};
            SchemeType::Pointer p_adjoint_scheme =
                Kratos::make_shared<ResidualBasedAdjointBossakScheme<SparseSpaceType, LocalSpaceType>>(
                    scheme_settings, pResponseFunction);
            mpSolver =
                Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
                    rModelPart, p_adjoint_scheme, p_linear_solver);
        }

        double Solve()
        {
            mpResultsData->LoadCurrentResult(mrModelPart);
            return mpSolver->Solve();
        }

    private:
        ModelPart& mrModelPart;
        Kratos::shared_ptr<ResultsData> mpResultsData;
        SolvingStrategyType::Pointer mpSolver;
    };
    }

    std::ostream& operator<<(std::ostream& os, const Matrix& m)
    {
        os << "\nPrinting Matrix:\n";
        os << std::scientific << std::fixed;
        for (std::size_t i = 0; i < m.size1(); ++i)
        {
            for (std::size_t j = 0; j < m.size2(); ++j)
                os << m(i, j) << ' ';
            os << std::endl;
        }
        return os;
}

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    os << "\nPrinting Vector:\n";
    os << std::scientific << std::fixed;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        os << v(i) << std::endl;
    }
    return os;
}

KRATOS_TEST_CASE_IN_SUITE(ResidualBasedAdjointBossak_TwoMassSpringDamperSystem, KratosCoreSchemesFastSuite)
{
    using namespace NonLinearMassSpringDamper;
    using namespace Solvers;
    ModelPart model_part("test");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
    model_part.SetBufferSize(2);
    for (auto& r_node : model_part.Nodes())
    {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
    }
    auto& node1 = model_part.GetNode(1);
    auto& node2 = model_part.GetNode(2);
    auto p_element = TwoMassSpringDamperSystem::Create(model_part.pGetNode(1),
                                                       model_part.pGetNode(2));
    model_part.AddElement(p_element);
    auto p_solver = CreateResidualBasedNewtonRaphsonStrategy(model_part);
    node2.FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    node1.FastGetSolutionStepValue(ACCELERATION_X) = 2.0;
    node2.FastGetSolutionStepValue(ACCELERATION_X) =-2.0;
    Kratos::shared_ptr<TwoMassSpringDamperSystemResultsData> p_results_data =
        Kratos::make_shared<TwoMassSpringDamperSystemResultsData>();
    const double end_time = 0.1;
    const double start_time = 0.;
    const std::size_t N = 5;
    const double delta_time = (end_time - start_time) / N;
    model_part.CloneTimeStep(start_time - delta_time);
    model_part.CloneTimeStep(start_time);
    for (double current_time = start_time; current_time < end_time;)
    {
        current_time += delta_time;
        model_part.CloneTimeStep(current_time);
        p_solver->Solve();
        p_results_data->StoreCurrentResult(model_part);
    }
    ModelPart adjoint_model_part("test");
    adjoint_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    adjoint_model_part.AddNodalSolutionStepVariable(REACTION);
    adjoint_model_part.AddNodalSolutionStepVariable(VELOCITY);
    adjoint_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    adjoint_model_part.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
    adjoint_model_part.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
    adjoint_model_part.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
    adjoint_model_part.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);
    adjoint_model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY);
    adjoint_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    adjoint_model_part.CreateNewNode(2, 0.0, 0.0, 0.0);
    adjoint_model_part.SetBufferSize(2);
    for (auto& r_node : adjoint_model_part.Nodes())
    {
        r_node.AddDof(ADJOINT_FLUID_VECTOR_1_X, REACTION_X);
    }
    auto p_adjoint_element = TwoMassSpringDamperSystemAdjoint::Create(
        adjoint_model_part.pGetNode(1), adjoint_model_part.pGetNode(2));
    adjoint_model_part.AddElement(p_adjoint_element);
    auto p_response_function =
        Kratos::make_shared<TwoMassSpringDamperSystemResponseFunction>(adjoint_model_part);
    AdjointBossakSolver adjoint_solver(adjoint_model_part, p_results_data, p_response_function);
    SensitivityBuilder sensitivity_builder(
        Parameters{R"({ "element_sensitivity_variables": ["NORMAL_SENSITIVITY"] })"},
        adjoint_model_part, p_response_function);
    sensitivity_builder.Initialize();
    adjoint_model_part.CloneTimeStep(end_time + 2. * delta_time);
    adjoint_model_part.CloneTimeStep(end_time + delta_time);
    for (double current_time = end_time + delta_time; current_time >= start_time + 1.5 * delta_time;)
    {
        current_time -= delta_time;
        adjoint_model_part.CloneTimeStep(current_time);
        adjoint_solver.Solve();
        sensitivity_builder.UpdateSensitivities();
    }
    const double adjoint_sensitivity = p_adjoint_element->GetValue(NORMAL_SENSITIVITY);
    const double fd_sensitivity = (0.1025147391 - 0.1025144844) / 1.e-4 /* = 0.00254700000007490*/;
    KRATOS_CHECK_NEAR(adjoint_model_part.GetNode(1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X), 0.0220776047, 1e-8);
    KRATOS_CHECK_NEAR(adjoint_model_part.GetNode(2).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X),-0.0141543232, 1e-8);
    KRATOS_CHECK_NEAR(adjoint_sensitivity, fd_sensitivity, 1e-7);
}
}
}