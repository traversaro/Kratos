// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    
//

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/shared_pointers.h"
#include "includes/kratos_components.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "utilities/sensitivity_builder.h"
#include "containers/model.h"

// Application includes
#include "custom_elements/total_lagrangian.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_response_functions/node_displacement_response_function.h"

namespace Kratos
{
namespace Testing
{
namespace
{

class PrimalModelPartFactory
{
    public:
    PrimalModelPartFactory(ModelPart& rModelPart) : mrModelPart(rModelPart) {}
    void Execute()
    {
        AddVariables();
        CreateNodes();
        mrModelPart.GetProcessInfo()[DOMAIN_SIZE] =
            KratosComponents<Element>::Get(mElementType).WorkingSpaceDimension();
        CreateConstitutiveLaw();
        CreateElements();
        mrModelPart.SetBufferSize(2);
        AddDofs();
        AssignBCs();
    }

private:
    void AddVariables()
    {
        mrModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        mrModelPart.AddNodalSolutionStepVariable(REACTION);
        mrModelPart.AddNodalSolutionStepVariable(VELOCITY);
        mrModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        mrModelPart.AddNodalSolutionStepVariable(DENSITY);
        mrModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        mrModelPart.AddNodalSolutionStepVariable(THICKNESS);
    }

    void CreateNodes()
    {
        mrModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        mrModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        mrModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    }

    void CreateConstitutiveLaw()
    {
        auto p_prop = mrModelPart.pGetProperties(mPropertiesId);
        (*p_prop)[CONSTITUTIVE_LAW] =
            LinearPlaneStrain::Pointer(new LinearPlaneStrain());
        (*p_prop)[DENSITY] = 1000.0;
        (*p_prop)[YOUNG_MODULUS] = 1400000.0;
        (*p_prop)[POISSON_RATIO] = 0.2;
        (*p_prop)[RAYLEIGH_ALPHA] = 0.; // 0.02;
        (*p_prop)[RAYLEIGH_BETA] = 0.;  // 0.03;
    }

    void CreateElements()
    {
        const auto& r_elem = KratosComponents<Element>::Get(mElementType);
        std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().PointsNumber());
        KRATOS_ERROR_IF(node_ids.size() != 3)
            << "Inconsistent number of nodes.\n";
        auto p_prop = mrModelPart.pGetProperties(mPropertiesId);
        (*p_prop)[VOLUME_ACCELERATION] = ZeroVector(3);
        (*p_prop)[VOLUME_ACCELERATION](1) = 100.0;
        // Element 1.
        node_ids.at(0) = 1;
        node_ids.at(1) = 2;
        node_ids.at(2) = 3;
        mrModelPart.CreateNewElement("TotalLagrangianElement2D3N", 1, node_ids, p_prop);
    }

    void AddDofs()
    {
        for (auto& r_node : mrModelPart.Nodes())
        {
            r_node.AddDof(DISPLACEMENT_X, REACTION_X);
            r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        }
    }

    void AssignBCs()
    {
        mrModelPart.GetNode(1).Fix(DISPLACEMENT_X);
        mrModelPart.GetNode(1).Fix(DISPLACEMENT_Y);
        mrModelPart.GetNode(3).Fix(DISPLACEMENT_X);
        mrModelPart.GetNode(3).Fix(DISPLACEMENT_Y);
    }

    ModelPart& mrModelPart;
    const std::string mElementType = "TotalLagrangianElement2D3N";
    std::size_t mPropertiesId = 1;
};

class AdjointModelPartFactory
{
    public:
        AdjointModelPartFactory(ModelPart const& rPrimalModelPart, ModelPart& rAdjointModelPart)
            : mrPrimalModelPart(rPrimalModelPart), mrAdjointModelPart(rAdjointModelPart)
        {
        }

        void Execute()
        {
            AddVariables();
            CreateNodes();
            mrAdjointModelPart.GetProcessInfo() = mrPrimalModelPart.GetProcessInfo();
            SetProperties();
            CreateElements();
            mrAdjointModelPart.SetBufferSize(mrPrimalModelPart.GetBufferSize());
            AddDofs();
            AssignBCs();
            for (auto& r_node : mrAdjointModelPart.Nodes())
                r_node.SetValue(UPDATE_SENSITIVITIES, true);
            SetPrimalSolutionStepData();
        }

    private:
        void AddVariables()
        {
            mrAdjointModelPart.GetNodalSolutionStepVariablesList() =
                mrPrimalModelPart.GetNodalSolutionStepVariablesList();
            mrAdjointModelPart.AddNodalSolutionStepVariable(ADJOINT_DISPLACEMENT);
            mrAdjointModelPart.AddNodalSolutionStepVariable(ADJOINT_VECTOR_2);
            mrAdjointModelPart.AddNodalSolutionStepVariable(ADJOINT_VECTOR_3);
            mrAdjointModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_VECTOR_1);
            mrAdjointModelPart.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY);
        }

        void CreateNodes()
        {
            for (const auto& r_node : mrPrimalModelPart.Nodes())
            {
                auto p_node = mrAdjointModelPart.CreateNewNode(
                    r_node.Id(), r_node.X0(), r_node.Y0(), r_node.Z0());
                p_node->Coordinates() = r_node.Coordinates();
            }
        }

        void SetProperties()
        {
            auto p_primal_prop = mrPrimalModelPart.ElementsBegin()->pGetProperties();
            auto p_adjoint_prop = mrAdjointModelPart.pGetProperties(p_primal_prop->Id());
            *p_adjoint_prop = *p_primal_prop;
        }

        void CreateElements()
        {
            std::vector<ModelPart::IndexType> node_ids;
            for (const auto& r_elem : mrPrimalModelPart.Elements())
            {
                auto p_prop = mrAdjointModelPart.pGetProperties(r_elem.pGetProperties()->Id());
                const auto& r_points = r_elem.GetGeometry().Points();
                node_ids.resize(r_points.size());
                for (std::size_t i = 0; i < r_points.size(); ++i)
                    node_ids[i] = r_points[i].Id();
                mrAdjointModelPart.CreateNewElement(mElementType, r_elem.Id(),
                                                    node_ids, p_prop);
            }
        }

        void AddDofs()
        {
            for (auto& r_node : mrAdjointModelPart.Nodes())
            {
                r_node.AddDof(ADJOINT_DISPLACEMENT_X);
                r_node.AddDof(ADJOINT_DISPLACEMENT_Y);
                r_node.AddDof(ADJOINT_DISPLACEMENT_Z);
            }
        }

        void AssignBCs()
        {
            for (auto& r_adjoint_node : mrAdjointModelPart.Nodes())
            {
                const auto& r_primal_node =
                    mrPrimalModelPart.GetNode(r_adjoint_node.Id());
                if (r_primal_node.IsFixed(DISPLACEMENT_X))
                    r_adjoint_node.Fix(ADJOINT_DISPLACEMENT_X);
                if (r_primal_node.IsFixed(DISPLACEMENT_Y))
                    r_adjoint_node.Fix(ADJOINT_DISPLACEMENT_Y);
            }
        }

        void SetPrimalSolutionStepData()
        {
            KRATOS_TRY;
            for (const auto& r_primal_node : mrPrimalModelPart.Nodes())
            {
                const auto p_variables_list = r_primal_node.pGetVariablesList();
                auto& r_adjoint_node =
                    mrAdjointModelPart.GetNode(r_primal_node.Id());
                const std::size_t queue_size = r_primal_node.SolutionStepData().QueueSize();

                for (const auto& r_variable_data : *p_variables_list)
                {
                    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(
                            r_variable_data.Name()))
                    {
                        const auto& r_variable =
                            KratosComponents<Variable<array_1d<double, 3>>>::Get(
                                r_variable_data.Name());
                        for (std::size_t i = 0; i < queue_size; ++i)
                            r_adjoint_node.FastGetSolutionStepValue(r_variable, i) =
                                r_primal_node.FastGetSolutionStepValue(r_variable, i);
                    }
                    else if (KratosComponents<Variable<double>>::Has(
                                 r_variable_data.Name()))
                    {
                        const auto& r_variable = KratosComponents<Variable<double>>::Get(
                            r_variable_data.Name());
                        for (std::size_t i = 0; i < queue_size; ++i)
                            r_adjoint_node.FastGetSolutionStepValue(r_variable, i) =
                                r_primal_node.FastGetSolutionStepValue(r_variable, i);
                    }
                    else
                        KRATOS_ERROR
                            << "Unsupported variable: " << r_variable_data.Name()
                            << std::endl;
                }
            }
            KRATOS_CATCH("");
        }

        const ModelPart& mrPrimalModelPart;
        ModelPart& mrAdjointModelPart;
        const std::string mElementType = "TotalLagrangianAdjointElement2D3N";
};

struct PrimalSolverFactory
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

    SolvingStrategyType::Pointer Execute(ModelPart& rModelPart)
    {
        LinearSolverType::Pointer p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        SchemeType::Pointer p_scheme =
            Kratos::make_shared<ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>>(0.0);
        ConvergenceCriteriaType::Pointer p_conv_criteria =
            Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(
                1e-24, 1e-25);
        return Kratos::make_shared<ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            rModelPart, p_scheme, p_linear_solver, p_conv_criteria, 30, true, false, true);
    }
};

AdjointResponseFunction::Pointer ResponseFunctionFactory(ModelPart& rModelPart)
{
    return Kratos::make_shared<NodeDisplacementResponseFunction>(Parameters(R"({"node_id": 2,
            "direction": [0.0, 1.0, 0.0]})"), rModelPart);
}

struct AdjointSolverFactory
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType>
    ConvergenceCriteriaType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType,
    LinearSolverType> SolvingStrategyType;

    static SolvingStrategyType::Pointer Execute(ModelPart& rModelPart,
                                                AdjointResponseFunction::Pointer pResponseFunction)
    {
        LinearSolverType::Pointer p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        Parameters settings(R"({ "alpha_bossak": 0.0 })");
        SchemeType::Pointer p_adjoint_scheme =
            Kratos::make_shared<ResidualBasedAdjointBossakScheme<SparseSpaceType, LocalSpaceType>>(
                settings, pResponseFunction);
        return Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            rModelPart, p_adjoint_scheme, p_linear_solver);
    }
};
}

double CalculateResponseValue(unsigned NodeToPerturb, char Direction, double Perturbation)
{
    Model this_model;
    ModelPart& primal_model_part = this_model.CreateModelPart("primal");
    PrimalModelPartFactory(primal_model_part).Execute();
    if (Direction == 'x')
        primal_model_part.GetNode(NodeToPerturb).X0() += Perturbation;
    else if (Direction == 'y')
        primal_model_part.GetNode(NodeToPerturb).Y0() += Perturbation;
    auto p_response_function = ResponseFunctionFactory(primal_model_part);
    auto p_solver = PrimalSolverFactory().Execute(primal_model_part);
    p_response_function->Initialize();
    p_solver->Initialize();
    double response_value = 0.;
    const double start_time = 0.0;
    const double delta_time = 0.008;
    const double end_time = 0.015;
    primal_model_part.CloneTimeStep(start_time - delta_time);
    primal_model_part.CloneTimeStep(start_time);
    for (double current_time = start_time;
         current_time < end_time;) // approx. 1/4 of a period.
    {
        current_time += delta_time;
        primal_model_part.CloneTimeStep(current_time);
        p_solver->Solve();
        response_value += delta_time * p_response_function->CalculateValue();
    }
    return response_value;
}

double CalculateSensitivity(unsigned NodeToPerturb, char Direction)
{
    Model this_model;
    ModelPart& primal_model_part = this_model.CreateModelPart("primal");
    PrimalModelPartFactory(primal_model_part).Execute();
    auto p_solver = PrimalSolverFactory().Execute(primal_model_part);
    p_solver->Initialize();
    const double start_time = 0.0;
    const double delta_time = 0.008;
    const double end_time = 0.015;
    primal_model_part.CloneTimeStep(start_time - delta_time);
    primal_model_part.CloneTimeStep(start_time);
    for (double current_time = start_time;
         current_time < end_time;) // approx. 1/4 of a period.
    {
        current_time += delta_time;
        primal_model_part.CloneTimeStep(current_time);
        p_solver->Solve();
    }

    ModelPart& adjoint_model_part = this_model.CreateModelPart("adjoint");
    AdjointModelPartFactory(primal_model_part, adjoint_model_part).Execute();
    auto p_response_function = ResponseFunctionFactory(adjoint_model_part);
    auto p_adjoint_solver =
        AdjointSolverFactory::Execute(adjoint_model_part, p_response_function);
    p_adjoint_solver->Initialize();
    SensitivityBuilder sensitivity_builder(
        Parameters(R"({"nodal_solution_step_sensitivity_variables": ["SHAPE_SENSITIVITY"]})"),
        adjoint_model_part, p_response_function);
    sensitivity_builder.Initialize();
    adjoint_model_part.CloneTimeStep(end_time + delta_time);
    adjoint_model_part.CloneTimeStep(end_time);
    for (auto& r_adjoint_node : adjoint_model_part.Nodes())
    {
        const auto& r_primal_node = primal_model_part.GetNode(r_adjoint_node.Id());
        r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT) = r_primal_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_adjoint_node.X() = r_adjoint_node.X0() + r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_X);
        r_adjoint_node.Y() = r_adjoint_node.Y0() + r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
        r_adjoint_node.FastGetSolutionStepValue(VELOCITY) = r_primal_node.FastGetSolutionStepValue(VELOCITY);
        r_adjoint_node.FastGetSolutionStepValue(ACCELERATION) = r_primal_node.FastGetSolutionStepValue(ACCELERATION);
    }
    p_adjoint_solver->Solve();
    sensitivity_builder.UpdateSensitivities();
    adjoint_model_part.CloneTimeStep(end_time - delta_time);
    for (auto& r_adjoint_node : adjoint_model_part.Nodes())
    {
        const auto& r_primal_node = primal_model_part.GetNode(r_adjoint_node.Id());
        r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT) = r_primal_node.FastGetSolutionStepValue(DISPLACEMENT,1);
        r_adjoint_node.X() = r_adjoint_node.X0() + r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_X);
        r_adjoint_node.Y() = r_adjoint_node.Y0() + r_adjoint_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
        r_adjoint_node.FastGetSolutionStepValue(VELOCITY) = r_primal_node.FastGetSolutionStepValue(VELOCITY,1);
        r_adjoint_node.FastGetSolutionStepValue(ACCELERATION) = r_primal_node.FastGetSolutionStepValue(ACCELERATION,1);
    }
    p_adjoint_solver->Solve();
    sensitivity_builder.UpdateSensitivities();
    auto& r_node = adjoint_model_part.GetNode(NodeToPerturb);
    double sensitivity{};
    if (Direction == 'x')
        sensitivity = r_node.FastGetSolutionStepValue(SHAPE_SENSITIVITY_X);
    else if (Direction == 'y')
        sensitivity = r_node.FastGetSolutionStepValue(SHAPE_SENSITIVITY_Y);
    else
        KRATOS_ERROR << "Invalid Direction: " << Direction << std::endl;
    return sensitivity;
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_SensitivityTransientOneElement, KratosStructuralMechanicsFastSuite)
{
    std::vector<double> fd_sensitivities, adjoint_sensitivities;
    const double delta = 1e-5;
    const double response_value0 = CalculateResponseValue(1, 'x', 0.);
    KRATOS_WATCH(response_value0);
    for (unsigned i_node : {1, 2, 3})
    {
        for (char d : {'x', 'y'})
        {
            const double response_value1 = CalculateResponseValue(i_node, d, delta);
            const double finite_diff_sensitivity =
                (response_value1 - response_value0) / delta;
            const double adjoint_sensitivity = CalculateSensitivity(i_node, d);
            fd_sensitivities.push_back(finite_diff_sensitivity);
            adjoint_sensitivities.push_back(adjoint_sensitivity);
        }
    }
    for (std::size_t i = 0; i < fd_sensitivities.size(); ++i)
    {
        KRATOS_CHECK_NEAR(adjoint_sensitivities.at(i), fd_sensitivities.at(i), 1e-9);
    }
}
}
}
