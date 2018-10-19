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
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
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
        (*p_prop)[RAYLEIGH_ALPHA] = 0.02;
        (*p_prop)[RAYLEIGH_BETA] = 0.03;
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
            Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>>();
        ConvergenceCriteriaType::Pointer p_conv_criteria =
            Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(
                1e-12, 1e-14);
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
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

    SolvingStrategyType::Pointer Execute(ModelPart& rModelPart, AdjointResponseFunction::Pointer pResponseFunction)
    {
        LinearSolverType::Pointer p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        SchemeType::Pointer p_adjoint_scheme =
            Kratos::make_shared<ResidualBasedAdjointStaticScheme<SparseSpaceType, LocalSpaceType>>(
                pResponseFunction);
        return Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            rModelPart, p_adjoint_scheme, p_linear_solver);
    }
};
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_SensitivityOneElement, KratosStructuralMechanicsFastSuite)
{
    // Compute primal solution.
    Model this_model;
    ModelPart& primal_model_part = this_model.CreateModelPart("primal");
    PrimalModelPartFactory(primal_model_part).Execute();
    auto p_solver = PrimalSolverFactory().Execute(primal_model_part);
    p_solver->Initialize();
    p_solver->Solve();
    // Compute adjoint solution.
    ModelPart& adjoint_model_part = this_model.CreateModelPart("adjoint");
    AdjointModelPartFactory(primal_model_part, adjoint_model_part).Execute();
    auto p_response_function = ResponseFunctionFactory(primal_model_part);
    auto p_adjoint_response_function = ResponseFunctionFactory(adjoint_model_part);
    auto p_adjoint_solver =
        AdjointSolverFactory().Execute(adjoint_model_part, p_adjoint_response_function);
    p_adjoint_solver->Initialize();
    p_adjoint_solver->Solve();
    SensitivityBuilder sensitivity_builder(
        Parameters(R"({"nodal_solution_step_sensitivity_variables": ["SHAPE_SENSITIVITY"], "build_mode": "static"})"),
        adjoint_model_part, p_adjoint_response_function);
    sensitivity_builder.Initialize();
    sensitivity_builder.UpdateSensitivities();
    // Compare with finite difference sensitivity.
    const double delta = 1e-7;
    const double response_value0 = p_response_function->CalculateValue();
    for (auto& r_node : primal_model_part.Nodes())
    {
        for (std::size_t d = 0; d < 2; ++d)
        {
            r_node.GetInitialPosition()[d] += delta;
            p_solver->Solve();
            const double response_value1 = p_response_function->CalculateValue();
            const double finite_diff_sensitivity =
                (response_value1 - response_value0) / delta;
            const double adjoint_sensitivity =
                adjoint_model_part.GetNode(r_node.Id()).FastGetSolutionStepValue(SHAPE_SENSITIVITY)[d];
            KRATOS_CHECK_NEAR(finite_diff_sensitivity, adjoint_sensitivity, 1e-8);
            r_node.GetInitialPosition()[d] -= delta;
        }
    }
}
}
}
